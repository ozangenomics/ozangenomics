// Build the PCT draft .docx from the Markdown source using docx (npm).
// Usage: node build_docx.js <input.md> <output.docx>
// Supported Markdown subset: #/##/###/#### headings, paragraphs, pipe tables,
// "- " bullets, numbered claims ("1. ..." with indented sub-lines), images
// ![caption](path), **bold**, *italic*, and <<<PAGEBREAK>>> markers.
const fs = require("fs");
const path = require("path");
const docx = require("docx");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell,
  WidthType, AlignmentType, ImageRun, PageBreak, BorderStyle, ShadingType,
  LevelFormat, Footer, PageNumber, TabStopType,
} = docx;

const [,, inPath, outPath] = process.argv;
const srcDir = path.dirname(inPath);
const md = fs.readFileSync(inPath, "utf8").split(/\r?\n/);

const FONT = "Times New Roman";
const BODY_SIZE = 24; // 12 pt

function inlineRuns(text, base = {}) {
  // split on **bold** and *italic*
  const runs = [];
  const re = /(\*\*[^*]+\*\*|\*[^*]+\*)/g;
  let last = 0, m;
  while ((m = re.exec(text)) !== null) {
    if (m.index > last) runs.push(new TextRun({ text: text.slice(last, m.index), font: FONT, size: BODY_SIZE, ...base }));
    const tok = m[0];
    if (tok.startsWith("**")) runs.push(new TextRun({ text: tok.slice(2, -2), bold: true, font: FONT, size: BODY_SIZE, ...base }));
    else runs.push(new TextRun({ text: tok.slice(1, -1), italics: true, font: FONT, size: BODY_SIZE, ...base }));
    last = m.index + tok.length;
  }
  if (last < text.length) runs.push(new TextRun({ text: text.slice(last), font: FONT, size: BODY_SIZE, ...base }));
  return runs;
}

function para(text, opts = {}) {
  return new Paragraph({ children: inlineRuns(text), spacing: { after: 160, line: 300 }, alignment: AlignmentType.JUSTIFIED, ...opts });
}

function heading(text, level) {
  const map = { 1: HeadingLevel.HEADING_1, 2: HeadingLevel.HEADING_2, 3: HeadingLevel.HEADING_3, 4: HeadingLevel.HEADING_4 };
  return new Paragraph({ text, heading: map[level], spacing: { before: 280, after: 160 } });
}

function cellBorders() {
  const b = { style: BorderStyle.SINGLE, size: 4, color: "000000" };
  return { top: b, bottom: b, left: b, right: b };
}

function buildTable(rows) {
  const ncol = rows[0].length;
  const total = 9360; // 6.5in text width in DXA
  const widths = Array(ncol).fill(Math.floor(total / ncol));
  widths[ncol - 1] = total - widths.slice(0, -1).reduce((a, b) => a + b, 0);
  const trows = rows.map((r, ri) => new TableRow({
    tableHeader: ri === 0,
    children: r.map((c, ci) => new TableCell({
      width: { size: widths[ci], type: WidthType.DXA },
      borders: cellBorders(),
      shading: ri === 0 ? { type: ShadingType.CLEAR, fill: "E7E6E6", color: "auto" } : undefined,
      margins: { top: 60, bottom: 60, left: 80, right: 80 },
      children: [new Paragraph({ children: inlineRuns(c, { size: 18, bold: ri === 0 }), spacing: { after: 0 } })],
    })),
  }));
  return new Table({ rows: trows, columnWidths: widths, width: { size: total, type: WidthType.DXA } });
}

function image(caption, rel) {
  const p = path.resolve(srcDir, rel);
  const data = fs.readFileSync(p);
  // read PNG dimensions
  const w = data.readUInt32BE(16), h = data.readUInt32BE(20);
  const maxW = 620, maxH = 760;
  let dw = maxW, dh = Math.round(h * maxW / w);
  if (dh > maxH) { dh = maxH; dw = Math.round(w * maxH / h); }
  return [
    new Paragraph({ alignment: AlignmentType.CENTER, children: [new ImageRun({ type: "png", data, transformation: { width: dw, height: dh } })], spacing: { before: 200, after: 80 } }),
    new Paragraph({ alignment: AlignmentType.CENTER, children: inlineRuns(caption, { size: 20, italics: true }), spacing: { after: 240 } }),
  ];
}

const children = [];
let i = 0;
let tableBuf = [];
function flushTable() {
  if (tableBuf.length) {
    const rows = tableBuf.filter(l => !/^\|\s*:?-{2,}/.test(l)).map(l => l.replace(/^\||\|$/g, "").split("|").map(s => s.trim()));
    children.push(buildTable(rows));
    children.push(new Paragraph({ spacing: { after: 120 } }));
    tableBuf = [];
  }
}

while (i < md.length) {
  const line = md[i];
  if (/^\|/.test(line)) { tableBuf.push(line); i++; continue; }
  flushTable();
  if (line.trim() === "") { i++; continue; }
  if (line.trim() === "<<<PAGEBREAK>>>") { children.push(new Paragraph({ children: [new PageBreak()] })); i++; continue; }
  let m;
  if ((m = line.match(/^(#{1,4})\s+(.*)$/))) { children.push(heading(m[2], m[1].length)); i++; continue; }
  if ((m = line.match(/^!\[(.*?)\]\((.*?)\)$/))) { children.push(...image(m[1], m[2])); i++; continue; }
  if ((m = line.match(/^-\s+(.*)$/))) {
    children.push(new Paragraph({ children: inlineRuns(m[1]), numbering: { reference: "bullets", level: 0 }, spacing: { after: 80 } }));
    i++; continue;
  }
  if ((m = line.match(/^(\d+)\.\s+(.*)$/))) {
    // claim: gather indented continuation lines
    const num = m[1]; let text = m[2]; const subs = [];
    i++;
    while (i < md.length && /^\s{3,}\S/.test(md[i])) { subs.push(md[i].trim()); i++; }
    children.push(new Paragraph({ children: [new TextRun({ text: `${num}. `, bold: true, font: FONT, size: BODY_SIZE }), ...inlineRuns(text)], spacing: { before: 160, after: subs.length ? 40 : 160 }, alignment: AlignmentType.JUSTIFIED, keepNext: subs.length > 0 }));
    subs.forEach(s => children.push(new Paragraph({ children: inlineRuns(s), indent: { left: 720 }, spacing: { after: 40 }, alignment: AlignmentType.JUSTIFIED })));
    continue;
  }
  // paragraph: join wrapped lines until blank
  let text = line.trim(); i++;
  while (i < md.length && md[i].trim() !== "" && !/^(#|\||!\[|-\s|\d+\.\s|<<<)/.test(md[i])) { text += " " + md[i].trim(); i++; }
  children.push(para(text));
}
flushTable();

const doc = new Document({
  creator: "Dalipharma",
  title: "PCT Application Draft - Ozonated Olive Oil Nanoemulsion",
  styles: {
    default: { document: { run: { font: FONT, size: BODY_SIZE } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 30, bold: true, font: FONT }, paragraph: { spacing: { before: 360, after: 200 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 26, bold: true, font: FONT, allCaps: true }, paragraph: { spacing: { before: 320, after: 160 }, outlineLevel: 1 } },
      { id: "Heading3", name: "Heading 3", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 24, bold: true, italics: true, font: FONT }, paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 2 } },
      { id: "Heading4", name: "Heading 4", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 24, bold: true, font: FONT }, paragraph: { spacing: { before: 200, after: 100 }, outlineLevel: 3 } },
    ],
  },
  numbering: { config: [{ reference: "bullets", levels: [{ level: 0, format: LevelFormat.BULLET, text: "•", alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 720, hanging: 360 } } } }] }] },
  sections: [{
    properties: { page: { size: { width: 11906, height: 16838 }, margin: { top: 1418, bottom: 1134, left: 1418, right: 1134 } } }, // A4, PCT Rule 11 margins (≥2.5 cm left, ≥2 cm others)
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ children: ["- ", PageNumber.CURRENT, " -"], font: FONT, size: 20 })] })] }) },
    children,
  }],
});

Packer.toBuffer(doc).then(buf => { fs.writeFileSync(outPath, buf); console.log("wrote", outPath, buf.length, "bytes"); });
