"""
Render docs/DNA_stress_repair_panel.md to a typeset PDF.

Uses DejaVu Sans so that the Greek letters, arrows and typographic dashes in the
report survive into the PDF - the Helvetica standard encoding does not carry them
and would silently emit black boxes.

Usage:  python3 src/build_report_pdf.py [input.md] [output.pdf]
"""

import os
import re
import sys

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, Frame, KeepTogether, PageBreak,
                               PageTemplate, Paragraph, Spacer, Table, TableStyle)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# Liberation Sans supplies a complete four-face family (DejaVu Sans ships no oblique
# in this image), and DejaVu Sans Mono covers the code spans. Both were checked to
# carry every Greek letter, arrow, superscript and dash the report uses.
FONTS = {
    "Body":     ("Body",   "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"),
    "BodyBd":   ("Body-Bd", "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"),
    "BodyIt":   ("Body-It", "/usr/share/fonts/truetype/liberation/LiberationSans-Italic.ttf"),
    "BodyBdIt": ("Body-BdIt", "/usr/share/fonts/truetype/liberation/LiberationSans-BoldItalic.ttf"),
    "Mono":     ("Mono",   "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf"),
    "MonoBd":   ("Mono-Bd", "/usr/share/fonts/truetype/dejavu/DejaVuSansMono-Bold.ttf"),
}

INK      = colors.HexColor("#1a1a1a")
MUTED    = colors.HexColor("#5b6470")
ACCENT   = colors.HexColor("#0b5d8f")
RULE     = colors.HexColor("#c9d1d9")
HEADBG   = colors.HexColor("#eef2f6")
CODEBG   = colors.HexColor("#f5f7f9")
QUOTEBAR = colors.HexColor("#0b5d8f")
QUOTEBG  = colors.HexColor("#f0f6fb")


def register_fonts():
    for _, (name, path) in FONTS.items():
        if not os.path.exists(path):
            raise SystemExit("missing font: %s" % path)
        pdfmetrics.registerFont(TTFont(name, path))
    pdfmetrics.registerFontFamily("Body", normal="Body", bold="Body-Bd",
                                  italic="Body-It", boldItalic="Body-BdIt")
    pdfmetrics.registerFontFamily("Mono", normal="Mono", bold="Mono-Bd",
                                  italic="Mono", boldItalic="Mono-Bd")


# ----------------------------------------------------------------------------------
# inline markdown -> reportlab inline markup
# ----------------------------------------------------------------------------------
def esc(text):
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def inline(text):
    """Convert inline markdown to reportlab's mini-HTML, escaping first.

    Code spans are extracted before escaping so that their contents are never
    interpreted as markdown emphasis - a span like `**kwargs` must stay literal.
    """
    spans = []

    def stash(m):
        spans.append(m.group(1))
        return "\x00%d\x00" % (len(spans) - 1)

    text = re.sub(r"`([^`]+)`", stash, text)
    text = esc(text)

    # Links before emphasis, so bracketed link text is not mangled. The URL pattern
    # allows one level of BALANCED parentheses: DOIs routinely contain them
    # (10.1016/s0300-483x(02)00480-8), and a naive [^)]+ truncates the URL at the
    # first inner ')' and spills the remainder into the body text.
    text = re.sub(r"\[([^\]]+)\]\(((?:[^()\s]|\([^()]*\))+)\)",
                  lambda m: '<link href="%s" color="#0b5d8f">%s</link>'
                            % (m.group(2), m.group(1)), text)
    text = re.sub(r"\*\*\*(.+?)\*\*\*", r"<b><i>\1</i></b>", text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*([^*]+?)\*(?!\*)", r"<i>\1</i>", text)

    def unstash(m):
        body = esc(spans[int(m.group(1))])
        return ('<font face="Mono" size="8.4" backColor="#f0f2f4">%s</font>' % body)

    return re.sub(r"\x00(\d+)\x00", unstash, text)


# ----------------------------------------------------------------------------------
# styles
# ----------------------------------------------------------------------------------
def build_styles():
    ss = getSampleStyleSheet()
    S = {}
    S["title"] = ParagraphStyle("title", parent=ss["Normal"], fontName="Body-Bd",
                                fontSize=21, leading=26, textColor=INK, spaceAfter=4)
    S["subtitle"] = ParagraphStyle("subtitle", parent=ss["Normal"], fontName="Body-It",
                                   fontSize=11, leading=15, textColor=MUTED, spaceAfter=16)
    S["h1"] = ParagraphStyle("h1", parent=ss["Normal"], fontName="Body-Bd",
                             fontSize=15.5, leading=19.5, textColor=ACCENT,
                             spaceBefore=17, spaceAfter=7)
    S["h2"] = ParagraphStyle("h2", parent=ss["Normal"], fontName="Body-Bd",
                             fontSize=11.8, leading=15, textColor=INK,
                             spaceBefore=12, spaceAfter=5)
    S["h3"] = ParagraphStyle("h3", parent=ss["Normal"], fontName="Body-Bd",
                             fontSize=10.2, leading=13, textColor=MUTED,
                             spaceBefore=9, spaceAfter=4)
    S["body"] = ParagraphStyle("body", parent=ss["Normal"], fontName="Body",
                               fontSize=9.3, leading=13.6, textColor=INK,
                               alignment=TA_LEFT, spaceAfter=6.5)
    S["bullet"] = ParagraphStyle("bullet", parent=S["body"], leftIndent=11,
                                 bulletIndent=2, spaceAfter=3.6)
    S["quote"] = ParagraphStyle("quote", parent=S["body"], leftIndent=8, rightIndent=6,
                                fontSize=9.1, leading=13.4, spaceAfter=0, spaceBefore=0)
    S["code"] = ParagraphStyle("code", parent=ss["Normal"], fontName="Mono",
                               fontSize=8.1, leading=11.2, textColor=INK, spaceAfter=0)
    S["th"] = ParagraphStyle("th", parent=ss["Normal"], fontName="Body-Bd",
                             fontSize=8.1, leading=10.4, textColor=INK)
    S["td"] = ParagraphStyle("td", parent=ss["Normal"], fontName="Body",
                             fontSize=8.1, leading=10.4, textColor=INK)
    return S


# ----------------------------------------------------------------------------------
# markdown parsing
# ----------------------------------------------------------------------------------
def split_row(line):
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    # split on | not inside a code span
    parts, buf, in_code = [], "", False
    for ch in line:
        if ch == "`":
            in_code = not in_code
        if ch == "|" and not in_code:
            parts.append(buf); buf = ""
        else:
            buf += ch
    parts.append(buf)
    return [p.strip() for p in parts]


def is_sep_row(line):
    return bool(re.fullmatch(r"\|?[\s:|-]+\|?", line.strip())) and "-" in line


def parse(md):
    """Return a list of block tuples."""
    lines = md.split("\n")
    blocks, i = [], 0
    while i < len(lines):
        raw = lines[i]
        line = raw.rstrip()
        stripped = line.strip()

        if not stripped:
            i += 1; continue

        if stripped == "---":
            blocks.append(("rule", None)); i += 1; continue

        if stripped.startswith("```"):
            i += 1
            buf = []
            while i < len(lines) and not lines[i].strip().startswith("```"):
                buf.append(lines[i]); i += 1
            i += 1
            blocks.append(("code", buf)); continue

        m = re.match(r"^(#{1,6})\s+(.*)$", stripped)
        if m:
            blocks.append(("h%d" % len(m.group(1)), m.group(2))); i += 1; continue

        # table
        if stripped.startswith("|") and i + 1 < len(lines) and is_sep_row(lines[i + 1]):
            header = split_row(stripped)
            i += 2
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                rows.append(split_row(lines[i])); i += 1
            blocks.append(("table", (header, rows))); continue

        # blockquote (may contain its own bullets / paragraphs)
        if stripped.startswith(">"):
            buf = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                buf.append(re.sub(r"^\s*>\s?", "", lines[i])); i += 1
            blocks.append(("quote", buf)); continue

        # bullet list
        if re.match(r"^\s*[-*]\s+", raw):
            items = []
            while i < len(lines):
                mm = re.match(r"^(\s*)[-*]\s+(.*)$", lines[i])
                if mm:
                    items.append(mm.group(2).strip()); i += 1
                elif lines[i].strip() and lines[i].startswith(("  ", "\t")) and items:
                    items[-1] += " " + lines[i].strip(); i += 1
                else:
                    break
            blocks.append(("ul", items)); continue

        # ordered list
        if re.match(r"^\s*\d+[.)]\s+", raw):
            items = []
            while i < len(lines):
                mm = re.match(r"^\s*(\d+)[.)]\s+(.*)$", lines[i])
                if mm:
                    items.append((mm.group(1), mm.group(2).strip())); i += 1
                elif lines[i].strip() and lines[i].startswith(("  ", "\t")) and items:
                    items[-1] = (items[-1][0], items[-1][1] + " " + lines[i].strip()); i += 1
                else:
                    break
            blocks.append(("ol", items)); continue

        # paragraph
        buf = []
        while i < len(lines) and lines[i].strip() and not re.match(
                r"^\s*([-*]\s|\d+[.)]\s|#{1,6}\s|>|\||```)", lines[i]) \
                and lines[i].strip() != "---":
            buf.append(lines[i].strip()); i += 1
        if buf:
            blocks.append(("p", " ".join(buf)))
        else:
            i += 1
    return blocks


# ----------------------------------------------------------------------------------
# flowable construction
# ----------------------------------------------------------------------------------
def hrule(width):
    t = Table([[""]], colWidths=[width], rowHeights=[0.7])
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), RULE),
                           ("TOPPADDING", (0, 0), (-1, -1), 0),
                           ("BOTTOMPADDING", (0, 0), (-1, -1), 0)]))
    return t


def code_block(buf, S, width):
    body = "<br/>".join(esc(l).replace(" ", "&nbsp;") for l in buf) or "&nbsp;"
    inner = Paragraph(body, S["code"])
    t = Table([[inner]], colWidths=[width])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), CODEBG),
        ("BOX", (0, 0), (-1, -1), 0.5, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 7), ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 6), ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
    ]))
    return t


def quote_block(buf, S, width):
    inner = []
    for blk in parse("\n".join(buf)):
        inner.extend(flow(blk, S, width - 16))
    t = Table([[inner]], colWidths=[width])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), QUOTEBG),
        ("LINEBEFORE", (0, 0), (0, -1), 2.4, QUOTEBAR),
        ("LEFTPADDING", (0, 0), (-1, -1), 9), ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 7), ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    return t


MD_SYNTAX = re.compile(r"(\*\*|\*|`)")


def plain(cell):
    """Cell text with markdown syntax and link targets removed, for measurement."""
    t = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", cell)
    return MD_SYNTAX.sub("", t)


def allocate_columns(header, rows, width, font, size, pad, head_font=None, head_size=None):
    """Give every column at least the width of its longest unbreakable word.

    Measuring character counts is not enough: a column whose longest entry is a
    single long token cannot wrap, so an under-allocated column silently CLIPS that
    token rather than reflowing it. Minimum widths therefore come from real string
    metrics, and only the leftover space is distributed by preference.
    """
    # Header cells are set in the bold header style, so they must be measured in that
    # font - measuring a bold header with the regular metrics under-allocates and the
    # header word wraps mid-letter.
    head_font = head_font or font
    head_size = head_size or size
    ncol = len(header)
    mins, prefs = [], []
    for c in range(ncol):
        cells = [(plain(header[c]), head_font, head_size)]
        cells += [(plain(r[c]), font, size) for r in rows]
        longest_word = 0.0
        longest_cell = 0.0
        for text, fnt, sz in cells:
            for word in text.split():
                longest_word = max(longest_word, pdfmetrics.stringWidth(word, fnt, sz))
            longest_cell = max(longest_cell, pdfmetrics.stringWidth(text, fnt, sz))
        mins.append(longest_word + 2 * pad + 1.0)
        prefs.append(min(longest_cell + 2 * pad, width * 0.62))

    total_min = sum(mins)
    if total_min > width:
        # pathological table: scale minimums down proportionally and accept wrapping
        k = width / total_min
        return [m * k for m in mins]

    slack = width - total_min
    head = [max(0.0, p - m) for p, m in zip(prefs, mins)]
    hs = sum(head)
    if hs <= 0:
        return [m + slack / ncol for m in mins]
    return [m + slack * h / hs for m, h in zip(mins, head)]


def md_table(header, rows, S, width):
    ncol = max(len(header), max((len(r) for r in rows), default=0))
    header = header + [""] * (ncol - len(header))
    norm = [r + [""] * (ncol - len(r)) for r in rows]

    pad = 4.5
    widths = allocate_columns(header, norm, width,
                              S["td"].fontName, S["td"].fontSize, pad,
                              S["th"].fontName, S["th"].fontSize)

    data = [[Paragraph(inline(h), S["th"]) for h in header]]
    for r in norm:
        data.append([Paragraph(inline(c), S["td"]) for c in r])

    t = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    style = [
        ("BACKGROUND", (0, 0), (-1, 0), HEADBG),
        ("LINEBELOW", (0, 0), (-1, 0), 0.8, RULE),
        ("GRID", (0, 0), (-1, -1), 0.3, RULE),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), pad),
        ("RIGHTPADDING", (0, 0), (-1, -1), pad),
        ("TOPPADDING", (0, 0), (-1, -1), 3.6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3.6),
    ]
    for ri in range(1, len(data)):
        if ri % 2 == 0:
            style.append(("BACKGROUND", (0, ri), (-1, ri), colors.HexColor("#fafbfc")))
    t.setStyle(TableStyle(style))
    return t


def flow(block, S, width):
    kind, payload = block
    if kind == "rule":
        return [Spacer(1, 5), hrule(width), Spacer(1, 7)]
    if kind == "code":
        return [code_block(payload, S, width), Spacer(1, 7)]
    if kind == "quote":
        return [Spacer(1, 2), quote_block(payload, S, width), Spacer(1, 8)]
    if kind == "table":
        header, rows = payload
        return [Spacer(1, 3), md_table(header, rows, S, width), Spacer(1, 9)]
    if kind == "h1":
        return []          # document title handled separately
    if kind == "h2":
        return [Paragraph(inline(payload), S["h1"])]
    if kind == "h3":
        return [Paragraph(inline(payload), S["h2"])]
    if kind in ("h4", "h5", "h6"):
        return [Paragraph(inline(payload), S["h3"])]
    if kind == "ul":
        return [Paragraph(inline(t), S["bullet"], bulletText="•") for t in payload]
    if kind == "ol":
        return [Paragraph(inline(t), S["bullet"], bulletText="%s." % n) for n, t in payload]
    if kind == "p":
        return [Paragraph(inline(payload), S["body"])]
    return []


# ----------------------------------------------------------------------------------
# document
# ----------------------------------------------------------------------------------
class Doc(BaseDocTemplate):
    def __init__(self, path, **kw):
        BaseDocTemplate.__init__(self, path, **kw)
        frame = Frame(self.leftMargin, self.bottomMargin,
                      self.width, self.height, id="main",
                      leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
        self.addPageTemplates([PageTemplate(id="std", frames=[frame],
                                            onPage=self.decorate)])
        self.running_title = ""

    def decorate(self, canv, doc):
        canv.saveState()
        canv.setFont("Body", 7.4)
        canv.setFillColor(MUTED)
        if doc.page > 1:
            canv.drawString(doc.leftMargin, A4[1] - 13 * mm, self.running_title)
            canv.setStrokeColor(RULE); canv.setLineWidth(0.4)
            canv.line(doc.leftMargin, A4[1] - 15.4 * mm,
                      A4[0] - doc.rightMargin, A4[1] - 15.4 * mm)
        canv.drawRightString(A4[0] - doc.rightMargin, 11 * mm, "%d" % doc.page)
        canv.drawString(doc.leftMargin, 11 * mm,
                        "DNA stress and repair capacity panel")
        canv.restoreState()


def main():
    md_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        ROOT, "docs", "DNA_stress_repair_panel.md")
    out_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        ROOT, "docs", "DNA_stress_repair_panel.pdf")

    register_fonts()
    S = build_styles()
    md = open(md_path, encoding="utf-8").read()
    blocks = parse(md)

    doc = Doc(out_path, pagesize=A4,
              leftMargin=19 * mm, rightMargin=17 * mm,
              topMargin=19 * mm, bottomMargin=17 * mm,
              title="A multiplex panel for DNA stress and DNA repair capacity",
              author="ozangenomics", subject="DNA damage and repair biomarker panel design")
    doc.running_title = "A multiplex panel for DNA stress and DNA repair capacity"
    width = doc.width

    story = []
    # title page block
    h1 = next((p for k, p in blocks if k == "h1"), "Report")
    story.append(Paragraph(inline(h1), S["title"]))
    story.append(hrule(width))
    story.append(Spacer(1, 6))
    # the bold lead paragraph directly after the title becomes the standfirst
    first_p = next((p for k, p in blocks if k == "p"), None)
    if first_p:
        story.append(Paragraph(inline(first_p.strip("*")), S["subtitle"]))

    skipped_lead = False
    for idx, blk in enumerate(blocks):
        if blk[0] == "h1":
            continue
        if blk[0] == "p" and not skipped_lead and blk[1] == first_p:
            skipped_lead = True
            continue
        # A rule immediately before a section heading adds nothing - the heading already
        # separates the sections - and if it lands at a page top it reads as a stray
        # duplicate of the running-header rule.
        if blk[0] == "rule":
            nxt = blocks[idx + 1][0] if idx + 1 < len(blocks) else None
            if nxt == "h2":
                story.append(Spacer(1, 10))
                continue
        story.extend(flow(blk, S, width))

    doc.build(story)
    size = os.path.getsize(out_path)
    print("wrote %s (%.1f KB)" % (out_path, size / 1024.0))
    print("blocks rendered: %d" % len(blocks))
    from collections import Counter
    print("block mix:", dict(Counter(b[0] for b in blocks)))


if __name__ == "__main__":
    main()
