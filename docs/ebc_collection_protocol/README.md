# Exhaled breath condensate (EBC) collection protocol · Lanserhof Lans

Two-page A4 document for the breath phenotyping study, styled in Bentley green, dark brown and
champagne gold on warm cream, with a serif display face.

- **Page 1 · Collection protocol** for nurses and physicians: the five steps from the pre-cooled
  condenser to the analysis laboratory, critical points, and a sample-form block.
- **Page 2 · Participant preparation**: countdown timeline, "please avoid" and "please do" lists,
  what to tell the team, and what happens during sampling.

| File | Use |
|------|-----|
| `EBC_collection_protocol_A4.pdf` | Print-ready two-page A4 PDF, fonts embedded |
| `EBC_collection_protocol_A4_page1.png`, `..._page2.png` | The pages as images (1588 × 2246 px, 2× resolution) |
| `EBC_collection_protocol_A4.html` | Editable source for both pages; print with margins set to none |
| `fonts/` | Cormorant Garamond and Source Sans 3 (Google Fonts, OFL) |

The footer of page 1 has blank lines for the study coordinator and the date to be filled in by hand.
The address used is Lanserhof Lans, Kochholzweg 153, 6072 Lans, Austria, +43 512 386660.
Colour tokens (change in the HTML or the canvas Tweaks panel): Bentley green `#0F3B2C` (headings,
illustration accents), dark brown `#2A1F17` (body ink) and espresso `#2E211B` (dark bands, step
numbers), champagne gold `#C9B37E` (hairlines, rings) and dark gold `#8C7449` (small text), cream
`#F4EEE1`, sand `#EFE6D2`, paper `#FBF9F4`, burgundy for warnings `#7A3B2E`.

## Page 1: the five steps

1. **Retrieve the pre-cooled condenser.** Take the −20 °C breath condenser from the freezer in the
   sample collection room. Check that the Eppendorf tube is seated at the bottom of the condensation
   tube, fit a clean mouthpiece and start at once while the tube is cold. Gloves on; label the tube
   with participant ID, date and time.
2. **Exhale strongly for 2 minutes.** The participant breathes in through the nose and exhales strongly
   and steadily through the mouthpiece into the −20 °C tube for 2 timed minutes. Keep the device
   upright so the condensate runs down into the Eppendorf tube. No coughing, talking or saliva into
   the mouthpiece.
3. **Recover the condensed sample.** Detach the Eppendorf tube from the bottom of the condenser.
   Expected yield is about 200–400 µL of condensate, which may be partly frozen. Cap immediately and
   record the approximate volume on the sample form. Low volume is recorded, never topped up or pooled.
4. **Store the Eppendorf tube at −20 °C.** Place the capped, labelled tube upright in the cryobox in the
   −20 °C freezer immediately after collection and log its position and time. Return the condenser to
   the freezer to re-cool for the next participant. Once frozen, it stays frozen: no thaw and refreeze.
5. **Ship to Gipfel Life Sciences on dry ice.** Send frozen samples to the Gipfel Life Sciences
   laboratory in an insulated box packed with dry ice so they stay at −20 °C or colder throughout.
   Enclose the sample manifest (ID, date, volume) and notify the lab before dispatch. Samples must
   never thaw in transit.

Critical points: cold chain unbroken; clean sample (new mouthpiece, no saliva or coughing);
traceability (label and log every tube).

## Page 2: participant preparation

Countdown: 2 hours before, last food · 1 hour before, last drink other than water · sample time,
rinse mouth with water, sit down, breathe out for 2 minutes.

Please avoid: eating for 2 hours (including gum, sweets, lozenges); any drink except water for at
least 1 hour; smoking, vaping or nicotine for at least 2 hours; mouthwash, throat sprays or lozenges
within 1 hour; sexual intercourse for at least 2 hours; strenuous exercise within 1 hour.

Please do: drink still water as usual; rinse the mouth with water just before sampling; arrive early
and rest seated; breathe in through the nose and out through the mouthpiece for the full 2 minutes;
take usual medication unless told otherwise and report it.

Tell the team before sampling about a respiratory infection in the last 2 weeks, asthma, COPD or
inhaler use, smoking or vaping, or any rule that could not be followed.

## Regenerating the PDF and PNGs

```bash
chromium --headless=new --no-pdf-header-footer \
  --print-to-pdf=EBC_collection_protocol_A4.pdf EBC_collection_protocol_A4.html
chromium --headless=new --force-device-scale-factor=2 --window-size=794,2500 \
  --screenshot=pages.png EBC_collection_protocol_A4.html   # crop 1588 × 2246 per page
```
