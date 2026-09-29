# Exhaled breath condensate (EBC) collection protocol

Step-by-step, one-page (A4) illustrated protocol for nurses and physicians collecting
exhaled breath condensate for the breath phenotyping study.

| File | Use |
|------|-----|
| `EBC_collection_protocol_A4.pdf` | Print-ready A4 page (vector, fonts embedded) |
| `EBC_collection_protocol_A4.png` | Same page as an image (1588 × 2246 px, 2× resolution) for slides, e-mail or a wall poster |
| `EBC_collection_protocol_A4.html` | Editable source; open in any browser, print with margins set to none |
| `fonts/` | Manrope and Source Sans 3 (Google Fonts, OFL) used by the HTML page |

Fill in the placeholders in the footer before distribution: `[DATE]` and `[STUDY COORDINATOR]`.

## The five steps

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

## Critical points

- **Cold chain:** condenser pre-cooled, tube frozen at once, dry ice in transit.
- **Clean sample:** new mouthpiece each time; no saliva, no coughing.
- **Traceability:** label and log every tube: ID, date, time, volume.

## Regenerating the PDF and PNG

```bash
chromium --headless=new --no-pdf-header-footer \
  --print-to-pdf=EBC_collection_protocol_A4.pdf EBC_collection_protocol_A4.html
chromium --headless=new --force-device-scale-factor=2 --window-size=794,1400 \
  --screenshot=page.png EBC_collection_protocol_A4.html   # then crop to 1588 × 2246
```
