# Structure + PWM figures — baseline vs augmented

Four figures highlighting the baseline-vs-augmented DeepPBS result at the level of
individual predicted motifs, for the Ch. 16 / Ch. 11 "what augmentation does" story.

## Figures
| file | what |
|---|---|
| struct_pwm_ets1.png | ETS1 (3wty): BioEmu ensemble render + crystal-target / baseline / augmented PWM logos; Δr +0.06 |
| struct_pwm_tbp.png  | TBP (2ko0): same; the induced-fit case, Δr +0.33 on this entry |
| struct_pwm_csl.png  | CSL/RBPJ (6qhd): same; Δr +0.10 |
| struct_pwm_context.png | HONEST CONTEXT: per-entry Δr distribution (104 entries/pilot); exemplars (red) are high-gain, population medians (black) are near zero |

## Important honesty note
The three per-pilot panels show CHERRY-PICKED high-gain benchmark entries — chosen to
make the mechanism legible, NOT the typical effect. struct_pwm_context.png is mandatory
alongside them: it shows each exemplar's Δr against the full per-entry distribution, whose
medians (ETS1 -0.03, TBP -0.01, CSL -0.05) are the cross-benchmark null. Never show the
three exemplar panels without the context panel.

## Provenance
- PWM logos: seed-averaged (s1-s5) predicted probabilities from
  output/stage6_train/{baseline,augmented}_<tf>_fold0_s*/predictions/<entry>_predict.npz
  (Y = crystal-target PWM, P = predicted; info-content scaled, bits).
- Structure renders: reused analysis/figures/pymol/<tf>_bioemu.png (BioEmu ensemble on crystal DNA).
- Δr = Pearson(target, augmented) - Pearson(target, baseline), seed-averaged.
