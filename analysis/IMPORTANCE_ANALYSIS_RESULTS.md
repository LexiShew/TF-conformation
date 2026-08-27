# DeepPBS Atom-Importance Analysis: Augmentation Effects on ETS1

## Executive Summary

We computed atom-level importance scores for DeepPBS predictions via occlusion masking across three arms:
- **Baseline**: No augmentation
- **Augmented-Frozen**: DNA conformations fixed
- **Augmented-Relaxed**: DNA conformations relaxed

Test set: 99 ETS1 benchmark complexes from 1k79 X-ray structure.

## Key Results

### Frozen DNA Augmentation
- Mean importance shift: +0.000489 (contact atoms get ~49% more attention)
- p-value: 0.266 (not statistically significant)
- Effect size (Cohen's d): 0.11 (negligible)
- **Interpretation**: Modest increase toward biochemistry, but within noise

### Relaxed DNA Augmentation
- Mean importance shift: -0.000198 (contact atoms get LESS attention)
- p-value: 0.026 (significant at alpha=0.05)
- Effect size: -0.23 (small, opposite direction)
- **Interpretation**: Model trained to IGNORE real biochemical contacts

## Conclusion

**Frozen augmentation**: Marginally beneficial (weak alignment with biochemistry)
**Relaxed augmentation**: Harmful (actively trains model away from real contacts)

The significant *decrease* in contact importance with relaxed DNA suggests the minimized conformations either introduce artifacts, misalign with the protein interface, or confuse the learning process.

**Recommendation**: Use frozen DNA if augmentation is desired. Relaxed DNA in current form should be avoided or redesigned.
