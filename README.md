# Sex Differences in Human Subcortical Anatomy

Analysis code for:

> Lee HM, Liu S, Lalonde FM, Clasen L, Guma E, DeCasien A, Moraczewski D, Thomas A, Lerch JP, Devenyi GA, Chakravarty MM, Berman KF, Delaney A, Raznahan A. **Sex Differences in Human Subcortical Anatomy are Robust and Replicable in Adulthood, Already Evident Prepubertally, and Aligned with Sex Chromosomal and Gonadal Influences.** *bioRxiv* (2026). https://doi.org/10.64898/2026.10.05.756728

## Overview

We mapped sex differences in the volume and shape of the striatum, pallidum, thalamus, hippocampus and amygdala in three cohorts (HCP-YA, UK Biobank, ABCD; combined N = 2,886). We then compared these maps with the anatomical effects of X- and Y-chromosome dosage and testicular hormone production from three clinical cohorts (combined N = 313). This repository contains the scripts for the statistical analyses and main figures.

## Scripts

| Script | Analysis |
| --- | --- |
| `hcp_stats.m` | Sex differences in subcortical volume and shape in early adulthood (HCP-YA); comparison of size-correction methods; Moran spectral randomization null maps |
| `ukb_stats.m` | The same models in an independent mid-adulthood cohort (UK Biobank, ages 45–50) |
| `replication.m` | Within-sample reproducibility (sex-balanced split-half subsampling, n = 50–450) and between-sample replicability (HCP-YA vs. UKB) |
| `xyt_stats.m` | Effects of X-chromosome dosage (XXY), Y-chromosome dosage (XYY) and testicular hormones (IGD), and their spatial congruence with sex differences |
| `abcd_stats.m` | Subcortical and cortical sex differences in pre/early puberty (ABCD, age 9) and comparison with adulthood |

## Dependencies

- MATLAB R2016b or later, with the Statistics and Machine Learning Toolbox and the Parallel Computing Toolbox
- [SurfStat](https://www.math.mcgill.ca/keith/surfstat/) (`SurfStatLinMod`, `SurfStatT`, `term`)
- [BrainSpace](https://github.com/MICA-MNI/BrainSpace) (`surface_to_graph`, `compute_mem`, `moran_randomization`)
- `fdr_bh` (Benjamini–Hochberg FDR) and `hexscatter` (MATLAB File Exchange)
- `daviolinplot` ([DataViz](https://github.com/povilaskarvelis/DataViz))
- Subcortical segmentation and surface generation used [MAGeT-Brain](https://github.com/CoBrALab/MAGeTbrain); you need the template surfaces to compute vertex-wise surface area.

Some custom visualization functions (`loadColorMaps`, `ThreshColorMap`, `HyoMagetView*`, `HyoSurfacePlotFunction`) are not included. They only affect figure rendering, not the statistics.

## Data

Participant-level data are not included in this repository.

- **HCP-YA** (S1200 release): https://db.humanconnectome.org/, available to registered users who accept the data use terms
- **UK Biobank**: https://www.ukbiobank.ac.uk, available to approved researchers via application
- **ABCD** (Release 2.0.0): NIMH Data Archive, https://nda.nih.gov/
- **NIH clinical cohorts** (sex chromosome aneuploidies, isolated gonadotropin-releasing hormone deficiency): not publicly available because of participant privacy; available from the corresponding author on reasonable request with a data use agreement

Group-level statistical maps will be provided as source data on publication.

## License

Released under the MIT License (see `LICENSE`).

## Acknowledgments

Supported by the Intramural Research Program of the National Institute of Mental Health (1ZIAMH002949).

## Contact

Hyo Min Lee — hyo.lee@nih.gov
