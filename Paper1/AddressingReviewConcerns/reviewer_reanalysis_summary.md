# Reviewer Reanalysis Summary

## Scope
This separate reanalysis preserves the validated dataset and prior outputs. It tests window nonstationarity, full-window metabolic alignment, shorter windows, entropy detrending, respiration, and participant-aware inference.

## VMax inventory
The processed VMax file contains 10-s epoch medians for VO2, VCO2, VE (STPD/BTPS), respiratory rate, tidal volume, RQ, VE/VO2, VE/VCO2, oxygen pulse, end-tidal gases, speed, and grade. No adjudicated ventilatory/anaerobic-threshold time marker was found.

## Within-window nonstationarity
Across 300-s windows, VO2 end-minus-start median was 21.55 ml kg-1 min-1 (IQR 17.15 to 24.55); mean 20.74 (SD 5.74). 99.6% of windows changed by at least 5 ml kg-1 min-1. This quantifies material drift rather than local stationarity for many segments.

## Alignment
Start, midpoint, mean, median, and end VO2 models were separately refit; inspect `alignment_sensitivity_models.csv` before retaining any claim. Whole-window mean/median VO2 are the scientifically preferable exposures for a window summary.

## Short windows and SampEn
Short-window time-domain and SampEn analyses are in `window_length_sensitivity_models.csv` and `sampen_detrending_sensitivity.csv`. Frequency-domain LF/HF is flagged as defensible only for the 300-s window because low-frequency estimation requires longer segments. Interpret 30--180-s spectral measures as exploratory only.

## Respiration and HF
81.5% of 300-s windows had mean respiratory frequency above 0.40 Hz. The median participant-specific VO2 at first crossing was 25.03 ml kg-1 min-1 (N=22). This means the fixed 0.15--0.40 Hz HF band may omit respiratory-linked variability at higher intensity; conventional HF and LF/HF cannot be interpreted without this limitation.

## Inference
GEE fits use participant clusters and grid-estimated AR(1) working correlation with robust sandwich covariance provided by statsmodels GEE. With 22 clusters, asymptotic sandwich standard errors may still be optimistic; mixed models and thinning are included as sensitivity analyses, not confirmation from 923 independent units.

## Scientific decision
1. A 300-s window is not uniformly locally stationary during the GXT; drift magnitude must be reported and whole-window exposures should be preferred.
2. Alignment sensitivity results, not prior start-time-only results, determine whether conclusions survive.
3. SampEn claims are conditional on detrended and shorter-window results; do not call it physiological complexity unless both are robust.
4. Respiratory frequency above 0.40 Hz limits conventional HF/LF-HF interpretation at high intensity.
5. The final Paper 1 claims should be limited to findings stable across whole-window alignment, detrending, window length, and participant-aware sensitivity models.
