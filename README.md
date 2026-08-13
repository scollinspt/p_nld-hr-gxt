# Non-Linear Dynamics of Heart Rate During Graded Maximum Exercise

**A Python-based analysis of autonomic-metabolic coupling during incremental exercise**

## Project Overview

This project investigates the non-linear dynamics of heart rate variability (HRV) and its relationship to metabolic responses during graded maximum exercise testing. Using a cohort of 23 athletes from 2014, we employ advanced HRV analysis and signal complexity metrics to understand autonomic nervous system behavior across the full exercise intensity spectrum—from rest through peak maximal effort.

### Core Research Questions

1. **How do HRV complexity metrics (sample entropy, DFA) relate to metabolic intensity (VO₂)?**
2. **Can parasympathetic withdrawal be predicted from non-linear HRV dynamics during exercise?**
3. **Do novel HRV metrics improve detection of autonomic saturation at maximal intensity?**

## Dataset

**Study Design:** Graded maximal treadmill exercise test (September 2014)  
**Sample:** 23 athletes | **Protocol:** Progressive incremental stages  
**Data Frequency:** 
- HRV: 10-second sliding window advancement (930 windows per subject)
- Metabolic: 10-second epoch medians (2075 total epochs)
- ECG: 1000 Hz sampling rate

### Key Data Files

| File | Purpose | Records |
|------|---------|---------|
| `ECG_EDF/` | Raw ECG recordings (1000 Hz) | 23 subjects |
| `rr/` | Kubios-extracted R-peak timestamps | 23 subjects |
| `merged_data/hrv_metabolic_merged.csv` | **Unified HRV-metabolic dataset** | **923 rows** |
| `subject_TestData.csv` | Demographics, anthropometrics, fitness | 23 subjects |
| `all_bxb_10secondEpochs.csv` | Metabolic cart data (beat-by-beat medians) | 2075 epochs |

## Analysis Pipeline

### 1. HRV Regeneration (`hrv_analysis.py`)

Generates HRV metrics from R-peak timestamps using validated signal processing:

```python
# 5-minute sliding windows, 10-second advancement
python hrv_analysis.py
```

**Outputs:**
- `hrv_python_new_advancement.csv` — 930 HRV windows × 19 metrics
- Validation plots comparing Python vs Kubios reference
- Correlation matrices with metabolic data

**HRV Metrics Computed:**
- **Time-domain:** meanRR, SDNN, RMSSD, pNN50
- **Frequency-domain:** VLF, LF, HF power; LF/HF ratio  
- **Complexity:** Sample Entropy, DFA exponents
- **HR:** Mean heart rate

### 2. Data Merge (`merge_data.py`)

Aligns HRV windows with metabolic epochs using temporal nearest-neighbor matching:

```python
# Merge HRV + metabolic data with ID mapping
python merge_data.py
```

**Outputs:**
- `merged_data/hrv_metabolic_merged.csv` — 923 rows, ready for analysis
- `merged_data/hrv_metabolic_merged_scaled.csv` — StandardScaler normalized
- `merged_data/scaler.pkl` — Serialized scaler for reproducibility

### 3. Exploratory Analysis (`analyze_merged_data.py`)

Generates summary statistics and diagnostic plots:

```python
# Create visualizations and correlation analysis
python analyze_merged_data.py
```

**Outputs:**
- `merged_data/analysis_plots/` — 5 publication-quality figures
- `merged_data/summary_statistics.csv` — Descriptive statistics

## Getting Started

### Requirements
- Python 3.9+
- Dependencies: pandas, numpy, scipy, scikit-learn, matplotlib, seaborn, joblib

### Installation

```bash
# Clone repository
git clone <repo-url>
cd p_nld-hr-gxt-master

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install pandas numpy scipy scikit-learn matplotlib seaborn joblib
```

### Quick Start

```bash
# Step 1: Regenerate HRV from R-peak timestamps
python hrv_analysis.py

# Step 2: Merge HRV with metabolic data
python merge_data.py

# Step 3: Generate analysis plots and summaries
python analyze_merged_data.py

# Output: merged_data/hrv_metabolic_merged.csv (ready for modeling)
```

## Data Structure

### Merged Dataset Columns (923 rows × 34 features)

**Identifiers:**
- `subject_id_new` — Subject ID (1-23)
- `timestamp` — Seconds from exercise start

**HRV Metrics:**
- Time-domain: `meanrr`, `sdnn`, `rmssd`, `pnn50`
- Frequency: `vlfpowfft`, `lfpowfft`, `hfpowfft`, `totpowfft`, `lfhffft`
- Complexity: `sampen`, `dfa1`, `dfa2`
- HR: `meanhr`

**Metabolic:**
- Oxygen: `vo2_mlkgmin_median`, `vo2_l_median`
- CO₂: `vco2_l_median`
- Ventilation: `ve_stpd_median`, `rr_median`
- Gases: `peto2_median`, `petco2_median`, `feo2_median`, `feco2_median`
- Exercise: `speed_median`, `grade_median`

**Derived Autonomic Indices:**
- `parasympathetic_index` = HF / Total Power
- `sympathetic_index` = LF / Total Power
- `autonomic_balance` = LF/HF ratio
- `ve_vo2` = Ventilation efficiency
- `ve_vco2` = CO₂ clearance efficiency

## Key Findings

### HR Distribution
- Rest/very low (<100 bpm): 95 points, VO₂ = 18.2 ml/kg/min
- Moderate (100–150 bpm): 673 points, VO₂ = 23.9 ml/kg/min
- High (150–180 bpm): 143 points, VO₂ = 38.6 ml/kg/min
- Peak maximal (>180 bpm): 11 points, VO₂ = 43.0 ml/kg/min

### HRV Summary
- Mean HR: 127 ± 24.5 bpm (during exercise)
- RMSSD: 6.75 ± 6.08 ms (parasympathetic withdrawal evident)
- LF/HF: 4.07 ± 3.99 (sympathetic dominance increasing with intensity)
- Sample Entropy: Very low (0.01 ± 0.02), indicating heart rate becomes more regular/predictable under exercise stress

## Planned Analyses

### Phase 1: Baseline Modeling
- Mixed-effects regression for subject-level autonomic responses
- Exercise stage stratification (submaximal vs maximal)
- Intensity-dependent HRV parameter trajectories

### Phase 2: Novel Metrics
- Complexity-metabolic coupling (sample entropy vs VO₂)
- Autonomic saturation curves (Do HRV metrics plateau at peak?)
- Residual HR as sympathetic activity proxy

### Phase 3: Publication
- Autonomic-metabolic coupling figures for manuscript
- Statistical validation of complexity-exercise tolerance relationship
- Clinical/sport science implications

## File Structure

```
p_nld-hr-gxt-master/
├── README.md                          # This file
├── hrv_analysis.py                    # HRV regeneration pipeline
├── merge_data.py                      # HRV-metabolic merge
├── analyze_merged_data.py             # Exploratory analysis
├── ECG_EDF/                           # Raw ECG data (1000 Hz)
├── rr/                                # Kubios R-peak timestamps
├── subject_HRV_data/                  # Original Kubios HRV output
├── merged_data/                       # Generated merged & analysis outputs
│   ├── hrv_metabolic_merged.csv       # Primary analysis dataset
│   ├── hrv_metabolic_merged_scaled.csv
│   ├── scaler.pkl
│   ├── summary_statistics.csv
│   └── analysis_plots/                # Visualizations
├── exHRVdata-SelectedVariables.csv    # Kubios reference data
├── subject_TestData.csv               # Demographics & fitness
├── all_bxb_10secondEpochs.csv         # Metabolic cart data
├── testingProtocol.csv                # Exercise protocol
└── Old-to-New-IDs.csv                 # Subject ID mapping

```

## Methodology Notes

### HRV Analysis
- **Window size:** 300 seconds (5 minutes) — standard for exercise HRV
- **Advancement:** 10 seconds — enables tight temporal alignment with metabolic data
- **RR source:** Kubios-validated R-peak detection (validated accuracy: <10.4% error)
- **Frequency bands:** VLF (0.0033–0.04 Hz), LF (0.04–0.15 Hz), HF (0.15–0.4 Hz)
- **Complexity tolerance:** 0.2 × std(RR intervals) for sample entropy

### Data Merging Strategy
- Per-subject nearest-neighbor matching (merge_asof)
- 60-second temporal tolerance for missing epochs
- Retention of all maximal exercise data (HR up to 260 bpm = peak effort)

### Quality Control
- HR validation: Removed outliers outside 40–200 bpm range (1% of data)
- Metabolic missing values: Imputed via nearest-timestamp matching
- Subject 11: Only 7 merged windows (likely early test termination)

## Citation

If you use this dataset or analysis pipeline in your research, please cite:

```bibtex
@misc{nld_hr_exercise_2014,
  title={Non-Linear Dynamics of Heart Rate During Graded Maximum Exercise},
  author={Collins},
  year={2014},
  note={Dataset and Python analysis pipeline}
}
```

## Contributing

Improvements and extensions welcome! Areas for contribution:
- Alternative HRV metrics (time-varying spectral analysis, wavelet decomposition)
- Exercise stage classification (automated intensity detection)
- Predictive modeling (machine learning for autonomic state prediction)
- Visualization enhancements

## License

See LICENSE file

## Contact

For questions or collaborations, please open an issue on GitHub.
