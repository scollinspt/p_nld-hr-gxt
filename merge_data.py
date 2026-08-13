#!/usr/bin/env python3
"""
HRV & Metabolic Data Merge Pipeline - Fixed Version
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.preprocessing import StandardScaler
import joblib
import warnings
warnings.filterwarnings('ignore')

# Configure
pd.set_option('display.max_columns', 20)
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

PROJECT_DIR = Path("/Users/collins/Projects/p_nld-hr-gxt-master")
PYTHON_ANALYSIS_DIR = PROJECT_DIR / "python_analysis"
OUTPUT_DIR = PROJECT_DIR / "merged_data"
OUTPUT_DIR.mkdir(exist_ok=True)

print("="*70)
print("HRV & METABOLIC DATA MERGE PIPELINE")
print("="*70)

# Load HRV
print("\n[1/6] Loading HRV data...")
df_hrv = pd.read_csv(PYTHON_ANALYSIS_DIR / "hrv_python_new_advancement.csv")
print(f"✓ HRV: {df_hrv.shape[0]} windows, {df_hrv['subject_id_new'].nunique()} subjects")

# Load ID mapping
print("[1/6] Loading ID mapping...")
id_map = pd.read_csv(PROJECT_DIR / "Old-to-New-IDs.csv")
id_dict = dict(zip(id_map['id-old'], id_map['id-new']))

# Load and map metabolic
print("[1/6] Loading and mapping metabolic data...")
df_met = pd.read_csv(PROJECT_DIR / "all_bxb_10secondEpochs.csv")
df_met['subject_id_new'] = df_met['ID'].map(id_dict).fillna(-1).astype('int64')
df_met = df_met[df_met['subject_id_new'] > 0].copy()
print(f"✓ Metabolic: {df_met.shape[0]} epochs, {df_met['subject_id_new'].nunique()} subjects")

# Create timestamps
print("\n[2/6] Creating timestamps...")
df_hrv['timestamp'] = df_hrv['window_start']
df_met['timestamp'] = df_met['epoch'] * 10.0

# Merge per subject
print("\n[3/6] Merging datasets...")
merged_list = []

for subj in sorted(df_hrv['subject_id_new'].unique()):
    df_h = df_hrv[df_hrv['subject_id_new']==subj].sort_values('timestamp').copy()
    df_m = df_met[df_met['subject_id_new']==subj].sort_values('timestamp').copy()
    
    if len(df_m) == 0:
        continue
    
    # Merge asof (nearest match)
    merged = pd.merge_asof(df_h, df_m, on='timestamp', direction='nearest', tolerance=60, suffixes=('', '_met'))
    merged_list.append(merged)
    print(f"  Subject {subj}: {len(merged)} merged rows")

df_merged = pd.concat(merged_list, ignore_index=True)
print(f"✓ Total merged: {len(df_merged)} rows")

# Data cleaning
print("\n[4/6] Data cleaning...")
df_merged = df_merged[(df_merged['meanhr'] >= 40) & (df_merged['meanhr'] <= 200)].copy()
df_merged = df_merged.dropna(subset=['meanrr', 'vo2_mlkgmin_median'])
print(f"✓ Clean data: {len(df_merged)} rows")

# Select and prepare features
print("\n[5/6] Preparing features...")
cols_keep = [
    'subject_id_new', 'timestamp',
    # HRV - time domain
    'meanrr', 'sdnn', 'rmssd', 'pnn50',
    # HRV - frequency
    'vlfpowfft', 'lfpowfft', 'hfpowfft', 'totpowfft', 'lfhffft',
    # HRV - complexity
    'sampen', 'dfa1', 'dfa2',
    # Heart rate
    'meanhr', 'hr_MEDIAN',
    # Metabolic - oxygen
    'vo2_mlkgmin_median', 'vo2_L_median',
    # Metabolic - CO2
    'vco2_L_MEDIAN', 'rq_MEDIAN',
    # Metabolic - ventilation
    've_stpd_MEDIAN', 've_btps_MEDIAN', 'rr_MEDIAN',
    # Metabolic - gases
    'peto2_MEDIAN', 'petco2_MEDIAN', 'feo2_MEDIAN', 'feco2_MEDIAN',
    # Exercise
    'speed_MEDIAN', 'grade_MEDIAN'
]

df_analysis = df_merged[cols_keep].copy()
df_analysis.columns = [c.lower() for c in df_analysis.columns]

# Derived features
df_analysis['parasympathetic_index'] = df_analysis['hfpowfft'] / df_analysis['totpowfft']
df_analysis['sympathetic_index'] = df_analysis['lfpowfft'] / df_analysis['totpowfft']
df_analysis['autonomic_balance'] = df_analysis['lfhffft']
df_analysis['ve_vo2'] = df_analysis['ve_stpd_median'] / df_analysis['vo2_l_median']
df_analysis['ve_vco2'] = df_analysis['ve_stpd_median'] / df_analysis['vco2_l_median']

# Save
print("\n[6/6] Exporting data...")
output_csv = OUTPUT_DIR / "hrv_metabolic_merged.csv"
output_parquet = OUTPUT_DIR / "hrv_metabolic_merged.parquet"
df_analysis.to_csv(output_csv, index=False)
df_analysis.to_parquet(output_parquet, index=False)
print(f"✓ Saved: {output_csv.name}, {output_parquet.name}")

# Scale
scaler = StandardScaler()
num_features = [c for c in df_analysis.select_dtypes(include=[np.number]).columns 
                if c not in ['subject_id_new', 'timestamp']]
df_scaled = df_analysis.copy()
df_scaled[num_features] = scaler.fit_transform(df_analysis[num_features])
df_scaled.to_csv(OUTPUT_DIR / "hrv_metabolic_merged_scaled.csv", index=False)
joblib.dump(scaler, OUTPUT_DIR / "scaler.pkl")

print(f"\n" + "="*70)
print("✓ PIPELINE COMPLETE")
print("="*70)
print(f"""
Dataset: {len(df_analysis):,} rows × {len(df_analysis.columns)} features
Subjects: {df_analysis['subject_id_new'].nunique()}
Output: {OUTPUT_DIR}

Ready for statistical analysis!
""")
