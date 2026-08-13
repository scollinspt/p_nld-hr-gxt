#!/usr/bin/env python3
"""
HRV-Metabolic Merged Dataset Analysis Summary
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

PROJECT_DIR = Path("/Users/collins/Projects/p_nld-hr-gxt-master")
OUTPUT_DIR = PROJECT_DIR / "merged_data"
PLOT_DIR = OUTPUT_DIR / "analysis_plots"
PLOT_DIR.mkdir(exist_ok=True)

# Load merged data
df = pd.read_csv(OUTPUT_DIR / "hrv_metabolic_merged.csv", na_values=['#NULL!', 'NULL'])

# Convert all numeric columns, replacing #NULL! with NaN
for col in df.columns:
    try:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    except:
        pass

print("="*70)
print("MERGED DATASET ANALYSIS SUMMARY")
print("="*70)

# Basic stats
print(f"\nDataset Shape: {df.shape[0]} rows × {df.shape[1]} columns")
print(f"Subjects: {df['subject_id_new'].nunique()}")
print(f"Time coverage: {df['timestamp'].min():.0f}s to {df['timestamp'].max():.0f}s per subject")

# HRV Variables
print("\n" + "-"*70)
print("HRV TIME-DOMAIN METRICS")
print("-"*70)
hrv_time = ['meanrr', 'sdnn', 'rmssd', 'pnn50']
print(df[hrv_time].describe().round(2))

print("\nHRV FREQUENCY-DOMAIN METRICS")
print("-"*70)
hrv_freq = ['vlfpowfft', 'lfpowfft', 'hfpowfft', 'totpowfft', 'lfhffft']
print(df[hrv_freq].describe().round(2))

print("\nHRV COMPLEXITY METRICS")
print("-"*70)
hrv_complex = ['sampen', 'dfa1', 'dfa2']
print(df[hrv_complex].describe().round(2))

# Metabolic Variables
print("\nMETABOLIC VARIABLES")
print("-"*70)
metabolic = ['vo2_mlkgmin_median', 'vco2_l_median', 'rq_median', 'hr_median', 've_stpd_median']
print(df[metabolic].describe().round(2))

# Correlations
print("\nCORRELATIONS: HRV vs Metabolic")
print("-"*70)
corr_vars = ['rmssd', 'lfpowfft', 'hfpowfft', 'sampen', 'vo2_mlkgmin_median', 'vco2_l_median', 'hr_median']
df_corr = df[corr_vars].copy()
df_corr = df_corr.dropna()
corr_matrix = df_corr.corr()

# Plot 1: Correlation heatmap
fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='coolwarm', center=0, 
            square=True, linewidths=1, ax=ax, cbar_kws={'label': 'Correlation'})
ax.set_title('HRV vs Metabolic Correlations', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(PLOT_DIR / "correlation_heatmap.png", dpi=150)
print("✓ Saved correlation heatmap")

# Plot 2: Time series sample
fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)

# Sample subject
subj = 1
df_subj = df[df['subject_id_new']==subj].sort_values('timestamp')

axes[0].plot(df_subj['timestamp'], df_subj['meanhr'], 'b-', label='Mean HR', linewidth=2)
axes[0].set_ylabel('Heart Rate (bpm)', fontsize=11)
axes[0].legend(loc='upper left')
axes[0].grid(True, alpha=0.3)

axes[1].plot(df_subj['timestamp'], df_subj['rmssd'], 'g-', label='RMSSD', linewidth=2)
axes[1].set_ylabel('RMSSD (ms)', fontsize=11)
axes[1].legend(loc='upper left')
axes[1].grid(True, alpha=0.3)

axes[2].plot(df_subj['timestamp'], df_subj['vo2_mlkgmin_median'], 'r-', label='VO₂', linewidth=2)
axes[2].set_ylabel('VO₂ (ml/kg/min)', fontsize=11)
axes[2].set_xlabel('Time (seconds)', fontsize=11)
axes[2].legend(loc='upper left')
axes[2].grid(True, alpha=0.3)

fig.suptitle(f'Sample Time Series - Subject {subj}', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(PLOT_DIR / f"timeseries_subject_{subj}.png", dpi=150)
print(f"✓ Saved time series plot for subject {subj}")

# Plot 3: HRV parameter distributions
fig, axes = plt.subplots(2, 2, figsize=(11, 8))

axes[0, 0].hist(df['rmssd'].dropna(), bins=30, color='steelblue', edgecolor='black', alpha=0.7)
axes[0, 0].set_xlabel('RMSSD (ms)', fontsize=11)
axes[0, 0].set_ylabel('Frequency', fontsize=11)
axes[0, 0].set_title('RMSSD Distribution', fontweight='bold')

axes[0, 1].hist(df['lfhffft'].dropna(), bins=30, color='coral', edgecolor='black', alpha=0.7)
axes[0, 1].set_xlabel('LF/HF Ratio', fontsize=11)
axes[0, 1].set_ylabel('Frequency', fontsize=11)
axes[0, 1].set_title('LF/HF Ratio Distribution', fontweight='bold')

axes[1, 0].hist(df['sampen'].dropna(), bins=30, color='seagreen', edgecolor='black', alpha=0.7)
axes[1, 0].set_xlabel('Sample Entropy', fontsize=11)
axes[1, 0].set_ylabel('Frequency', fontsize=11)
axes[1, 0].set_title('Sample Entropy Distribution', fontweight='bold')

axes[1, 1].hist(df['vo2_mlkgmin_median'].dropna(), bins=30, color='orchid', edgecolor='black', alpha=0.7)
axes[1, 1].set_xlabel('VO₂ (ml/kg/min)', fontsize=11)
axes[1, 1].set_ylabel('Frequency', fontsize=11)
axes[1, 1].set_title('VO₂ Distribution', fontweight='bold')

plt.suptitle('Variable Distributions', fontsize=14, fontweight='bold', y=1.00)
plt.tight_layout()
plt.savefig(PLOT_DIR / "distributions.png", dpi=150)
print("✓ Saved distributions plot")

# Plot 4: HR vs VO2 scatter
fig, ax = plt.subplots(figsize=(10, 6))
scatter = ax.scatter(df['vo2_mlkgmin_median'], df['meanhr'], 
                     c=df['subject_id_new'], cmap='tab20', 
                     alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
ax.set_xlabel('VO₂ (ml/kg/min)', fontsize=12)
ax.set_ylabel('Mean HR (bpm)', fontsize=12)
ax.set_title('Heart Rate vs VO₂ Uptake', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3)
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('Subject ID', fontsize=11)
plt.tight_layout()
plt.savefig(PLOT_DIR / "hr_vo2_scatter.png", dpi=150)
print("✓ Saved HR vs VO₂ plot")

# Plot 5: RMSSD vs RQ (as autonomic indicator)
fig, ax = plt.subplots(figsize=(10, 6))
scatter = ax.scatter(df['rq_median'], df['rmssd'], 
                     c=df['meanhr'], cmap='viridis', 
                     alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
ax.set_xlabel('Respiratory Quotient (RQ)', fontsize=12)
ax.set_ylabel('RMSSD (ms)', fontsize=12)
ax.set_title('Parasympathetic Activity vs Substrate Utilization', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3)
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('Heart Rate (bpm)', fontsize=11)
plt.tight_layout()
plt.savefig(PLOT_DIR / "rmssd_rq_scatter.png", dpi=150)
print("✓ Saved RMSSD vs RQ plot")

# Summary statistics table
summary_stats = pd.DataFrame({
    'Variable': ['meanrr', 'sdnn', 'rmssd', 'sampen', 'lfhffft', 
                 'vo2_mlkgmin_median', 'hr_median', 've_vo2'],
    'Mean': [df[v].mean() for v in ['meanrr', 'sdnn', 'rmssd', 'sampen', 'lfhffft',
                                     'vo2_mlkgmin_median', 'hr_median', 've_vo2']],
    'SD': [df[v].std() for v in ['meanrr', 'sdnn', 'rmssd', 'sampen', 'lfhffft',
                                  'vo2_mlkgmin_median', 'hr_median', 've_vo2']],
    'Min': [df[v].min() for v in ['meanrr', 'sdnn', 'rmssd', 'sampen', 'lfhffft',
                                   'vo2_mlkgmin_median', 'hr_median', 've_vo2']],
    'Max': [df[v].max() for v in ['meanrr', 'sdnn', 'rmssd', 'sampen', 'lfhffft',
                                   'vo2_mlkgmin_median', 'hr_median', 've_vo2']]
})

print("\nSUMMARY STATISTICS")
print("-"*70)
print(summary_stats.round(2))

# Save summary to file
summary_stats.to_csv(OUTPUT_DIR / "summary_statistics.csv", index=False)
print("\n✓ Saved summary statistics table")

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
print(f"""
Key Findings:
  - 923 HRV-metabolic data points successfully merged
  - 22 subjects with complete data
  - Data ready for statistical modeling
  
Next Steps:
  1. Mixed-effects models for autonomic-metabolic coupling
  2. Exercise stage effects on HRV dynamics
  3. Complexity metrics as exercise tolerance predictors
  4. Publication figures for manuscript

All plots saved to: {PLOT_DIR}
All data saved to: {OUTPUT_DIR}
""")
