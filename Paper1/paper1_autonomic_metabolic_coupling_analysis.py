"""
PAPER 1: Non-Linear Dynamics of Heart Rate Variability Reveal Autonomic-Metabolic Coupling 
During Incremental Exercise

Statistical Analysis Framework
- Exercise stage stratification
- Autonomic-metabolic coupling analysis
- Publication-quality figures and LaTeX-formatted tables

Author: [Your Name]
Date: 2026
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

# Configure plotting for publication quality
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")
plt.rcParams.update({
    'font.size': 10,
    'font.family': 'sans-serif',
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.dpi': 300,
})

class Paper1Analysis:
    """Complete statistical analysis for Paper 1: Autonomic-Metabolic Coupling"""
    
    def __init__(self, data_path='merged_data/hrv_metabolic_merged.csv'):
        """Initialize analysis with merged HRV-metabolic data"""
        self.df = pd.read_csv(data_path)
        self.results = {}
        self.figures = {}
        print(f"Loaded {len(self.df)} HRV windows from {len(self.df['subject_id_new'].unique())} subjects")
        
    def create_intensity_stages(self):
        """Stratify data into exercise intensity stages based on VO2"""
        # Define thresholds based on VO2 percentiles
        vo2_col = 'vo2_mlkgmin_median'
        vo2_quantiles = self.df[vo2_col].quantile([0.1, 0.3, 0.7])
        
        def classify_stage(vo2):
            if pd.isna(vo2):
                return None
            if vo2 < vo2_quantiles[0.1]:
                return 'Rest'
            elif vo2 < vo2_quantiles[0.3]:
                return 'Moderate'
            elif vo2 < vo2_quantiles[0.7]:
                return 'High'
            else:
                return 'Maximal'
        
        self.df['intensity_stage'] = self.df[vo2_col].apply(classify_stage)
        
        # Print intensity stage distribution
        print("\n=== EXERCISE INTENSITY STAGE DISTRIBUTION ===")
        for stage in ['Rest', 'Moderate', 'High', 'Maximal']:
            mask = self.df['intensity_stage'] == stage
            n = mask.sum()
            if n > 0:
                vo2_range = self.df[mask][vo2_col].agg(['min', 'max'])
                print(f"{stage:10s}: n={n:4d}  VO2={vo2_range['min']:5.1f}–{vo2_range['max']:5.1f} ml/kg/min")
        
        return self.df
    
    def autonomic_summary(self):
        """Print autonomic indices by stage"""
        print("\n=== AUTONOMIC NERVOUS SYSTEM SUMMARY BY INTENSITY ===\n")
        for stage in ['Rest', 'Moderate', 'High', 'Maximal']:
            stage_data = self.df[self.df['intensity_stage'] == stage]
            print(f"{stage}:")
            print(f"  Heart Rate:     {stage_data['meanhr'].mean():6.1f} ± {stage_data['meanhr'].std():5.1f} bpm")
            print(f"  RMSSD:          {stage_data['rmssd'].mean():6.1f} ± {stage_data['rmssd'].std():5.1f} ms")
            print(f"  LF/HF Ratio:    {stage_data['lfhffft'].mean():6.1f} ± {stage_data['lfhffft'].std():5.1f}")
            print(f"  Sample Entropy: {stage_data['sampen'].mean():6.4f} ± {stage_data['sampen'].std():5.4f}")
            print()
    
    def descriptive_statistics(self):
        """Generate comprehensive descriptive statistics by intensity stage"""
        print("\n=== DESCRIPTIVE STATISTICS BY INTENSITY STAGE ===\n")
        
        desc_vars = {
            'Heart Rate (bpm)': 'meanhr',
            'VO2 (ml/kg/min)': 'vo2_mlkgmin_median',
            'VCO2 (L/min)': 'vco2_l_median',
            'VE (L/min)': 've_stpd_median',
            'SDNN (ms)': 'sdnn',
            'RMSSD (ms)': 'rmssd',
            'LF/HF Ratio': 'lfhffft',
            'Sample Entropy': 'sampen',
            'DFA α1': 'dfa1',
        }
        
        results_list = []
        for stage in ['Rest', 'Moderate', 'High', 'Maximal']:
            stage_data = self.df[self.df['intensity_stage'] == stage]
            for label, var in desc_vars.items():
                if var in stage_data.columns:
                    valid = stage_data[var].dropna()
                    results_list.append({
                        'Stage': stage,
                        'Variable': label,
                        'N': len(valid),
                        'Mean': valid.mean(),
                        'SD': valid.std(),
                        'Min': valid.min(),
                        'Max': valid.max(),
                    })
        
        desc_df = pd.DataFrame(results_list)
        self._save_latex_table(desc_df, 'table_descriptive_statistics.tex', 
                              caption='Descriptive Statistics by Exercise Intensity Stage')
        
        print(desc_df.to_string(index=False))
        return desc_df
    
    def correlation_analysis(self):
        """Correlation analysis between HRV metrics and metabolic variables"""
        print("\n=== CORRELATION ANALYSIS: HRV vs METABOLIC VARIABLES ===\n")
        
        hrv_vars = ['meanhr', 'sdnn', 'rmssd', 'hfpowfft', 'lfpowfft', 'lfhffft', 'sampen', 'dfa1']
        metabolic_vars = ['vo2_mlkgmin_median', 'vco2_l_median', 've_stpd_median', 'rq_median']
        
        correlation_results = []
        for hrv in hrv_vars:
            for met in metabolic_vars:
                if hrv in self.df.columns and met in self.df.columns:
                    valid = self.df[[hrv, met]].dropna()
                    if len(valid) > 10:
                        r, p = pearsonr(valid[hrv], valid[met])
                        correlation_results.append({
                            'HRV Variable': hrv,
                            'Metabolic': met,
                            'r': r,
                            'p': p,
                            'Sig': 'Yes' if p < 0.05 else 'No',
                        })
        
        corr_df = pd.DataFrame(correlation_results).sort_values('p')
        sig_corr = corr_df[corr_df['p'] < 0.05]
        
        print(f"Significant correlations (p < 0.05): {len(sig_corr)}/{len(corr_df)}")
        print(sig_corr.head(12).to_string(index=False))
        
        self._save_latex_table(sig_corr.head(15), 'table_correlations_significant.tex',
                              caption='Significant HRV-Metabolic Correlations (p<0.05)')
        return corr_df
    
    def intensity_effects(self):
        """Test for significant effects across intensity stages"""
        print("\n=== INTENSITY STAGE EFFECTS (One-Way ANOVA) ===\n")
        
        test_vars = ['meanhr', 'sdnn', 'rmssd', 'hfpowfft', 'lfpowfft', 'lfhffft', 'sampen', 'dfa1']
        intensity_effects = []
        
        for var in test_vars:
            if var in self.df.columns:
                rest = self.df[self.df['intensity_stage'] == 'Rest'][var].dropna()
                moderate = self.df[self.df['intensity_stage'] == 'Moderate'][var].dropna()
                high = self.df[self.df['intensity_stage'] == 'High'][var].dropna()
                maximal = self.df[self.df['intensity_stage'] == 'Maximal'][var].dropna()
                
                if len(rest) > 2 and len(moderate) > 2 and len(high) > 2 and len(maximal) > 2:
                    f_stat, p_val = stats.f_oneway(rest, moderate, high, maximal)
                    
                    # Effect size (eta-squared)
                    all_data = np.concatenate([rest.values, moderate.values, high.values, maximal.values])
                    grand_mean = np.mean(all_data)
                    ss_between = (len(rest) * (np.mean(rest) - grand_mean)**2 + 
                                 len(moderate) * (np.mean(moderate) - grand_mean)**2 +
                                 len(high) * (np.mean(high) - grand_mean)**2 +
                                 len(maximal) * (np.mean(maximal) - grand_mean)**2)
                    ss_total = np.sum((all_data - grand_mean)**2)
                    eta_sq = ss_between / ss_total if ss_total > 0 else 0
                    
                    intensity_effects.append({
                        'Variable': var,
                        'F': f_stat,
                        'p': p_val,
                        'Sig': 'Yes' if p_val < 0.05 else 'No',
                        'η²': eta_sq
                    })
        
        anova_df = pd.DataFrame(intensity_effects).sort_values('p')
        print(anova_df[['Variable', 'F', 'p', 'η²']].to_string(index=False))
        
        self._save_latex_table(anova_df, 'table_anova_intensity.tex',
                              caption='Effects of Exercise Intensity (One-Way ANOVA)')
        return anova_df
    
    def plot_autonomic_metabolic_coupling(self):
        """Figure 1: Autonomic-Metabolic Coupling"""
        fig, axes = plt.subplots(2, 3, figsize=(14, 8))
        fig.suptitle('Autonomic-Metabolic Coupling Across Exercise Intensity Stages', 
                     fontsize=14, fontweight='bold')
        
        pairs = [
            ('lfhffft', 'vo2_mlkgmin_median', axes[0, 0]),
            ('sampen', 'vo2_mlkgmin_median', axes[0, 1]),
            ('rmssd', 'vco2_l_median', axes[0, 2]),
            ('dfa1', 'vo2_mlkgmin_median', axes[1, 0]),
            ('hfpowfft', 'vo2_mlkgmin_median', axes[1, 1]),
            ('meanhr', 'vo2_mlkgmin_median', axes[1, 2]),
        ]
        
        colors = {'Rest': '#2ecc71', 'Moderate': '#f39c12', 'High': '#e74c3c', 'Maximal': '#8e44ad'}
        
        for hrv_var, met_var, ax in pairs:
            if hrv_var in self.df.columns and met_var in self.df.columns:
                for stage in ['Rest', 'Moderate', 'High', 'Maximal']:
                    stage_data = self.df[self.df['intensity_stage'] == stage]
                    valid = stage_data[[hrv_var, met_var]].dropna()
                    
                    if len(valid) > 0:
                        ax.scatter(valid[met_var], valid[hrv_var], alpha=0.4, s=25, 
                                  color=colors[stage], label=stage)
                
                ax.set_xlabel('Metabolic Variable', fontsize=9)
                ax.set_ylabel(hrv_var, fontsize=9)
                ax.grid(True, alpha=0.2)
        
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.01), ncol=4)
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.98])
        plt.savefig('figure_1_autonomic_metabolic_coupling.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved: figure_1_autonomic_metabolic_coupling.png")
    
    def plot_intensity_effects(self):
        """Figure 2: HRV by Intensity"""
        fig, axes = plt.subplots(2, 3, figsize=(14, 8))
        fig.suptitle('Heart Rate Variability Across Exercise Intensity', fontsize=14, fontweight='bold')
        
        vars_to_plot = [
            ('rmssd', 'RMSSD (ms)', axes[0, 0]),
            ('sdnn', 'SDNN (ms)', axes[0, 1]),
            ('hfpowfft', 'HF Power (ms²)', axes[0, 2]),
            ('lfhffft', 'LF/HF Ratio', axes[1, 0]),
            ('sampen', 'Sample Entropy', axes[1, 1]),
            ('dfa1', 'DFA α1', axes[1, 2]),
        ]
        
        stage_order = ['Rest', 'Moderate', 'High', 'Maximal']
        colors = ['#2ecc71', '#f39c12', '#e74c3c', '#8e44ad']
        
        for var, ylabel, ax in vars_to_plot:
            if var in self.df.columns:
                data = [self.df[self.df['intensity_stage'] == s][var].dropna() for s in stage_order]
                bp = ax.boxplot(data, labels=stage_order, patch_artist=True)
                
                for patch, color in zip(bp['boxes'], colors):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.6)
                
                ax.set_ylabel(ylabel, fontsize=10)
                ax.grid(True, axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('figure_2_intensity_effects.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved: figure_2_intensity_effects.png")
    
    def plot_complexity_progression(self):
        """Figure 3: Complexity Metrics"""
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle('Non-Linear Heart Rate Dynamics During Exercise', fontsize=14, fontweight='bold')
        
        stage_order = ['Rest', 'Moderate', 'High', 'Maximal']
        colors = ['#2ecc71', '#f39c12', '#e74c3c', '#8e44ad']
        x_pos = np.arange(len(stage_order))
        
        # Sample Entropy
        ax = axes[0]
        if 'sampen' in self.df.columns:
            means = [self.df[self.df['intensity_stage'] == s]['sampen'].mean() for s in stage_order]
            sds = [self.df[self.df['intensity_stage'] == s]['sampen'].std() for s in stage_order]
            ax.bar(x_pos, means, yerr=sds, capsize=5, color=colors, alpha=0.7, edgecolor='black')
            ax.set_ylabel('Sample Entropy', fontsize=11)
            ax.grid(True, axis='y', alpha=0.3)
        
        # DFA α1
        ax = axes[1]
        if 'dfa1' in self.df.columns:
            means = [self.df[self.df['intensity_stage'] == s]['dfa1'].mean() for s in stage_order]
            sds = [self.df[self.df['intensity_stage'] == s]['dfa1'].std() for s in stage_order]
            ax.bar(x_pos, means, yerr=sds, capsize=5, color=colors, alpha=0.7, edgecolor='black')
            ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, label='White Noise (α=1.0)')
            ax.set_ylabel('DFA Scaling Exponent (α1)', fontsize=11)
            ax.legend()
            ax.grid(True, axis='y', alpha=0.3)
        
        for ax in axes:
            ax.set_xticks(x_pos)
            ax.set_xticklabels(stage_order)
            ax.set_xlabel('Exercise Intensity', fontsize=11)
        
        plt.tight_layout()
        plt.savefig('figure_3_complexity_progression.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved: figure_3_complexity_progression.png")
    
    def plot_subject_trajectories(self):
        """Figure 4: Subject Trajectories"""
        subject_counts = self.df.groupby('subject_id_new').size()
        top_subjects = subject_counts.nlargest(4).index.tolist()
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Individual HRV Trajectories During Incremental Exercise', fontsize=14, fontweight='bold')
        axes = axes.flatten()
        
        for idx, subject in enumerate(top_subjects):
            subject_data = self.df[self.df['subject_id_new'] == subject].sort_values('vo2_mlkgmin_median')
            
            ax = axes[idx]
            ax2 = ax.twinx()
            
            line1 = ax.plot(subject_data['vo2_mlkgmin_median'], subject_data['rmssd'], 
                           'o-', color='#3498db', linewidth=2, markersize=4, label='RMSSD')
            ax.set_ylabel('RMSSD (ms)', color='#3498db')
            ax.tick_params(axis='y', labelcolor='#3498db')
            
            line2 = ax2.plot(subject_data['vo2_mlkgmin_median'], subject_data['lfhffft'], 
                            's-', color='#e74c3c', linewidth=2, markersize=4, label='LF/HF')
            ax2.set_ylabel('LF/HF Ratio', color='#e74c3c')
            ax2.tick_params(axis='y', labelcolor='#e74c3c')
            
            ax.set_xlabel('VO₂ (ml/kg/min)', fontsize=10)
            ax.set_title(f'Subject {subject}', fontweight='bold')
            ax.grid(True, alpha=0.3)
            
            lines = line1 + line2
            labels = [l.get_label() for l in lines]
            ax.legend(lines, labels, loc='upper left', fontsize=8)
        
        plt.tight_layout()
        plt.savefig('figure_4_subject_trajectories.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved: figure_4_subject_trajectories.png")
    
    def plot_distribution_by_intensity(self):
        """Figure 5: Distributions by Intensity"""
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle('Distribution of Variables by Exercise Intensity', fontsize=14, fontweight='bold')
        
        vars_to_plot = [
            ('meanhr', 'Heart Rate (bpm)', axes[0, 0]),
            ('vo2_mlkgmin_median', 'VO₂ (ml/kg/min)', axes[0, 1]),
            ('lfhffft', 'LF/HF Ratio', axes[1, 0]),
            ('hfpowfft', 'HF Power (ms²)', axes[1, 1]),
        ]
        
        stage_order = ['Rest', 'Moderate', 'High', 'Maximal']
        colors = ['#2ecc71', '#f39c12', '#e74c3c', '#8e44ad']
        
        for var, ylabel, ax in vars_to_plot:
            if var in self.df.columns:
                data = [self.df[self.df['intensity_stage'] == s][var].dropna() for s in stage_order]
                parts = ax.violinplot(data, positions=range(len(stage_order)), showmeans=True)
                
                for i, pc in enumerate(parts['bodies']):
                    pc.set_facecolor(colors[i])
                    pc.set_alpha(0.7)
                
                ax.set_xticks(range(len(stage_order)))
                ax.set_xticklabels(stage_order)
                ax.set_ylabel(ylabel, fontsize=10)
                ax.grid(True, axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('figure_5_distributions.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved: figure_5_distributions.png")
    
    def _save_latex_table(self, df, filename, caption=''):
        """Save DataFrame as LaTeX table"""
        latex_str = df.to_latex(index=False, float_format='%.3f')
        full_latex = f"""\\begin{{table}}[htbp]
\\centering
\\caption{{{caption}}}
\\label{{tab:{filename.replace('.tex', '')}}}
\\small
{latex_str}
\\end{{table}}
"""
        with open(filename, 'w') as f:
            f.write(full_latex)
        print(f"✓ Saved: {filename}")
    
    def run_complete_analysis(self):
        """Execute complete analysis pipeline"""
        print("\n" + "="*80)
        print("PAPER 1: AUTONOMIC-METABOLIC COUPLING ANALYSIS")
        print("="*80)
        
        self.create_intensity_stages()
        self.autonomic_summary()
        self.descriptive_statistics()
        self.correlation_analysis()
        self.intensity_effects()
        
        print("\n=== GENERATING PUBLICATION FIGURES ===\n")
        self.plot_autonomic_metabolic_coupling()
        self.plot_intensity_effects()
        self.plot_complexity_progression()
        self.plot_subject_trajectories()
        self.plot_distribution_by_intensity()
        
        print("\n" + "="*80)
        print("✓ ANALYSIS COMPLETE - Ready for Overleaf")
        print("="*80)


if __name__ == '__main__':
    analysis = Paper1Analysis('merged_data/hrv_metabolic_merged.csv')
    analysis.run_complete_analysis()
