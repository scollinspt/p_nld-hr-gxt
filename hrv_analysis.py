"""
Non-Linear Dynamics of Heart Rate during Graded Exercise
HRV Analysis Pipeline - Python Implementation

Purpose:
    1. Load and validate RR intervals from Kubios export
    2. Implement HRV metric calculations
    3. Validate against existing Kubios HRV values
    4. Regenerate HRV with 10-second window advancement (vs. original 30-second)
    5. Export aligned HRV + metabolic data for statistical analysis

Author: Sean M. Collins
Institution: Physiological Dynamics Lab, Plymouth State University
Date: 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import signal, stats
from scipy.fft import fft, fftfreq
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# Configuration
# ============================================================================

PROJECT_DIR = Path("/Users/collins/Projects/p_nld-hr-gxt-master")
RR_DIR = PROJECT_DIR / "rr"
DATA_DIR = PROJECT_DIR
OUTPUT_DIR = PROJECT_DIR / "python_analysis"
OUTPUT_DIR.mkdir(exist_ok=True)

# Subject ID mapping
ID_MAPPING = pd.read_csv(DATA_DIR / "Old-to-New-IDs.csv")
ID_MAPPING_DICT = dict(zip(ID_MAPPING['id-old'], ID_MAPPING['id-new']))

# HRV Analysis Parameters
WINDOW_DURATION = 300.0  # 5 minutes in seconds
ORIGINAL_ADVANCEMENT = 30.0  # Original Kubios advancement
NEW_ADVANCEMENT = 10.0  # New advancement for tighter temporal resolution
MIN_WINDOW_RR_COUNT = 200  # Minimum RR intervals for valid HRV window
HRV_FREQ_BANDS = {
    'VLF': (0.0033, 0.04),      # Very Low Frequency (0.0033-0.04 Hz)
    'LF': (0.04, 0.15),          # Low Frequency (0.04-0.15 Hz)
    'HF': (0.15, 0.4)            # High Frequency (0.15-0.4 Hz)
}

# ============================================================================
# HRV Metric Calculations
# ============================================================================

class HRVAnalyzer:
    """Calculate HRV metrics from RR interval time series."""
    
    def __init__(self):
        """Initialize HRV analyzer."""
        pass
    
    @staticmethod
    def rr_to_intervals(rr_times):
        """
        Convert RR peak timestamps to RR interval durations.
        
        Parameters
        ----------
        rr_times : array-like
            R-peak times in seconds
            
        Returns
        -------
        rr_intervals : ndarray
            RR intervals in milliseconds
        """
        if len(rr_times) < 2:
            return np.array([])
        
        # Calculate intervals between consecutive peaks and convert to ms
        intervals_sec = np.diff(rr_times)
        intervals_ms = intervals_sec * 1000.0
        
        return intervals_ms
    
    @staticmethod
    def remove_ectopic_beats(rr_intervals, threshold_pct=30):
        """
        Remove likely ectopic beats using threshold on interval duration.
        
        Parameters
        ----------
        rr_intervals : ndarray
            RR intervals in milliseconds
        threshold_pct : float
            Percentage threshold for deviation from median
            
        Returns
        -------
        cleaned_intervals : ndarray
            RR intervals with ectopic beats removed
        """
        if len(rr_intervals) < 3:
            return rr_intervals
        
        median_rr = np.median(rr_intervals)
        lower_bound = median_rr * (1 - threshold_pct/100)
        upper_bound = median_rr * (1 + threshold_pct/100)
        
        mask = (rr_intervals >= lower_bound) & (rr_intervals <= upper_bound)
        return rr_intervals[mask]
    
    @staticmethod
    def time_domain_metrics(rr_intervals):
        """
        Calculate time-domain HRV metrics.
        
        Parameters
        ----------
        rr_intervals : ndarray
            RR intervals in milliseconds
            
        Returns
        -------
        dict : Time-domain metrics
        """
        if len(rr_intervals) < 2:
            return {m: np.nan for m in ['meanrr', 'sdnn', 'rmssd', 'pnn50', 'nn50']}
        
        rr_diff = np.diff(rr_intervals)
        
        metrics = {
            'meanrr': np.mean(rr_intervals),
            'sdnn': np.std(rr_intervals, ddof=1),
            'rmssd': np.sqrt(np.mean(rr_diff**2)),
            'nn50': np.sum(np.abs(rr_diff) > 50),
            'pnn50': 100.0 * np.sum(np.abs(rr_diff) > 50) / len(rr_diff) if len(rr_diff) > 0 else 0,
            'meanhr': 60000.0 / np.mean(rr_intervals),
        }
        
        return metrics
    
    @staticmethod
    def frequency_domain_metrics(rr_intervals, sampling_rate=4.0):
        """
        Calculate frequency-domain HRV metrics using FFT.
        
        Parameters
        ----------
        rr_intervals : ndarray
            RR intervals in milliseconds
        sampling_rate : float
            Interpolation sampling rate in Hz (default: 4 Hz)
            
        Returns
        -------
        dict : Frequency-domain metrics
        """
        if len(rr_intervals) < 2:
            return {m: np.nan for m in ['vlfpowfft', 'lfpowfft', 'hfpowfft', 'lfhffft']}
        
        # Detrend and normalize RR intervals
        rr_norm = signal.detrend(rr_intervals - np.mean(rr_intervals))
        
        # Cubic spline interpolation to uniform sampling
        time = np.cumsum(rr_intervals) / 1000.0  # Convert to seconds
        time_uniform = np.arange(0, time[-1], 1.0/sampling_rate)
        
        try:
            rr_interp = np.interp(time_uniform, time, rr_norm)
        except:
            return {m: np.nan for m in ['vlfpowfft', 'lfpowfft', 'hfpowfft', 'lfhffft']}
        
        # Apply Hanning window and compute FFT
        window = signal.windows.hann(len(rr_interp))
        fft_vals = np.abs(fft(rr_interp * window))
        freqs = fftfreq(len(rr_interp), 1.0/sampling_rate)
        
        # Keep only positive frequencies
        pos_mask = freqs >= 0
        freqs = freqs[pos_mask]
        power = (fft_vals[pos_mask] ** 2) / len(rr_interp)
        
        # Band power integration
        metrics = {}
        for band_name, (f_low, f_high) in HRV_FREQ_BANDS.items():
            band_mask = (freqs >= f_low) & (freqs <= f_high)
            band_power = np.sum(power[band_mask])
            metrics[f'{band_name.lower()}powfft'] = band_power
        
        total_power = metrics['vlfpowfft'] + metrics['lfpowfft'] + metrics['hfpowfft']
        metrics['totpowfft'] = total_power
        metrics['lfhffft'] = metrics['lfpowfft'] / metrics['hfpowfft'] if metrics['hfpowfft'] > 0 else np.nan
        
        return metrics
    
    @staticmethod
    def sample_entropy(rr_intervals, m=2, r=None):
        """
        Calculate sample entropy (complexity measure).
        
        Parameters
        ----------
        rr_intervals : ndarray
            RR intervals in milliseconds
        m : int
            Embedding dimension (default: 2)
        r : float
            Tolerance (default: 0.2 * std(RR))
            
        Returns
        -------
        float : Sample entropy value
        """
        if len(rr_intervals) < m + 1:
            return np.nan
        
        if r is None:
            r = 0.2 * np.std(rr_intervals, ddof=1)
        
        # Normalize
        rr_norm = (rr_intervals - np.mean(rr_intervals)) / np.std(rr_intervals, ddof=1)
        
        # Count template matches
        def count_patterns(x, m, r):
            n = len(x)
            count = 0
            for i in range(n - m):
                for j in range(i + 1, n - m + 1):
                    if np.max(np.abs(x[i:i+m] - x[j:j+m])) <= r:
                        count += 1
            return count
        
        phi_m = count_patterns(rr_norm, m, r)
        phi_m1 = count_patterns(rr_norm, m + 1, r)
        
        if phi_m1 == 0:
            return np.nan
        
        sampen = -np.log(phi_m1 / phi_m)
        return sampen
    
    @staticmethod
    def detrended_fluctuation_analysis(rr_intervals, min_scale=10, max_scale=500):
        """
        Calculate Detrended Fluctuation Analysis (DFA) exponents.
        
        Parameters
        ----------
        rr_intervals : ndarray
            RR intervals in milliseconds
        min_scale : int
            Minimum scale (samples)
        max_scale : int
            Maximum scale (samples)
            
        Returns
        -------
        dict : DFA exponents (dfa1, dfa2)
        """
        if len(rr_intervals) < max_scale * 2:
            return {'dfa1': np.nan, 'dfa2': np.nan}
        
        # Normalize RR intervals
        rr_norm = rr_intervals - np.mean(rr_intervals)
        
        # Cumulative sum (integration)
        y = np.cumsum(rr_norm)
        
        scales = np.logspace(np.log10(min_scale), np.log10(max_scale), 30).astype(int)
        scales = np.unique(scales)
        fluctuations = []
        
        for scale in scales:
            # Forward direction
            n_full = len(y) // scale
            y_forward = y[:n_full * scale].reshape(n_full, scale)
            
            # Fit polynomial trend
            fit_forward = np.polyfit(np.arange(scale), y_forward[0], 1)
            
            # Calculate fluctuation
            fluct = 0
            for i in range(n_full):
                trend = np.polyval(fit_forward, np.arange(scale))
                fluct += np.mean((y_forward[i] - trend) ** 2)
            
            fluctuations.append(np.sqrt(fluct / n_full))
        
        # Fit log-log regression
        log_scales = np.log10(scales)
        log_fluct = np.log10(fluctuations)
        
        # Split at middle scale for DFA1 and DFA2
        mid_idx = len(scales) // 2
        
        slope1 = np.polyfit(log_scales[:mid_idx], log_fluct[:mid_idx], 1)[0]
        slope2 = np.polyfit(log_scales[mid_idx:], log_fluct[mid_idx:], 1)[0]
        
        return {'dfa1': slope1, 'dfa2': slope2}
    
    def calculate_all_metrics(self, rr_times):
        """
        Calculate all HRV metrics for a given RR interval sequence.
        
        Parameters
        ----------
        rr_times : array-like
            R-peak times in seconds
            
        Returns
        -------
        dict : All HRV metrics
        """
        # Convert timestamps to intervals
        rr_intervals = self.rr_to_intervals(rr_times)
        
        if len(rr_intervals) < MIN_WINDOW_RR_COUNT:
            return None
        
        # Remove ectopic beats
        rr_clean = self.remove_ectopic_beats(rr_intervals)
        
        if len(rr_clean) < MIN_WINDOW_RR_COUNT:
            return None
        
        # Calculate all metrics
        metrics = {}
        
        # Time domain
        metrics.update(self.time_domain_metrics(rr_clean))
        
        # Frequency domain
        metrics.update(self.frequency_domain_metrics(rr_clean))
        
        # Complexity measures
        metrics['sampen'] = self.sample_entropy(rr_clean)
        dfa_metrics = self.detrended_fluctuation_analysis(rr_clean)
        metrics.update(dfa_metrics)
        
        return metrics


# ============================================================================
# Data Loading and Processing
# ============================================================================

def load_rr_intervals(subject_id_old):
    """
    Load RR intervals from Kubios export file.
    
    Parameters
    ----------
    subject_id_old : int
        Original subject ID (e.g., 130, 146)
        
    Returns
    -------
    array : RR peak times in seconds
    """
    rr_file = RR_DIR / f"{subject_id_old}.txt"
    
    if not rr_file.exists():
        raise FileNotFoundError(f"RR file not found: {rr_file}")
    
    rr_times = np.loadtxt(rr_file)
    return rr_times


def create_windows(rr_times, window_duration=WINDOW_DURATION, advancement=ORIGINAL_ADVANCEMENT):
    """
    Create sliding windows of RR intervals.
    
    Parameters
    ----------
    rr_times : array
        RR peak times in seconds
    window_duration : float
        Window length in seconds
    advancement : float
        Window advancement in seconds
        
    Returns
    -------
    list : List of (window_start, window_end, rr_subset) tuples
    """
    windows = []
    test_duration = rr_times[-1]
    
    window_start = 0
    window_count = 0
    
    while window_start + window_duration <= test_duration:
        window_end = window_start + window_duration
        
        # Extract RR intervals within this window
        mask = (rr_times >= window_start) & (rr_times <= window_end)
        rr_subset = rr_times[mask]
        
        if len(rr_subset) >= MIN_WINDOW_RR_COUNT:
            windows.append({
                'window_num': window_count,
                'window_start': window_start,
                'window_end': window_end,
                'rr_times': rr_subset
            })
            window_count += 1
        
        window_start += advancement
    
    return windows


# ============================================================================
# Main Analysis
# ============================================================================

def analyze_subject(subject_id_old, advancement=ORIGINAL_ADVANCEMENT):
    """
    Analyze single subject: calculate HRV for all windows.
    
    Parameters
    ----------
    subject_id_old : int
        Original subject ID
    advancement : float
        Window advancement in seconds
        
    Returns
    -------
    DataFrame : HRV metrics by window
    """
    # Load RR times
    rr_times = load_rr_intervals(subject_id_old)
    
    # Create windows
    windows = create_windows(rr_times, window_duration=WINDOW_DURATION, advancement=advancement)
    
    # Calculate HRV for each window
    analyzer = HRVAnalyzer()
    results = []
    
    for window_info in windows:
        metrics = analyzer.calculate_all_metrics(window_info['rr_times'])
        
        if metrics is not None:
            metrics['subject_id_old'] = subject_id_old
            metrics['subject_id_new'] = ID_MAPPING_DICT.get(subject_id_old)
            metrics['window_num'] = window_info['window_num']
            metrics['window_start'] = window_info['window_start']
            metrics['window_end'] = window_info['window_end']
            metrics['rr_count'] = len(window_info['rr_times'])
            
            results.append(metrics)
    
    return pd.DataFrame(results)


def validate_against_kubios(sample_num=1):
    """
    Validate calculated HRV against existing Kubios values.
    
    Parameters
    ----------
    sample_num : int
        Which sample (window) to compare (1-indexed)
        
    Returns
    -------
    DataFrame : Comparison of our values vs. Kubios
    """
    # Load our calculated HRV
    subject_id_old = 130  # Use first subject for validation
    subject_id_new = 1
    
    df_our_hrv = analyze_subject(subject_id_old, advancement=ORIGINAL_ADVANCEMENT)
    
    # Load Kubios HRV
    df_kubios = pd.read_csv(DATA_DIR / "exHRVdata-SelectedVariables.csv")
    df_kubios_subj = df_kubios[df_kubios['Subject'] == subject_id_new]
    
    # Compare specific window
    if sample_num <= len(df_our_hrv):
        our_row = df_our_hrv.iloc[sample_num - 1]
        kubios_row = df_kubios_subj[df_kubios_subj['Sample'] == sample_num].iloc[0]
        
        # Metrics to compare
        compare_metrics = ['meanrr', 'sdnn', 'rmssd', 'pnn50', 'meanhr']
        
        comparison = pd.DataFrame({
            'Metric': compare_metrics,
            'Our Value': [our_row.get(m, np.nan) for m in compare_metrics],
            'Kubios Value': [kubios_row.get(m, np.nan) for m in compare_metrics],
        })
        
        comparison['Difference'] = comparison['Our Value'] - comparison['Kubios Value']
        comparison['Pct Diff'] = 100 * comparison['Difference'] / comparison['Kubios Value'].abs()
        
        return comparison
    else:
        return None


def process_all_subjects_original():
    """
    Process all subjects with original 30-second advancement.
    Export for comparison with existing Kubios data.
    """
    all_results = []
    
    for subject_id_old in sorted(ID_MAPPING['id-old'].unique()):
        try:
            print(f"Processing subject {subject_id_old}...", end=" ")
            df = analyze_subject(subject_id_old, advancement=ORIGINAL_ADVANCEMENT)
            all_results.append(df)
            print(f"✓ ({len(df)} windows)")
        except Exception as e:
            print(f"✗ Error: {e}")
    
    df_all = pd.concat(all_results, ignore_index=True)
    
    # Save
    output_file = OUTPUT_DIR / "hrv_python_original_advancement.csv"
    df_all.to_csv(output_file, index=False)
    print(f"\n✓ Saved to {output_file}")
    
    return df_all


def process_all_subjects_new():
    """
    Process all subjects with NEW 10-second advancement.
    This is the enhanced dataset for metabolic coupling analysis.
    """
    all_results = []
    
    for subject_id_old in sorted(ID_MAPPING['id-old'].unique()):
        try:
            print(f"Processing subject {subject_id_old}...", end=" ")
            df = analyze_subject(subject_id_old, advancement=NEW_ADVANCEMENT)
            all_results.append(df)
            print(f"✓ ({len(df)} windows)")
        except Exception as e:
            print(f"✗ Error: {e}")
    
    df_all = pd.concat(all_results, ignore_index=True)
    
    # Save
    output_file = OUTPUT_DIR / "hrv_python_new_advancement.csv"
    df_all.to_csv(output_file, index=False)
    print(f"\n✓ Saved to {output_file}")
    
    return df_all


def plot_validation(sample_num=1):
    """
    Create visualization of validation comparison.
    """
    comparison = validate_against_kubios(sample_num=sample_num)
    
    if comparison is None:
        print("Sample not found")
        return
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x = np.arange(len(comparison))
    width = 0.35
    
    ax.bar(x - width/2, comparison['Our Value'], width, label='Python', alpha=0.8)
    ax.bar(x + width/2, comparison['Kubios Value'], width, label='Kubios', alpha=0.8)
    
    ax.set_ylabel('Value')
    ax.set_title(f'HRV Validation: Subject 1, Sample {sample_num}')
    ax.set_xticks(x)
    ax.set_xticklabels(comparison['Metric'])
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    output_file = OUTPUT_DIR / f"validation_sample{sample_num}.png"
    plt.savefig(output_file, dpi=300)
    print(f"✓ Saved validation plot to {output_file}")
    plt.close()


# ============================================================================
# Main Execution
# ============================================================================

if __name__ == "__main__":
    print("="*70)
    print("HRV Analysis Pipeline - Validation & Regeneration")
    print("="*70)
    
    # Step 1: Validation on first subject/sample
    print("\n[STEP 1] Validating calculation against Kubios...")
    print("-"*70)
    
    try:
        comparison = validate_against_kubios(sample_num=1)
        print("\nComparison (Sample 1 of Subject 1):")
        print(comparison.to_string(index=False))
        plot_validation(sample_num=1)
    except Exception as e:
        print(f"Validation error: {e}")
    
    # Step 2: Process all subjects with original advancement
    print("\n[STEP 2] Processing all subjects (original 30-sec advancement)...")
    print("-"*70)
    
    try:
        df_original = process_all_subjects_original()
        print(f"\nTotal windows processed: {len(df_original)}")
        print(f"Subjects: {df_original['subject_id_new'].nunique()}")
    except Exception as e:
        print(f"Processing error: {e}")
    
    # Step 3: Process all subjects with NEW advancement
    print("\n[STEP 3] Processing all subjects (NEW 10-sec advancement)...")
    print("-"*70)
    
    try:
        df_new = process_all_subjects_new()
        print(f"\nTotal windows processed: {len(df_new)}")
        print(f"Subjects: {df_new['subject_id_new'].nunique()}")
        print(f"Windows per subject (avg): {len(df_new) / df_new['subject_id_new'].nunique():.1f}")
    except Exception as e:
        print(f"Processing error: {e}")
    
    print("\n" + "="*70)
    print("Analysis complete!")
    print("="*70)
