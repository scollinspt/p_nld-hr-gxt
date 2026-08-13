# Contributing to Non-Linear Dynamics of Heart Rate Analysis

We welcome contributions, suggestions, and bug reports! Here's how you can help.

## Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/YOUR_USERNAME/p_nld-hr-gxt.git
   cd p_nld-hr-gxt
   ```
3. **Create a virtual environment**:
   ```bash
   python -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   ```
4. **Create a branch** for your feature:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## Development Workflow

### Running the Analysis Pipeline

```bash
# Regenerate HRV metrics
python hrv_analysis.py

# Merge HRV with metabolic data
python merge_data.py

# Generate exploratory plots
python analyze_merged_data.py
```

### Code Style

- Follow PEP 8 conventions
- Use meaningful variable names
- Add docstrings to functions
- Comment complex signal processing logic

### Testing

Please test your changes:
```bash
python -m pytest tests/
```

## Areas for Contribution

### High Priority
- **Exercise stage classification** — Automated detection of rest/submaximal/maximal phases
- **Statistical modeling** — Mixed-effects models for autonomic-metabolic coupling
- **Complexity metrics** — Alternative DFA implementations, wavelet analysis

### Medium Priority
- **Visualization** — Publication-quality figures for manuscript
- **Documentation** — Usage examples, method validation
- **Performance** — Optimization of sliding window calculations

### Lower Priority
- **Alternative HRV metrics** — Entropy variations, time-frequency analysis
- **Data visualization** — Interactive dashboards (Plotly/Dash)
- **CI/CD** — GitHub Actions for automated testing

## Submitting Changes

1. **Commit your changes** with descriptive messages:
   ```bash
   git commit -m "Add feature: autonomic saturation detection"
   ```
2. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```
3. **Create a Pull Request** on GitHub with:
   - Clear description of changes
   - Reference to any related issues
   - Tests demonstrating functionality

## Data Access

Raw data files (CSV, EDF) are excluded from the repository but available:
- Contact the project maintainer for access
- Place data files in project root; they'll be ignored by `.gitignore`

## Reporting Issues

Please report bugs using GitHub Issues with:
- **Title**: Concise description
- **Description**: Steps to reproduce, expected vs actual behavior
- **Environment**: Python version, OS, dependency versions

## Code of Conduct

This project adheres to the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md).

## Questions?

Open an issue or contact the project maintainer.

Thanks for contributing! 🙏
