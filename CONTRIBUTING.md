# Contributing to FDTCosmo

Thank you for your interest in contributing to FDTCosmo! This document provides guidelines and information for contributors.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [How to Contribute](#how-to-contribute)
- [Development Setup](#development-setup)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Documentation](#documentation)
- [Submitting Changes](#submitting-changes)

## Code of Conduct

This project adheres to a code of conduct. By participating, you are expected to uphold this code. Please be respectful and constructive in all interactions.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/your-username/FDTCosmo.git
   cd FDTCosmo
   ```
3. Add the upstream repository:
   ```bash
   git remote add upstream https://github.com/original-owner/FDTCosmo.git
   ```
4. Create a branch for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## How to Contribute

### Reporting Bugs

- Use the GitHub issue tracker
- Include a clear description of the bug
- Provide steps to reproduce
- Include relevant error messages and logs
- Specify your environment (OS, Python version, etc.)

### Suggesting Features

- Open an issue with the "enhancement" label
- Describe the feature and its use case
- Explain how it fits within FDT framework

### Contributing Code

1. Check existing issues and PRs to avoid duplicates
2. For significant changes, open an issue first for discussion
3. Follow the coding standards below
4. Include tests for new functionality
5. Update documentation as needed

## Development Setup

### Python Environment

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install development dependencies
pip install pytest black flake8
```

### Julia Environment

```julia
using Pkg
Pkg.develop(path="FDTEffort.jl")
Pkg.test("FDTEffort")
```

## Coding Standards

### Python

- Follow PEP 8 style guide
- Use type hints where appropriate
- Maximum line length: 100 characters
- Use meaningful variable names
- Add docstrings to all public functions

```python
def convert_to_fdt(
    name: str,
    H0: float,
    omega_m: float,
    **kwargs
) -> FDTCosmology:
    """
    Convert standard cosmological parameters to FDT framework.

    Parameters
    ----------
    name : str
        Identifier for the analysis
    H0 : float
        Hubble constant in km/s/Mpc
    omega_m : float
        Total matter density parameter

    Returns
    -------
    FDTCosmology
        FDT-native cosmological parameters
    """
    ...
```

### Julia

- Follow Julia style guide
- Use descriptive function and variable names
- Add docstrings to exported functions

```julia
"""
    α_outside(r, R_s)

Compute the FDT density parameter for r > R_s.

# Arguments
- `r`: Radial distance
- `R_s`: Schwarzschild radius

# Returns
- Density parameter α ∈ [0, 1)
"""
α_outside(r, R_s) = clamp(R_s / r, 0.0, 0.9999)
```

## Testing

### Python Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=.

# Run specific test file
pytest tests/test_fdt_analysis.py
```

### Julia Tests

```julia
using Pkg
Pkg.test("FDTEffort")
```

### Test Guidelines

- Write tests for all new functionality
- Maintain test coverage above 80%
- Use descriptive test names
- Test edge cases and error conditions

## Documentation

### Code Documentation

- All public functions need docstrings
- Include parameter descriptions and return values
- Provide examples where helpful

### LaTeX Reports

- Follow existing report structure
- Include proper citations
- Use consistent notation for FDT quantities

### README Updates

- Keep the README current with new features
- Update usage examples as needed
- Add new datasets to the datasets section

## Submitting Changes

### Pull Request Process

1. Update your branch with upstream changes:
   ```bash
   git fetch upstream
   git rebase upstream/main
   ```

2. Run tests and ensure they pass

3. Format your code:
   ```bash
   black *.py
   flake8 *.py
   ```

4. Commit with a clear message:
   ```bash
   git commit -m "Add feature: description of changes"
   ```

5. Push to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

6. Open a Pull Request on GitHub

### PR Guidelines

- Use a clear, descriptive title
- Reference any related issues
- Describe what the PR does and why
- Include any relevant test results
- Request review from maintainers

### Commit Message Format

```
<type>: <short description>

<longer description if needed>

<references to issues>
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation only
- `style`: Formatting, no code change
- `refactor`: Code restructuring
- `test`: Adding tests
- `chore`: Maintenance tasks

Example:
```
feat: Add eBOSS QSO sample analysis

Extends FDT analysis to quasar samples at z~1.5.
Includes new redshift-dependent alpha calculations.

Closes #42
```

## Questions?

Feel free to open an issue for any questions about contributing. We're happy to help!
