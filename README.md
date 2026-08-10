# SP2021: Non-Stationary Signal Decomposition and Time-Frequency Analysis

A MATLAB toolbox for non-stationary signal decomposition, instantaneous frequency estimation, and time-frequency analysis based on adaptive optimization and sparse regularization.

## 🏗️ Method Overview

This repository provides a set of MATLAB implementations for advanced signal processing algorithms, including:

- **Adaptive Signal Separation**: An optimization-based iterative method for decomposing multi-component non-stationary signals into individual components, with automatic estimation of instantaneous amplitude (IA) and instantaneous frequency (IF).
- **Empirical Mode Decomposition (EMD)**: Classic EMD and Ensemble EMD (EEMD) implementations for data-driven signal decomposition.
- **Synchrosqueezing Transform (SST)**: High-resolution time-frequency representation that sharpens the spectrogram by reallocating energy along the frequency axis.
- **Hilbert Spectral Analysis**: Instantaneous frequency and amplitude estimation via the Hilbert transform.
- **Ridge Extraction**: Multi-ridge detection algorithms for tracking frequency trajectories in time-frequency representations.

## 🔧 Requirements

- MATLAB R2018b or later
- Signal Processing Toolbox

## 📁 Project Structure

```text
SP2021/
|-- demo.m                      # Synthetic signal separation demo
|-- demo4real.m                 # Real-world signal analysis demo
|-- extract_one_component.m     # Extract single component from mixture
|-- nspWithTrueT.m              # Non-stationary signal separation (core)
|
|-- emd.m                       # Empirical Mode Decomposition
|-- eemd.m                      # Ensemble EMD
|-- isimf.m                     # Check IMF conditions
|-- ismonotonic.m               # Check monotonicity
|-- extrema.m                   # Find local extrema
|-- getspline.m                 # Spline interpolation for envelope
|-- findpeaks.m                 # Peak detection
|
|-- sst.m                       # Synchrosqueezing Transform
|-- HilbertAnalysis.m           # Hilbert spectral analysis
|-- estimateIA.m                # Instantaneous amplitude estimation
|
|-- brevridge.m                 # Ridge extraction (single)
|-- brevridge_mult.m            # Ridge extraction (multiple)
|-- wvdc.m                      # Wavelet-transform based decomposition
|
|-- TKoperator.m                # Teager-Kaiser energy operator
|-- Exterpolation.m             # Signal extrapolation
|-- filter_mean.m               # Mean filtering
|-- sgolayfiltDesign.m          # Savitzky-Golay filter design
|
|-- updateP.m                   # IF parameter update
|-- updateQ.m                   # IA parameter update
|-- updatePQ.m                  # Joint IF/IA parameter update
|-- snr_compute.m               # SNR computation
|
|-- zai.mat                     # Real-world test data
```

## 🚀 Quick Start

### Synthetic Signal Separation

```matlab
>> demo
```

Demonstrates separation of two overlapping non-stationary signals with close frequency components. Plots the input signal, recovered component with instantaneous amplitude, and estimated instantaneous frequency.

### Real Data Analysis

```matlab
>> demo4real
```

Loads real signal data from `zai.mat` and extracts up to 5 components using iterative adaptive decomposition. Each component is plotted for visual inspection.

## 📦 Core Functions

### Signal Separation

```matlab
[u, v, p, q] = extract_one_component(s, wd, n1, n2, beta, eps);
```

| Parameter | Description |
|-----------|-------------|
| `s` | Input multi-component signal |
| `wd` | Boundary padding width |
| `n1`, `n2` | Filter parameters |
| `beta` | Regularization weight |
| `eps` | Convergence threshold |
| `u` | Residual (remaining components) |
| `v` | Extracted single component |
| `p` | Estimated squared IF reciprocal |
| `q` | Estimated IA ratio |

### EMD

```matlab
imf = emd(s);
```

Decompose signal `s` into Intrinsic Mode Functions (IMFs).

### SST

```matlab
[Tf, freq] = sst(s, dt);
```

Compute synchrosqueezed time-frequency transform of signal `s` with sampling interval `dt`.

## 📝 Citation

```bibtex
@article{hu2021accurate,
  title={Accurate AM-FM signal demodulation and separation using nonparametric regularization method},
  author={Hu, Xiyuan and Peng, Silong and Guo, Baokui and Xu, Pengcheng},
  journal={Signal Processing},
  volume={186},
  pages={108131},
  year={2021},
  publisher={Elsevier}
}
```