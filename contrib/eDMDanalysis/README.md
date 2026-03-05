## Random Fourier features in eDMD for analyzing low beta, theta neural dynamics of working memory
David Chavez Huerta
## Introduction

This software provides a MATLAB implementation of the Extended Dynamic Mode Decomposition (eDMD) algorithm fully integrated with the FieldTrip software toolbox for MEG, EEG and iEEG analysis. The presented software offers several sets of observable functions: Random Fourier features (RFF), Hermite polynomials, polynomial basis, and the identity observable.
Through the mode decomposition inherent to eDMD, the user has the option to generate a reconstruction of the dynamics in the data from the initial state. 
Frequency content is deduced from the eDMD modes. Power associated to these frequencies is also calculated. By binning the frequency content into relevant neural bands (alpha, beta, gamma...), peak frequencies per band are calculated, as well as total power per band.  

The package integrates seamlessly with FieldTrip, enabling users to incorporate nonlinear methods into existing EEG workflows. The software also allows for the implementation of custom dictionaries, making it extendable to task-specific and model-specific nonlinear representations. The software is meant to support researchers studying neural dynamics, cognitive processes, and data-driven modeling in neuroscience.

## Statement of need

Analyzing human EEG and MEG datasets often requires methods capable of capturing nonlinear, time-varying dynamics that traditional linear models cannot fully represent. Extended Dynamic Mode Decomposition (eDMD) provides a framework based on Koopman theory for learning interpretable models in a feature space, enabling the extraction of coherent temporal patterns from high-dimensional, noisy neural recordings. However, existing implementations generally require substantial customization, miss sufficient documentation, or lack integration with standard neuroscience toolboxes.

This software addresses this gap by providing a neuroscience-focused eDMD pipeline built around nonlinear dictionaries, enabling efficient approximation of nonlinear dynamics in high-dimensional EEG data. The eDMD implementation includes parameter choices tailored to the characteristics of EEG data, while still allowing users to adjust them easily for different analysis goals or needs
Concise documentation explains the effect of each parameter, supporting operation and customization. This makes the software immediately usable for cognitive-neuroscience data analysis.

The package is fully integrated to FieldTrip, allowing researchers to incorporate eDMD-based modeling directly into well established pre-processing and analysis workflows. It supports the substitution or extension of dictionaries, facilitating exploration of task-specific nonlinear observables or alternative basis functions. By providing an integrated, accessible, and extensible implementation, the software fills a methodological need for researchers seeking to integrate Koopman-based techniques into their own analyses. The software is currently used to identify frontal midline theta and high-beta activity correlated with working memory encoding in a EEG dataset via the RFF dictionary, as it shows an excellent compromise of reconstruction capabilities and computational cost.

## Installation

### Requirements

* MATLAB (R2021a or later recommended)

* FieldTrip toolbox (recent stable release recommended)

### Install FieldTrip

* Download FieldTrip from the official FieldTrip website: https://www.fieldtriptoolbox.org/

* Add it to your MATLAB path and initalize:

`ft_defaults`

* To verify the installation:

`ft_version`

### Install This Package

* Clone the repository using Git:

`git clone https://github.com/yourusername/your-repo-name.git`

* Or download the ZIP file and extract it

* Add the repository folder to your MATLAB path:

`addpath(genpath('path_to_repo'))`

### Reproducibility

* For deterministic behavior when using random Fourier features:

`cfg.seed = 1;`

* Tested Configuration

MATLAB R2023b

FieldTrip (2024 stable release)

* Users are encouraged to report compatibility issues via the project repository.

NOTE:
The computational complexity of eDMD depends on dictionary size, data sampling, number of channels, and stacking depth. Users should adjust cfg.D, cfg.poly_degree, and cfg.nstacks carefully for large trials.

## Example Usage

The function ft_edmd2 performs Extended Dynamic Mode Decomposition (eDMD) on FieldTrip raw data structures.
The input data must be a valid FieldTrip raw structure containing time series data in datain.trial.

### Minimal Example
```
cfg = [];
cfg.dictionary = 'identity';
cfg.reconstruct = 'no';

dataout = ft_edmd2(cfg, datain);
```
This runs classical DMD on all trials.

### Using Random Fourier Features (default)
```
cfg = [];
cfg.dictionary = 'rff';
cfg.D = 900;
cfg.gamma = 4;
cfg.seed = 1;

dataout = ft_edmd2(cfg, datain);
```
This applies an RFF dictionary with 900 random features.

### Using a Polynomial Dictionary
```
cfg = [];
cfg.dictionary = 'poly';
cfg.poly_degree = 3;

dataout = ft_edmd2(cfg, datain);
```
This constructs a monomial polynomial dictionary up to degree 3.

### Using a Hermite Dictionary
```
cfg = [];
cfg.dictionary = 'hermite';
cfg.hermite_degree = 3;

dataout = ft_edmd2(cfg, datain);
```
This applies probabilists’ Hermite polynomials independently to each channel.

### State Reconstruction

To reconstruct the time-domain signal from Koopman modes:
```
cfg = [];
cfg.reconstruct = 'yes';

dataout = ft_edmd2(cfg, datain);
```
The reconstructed signal will be available in:
```
dataout.xrecon
```
### Output Structure

The returned structure contains:


`dataout.peakfeatures` – dominant frequency per bin and trial


`dataout.sumfeatures` – total power per frequency bin and trial


`dataout.modefreqs` – unique Koopman mode frequencies


`dataout.modepowers` – corresponding mode powers


`dataout.xrecon` – reconstructed signals (if enabled)


Each cell corresponds to one trial in the input data.

### API

The main entry point of the package is:

`ft_edmd2(cfg, datain)`


Inputs: refer to FieldTrip documentation for a complete description.

`cfg` – configuration structure controlling the decomposition

`datain` – FieldTrip raw data structure


Key Configuration Options


* dictionary – 'identity', 'rff', 'poly', or 'hermite'

* reconstruct – 'yes' or 'no'

* D – number of random Fourier features

* gamma – RFF scaling parameter

* poly_degree – polynomial degree

* hermite_degree – Hermite polynomial degree

* nstacks – Hankel stacking depth

* MA – moving average window size

* freqEdges – frequency bin edges


Output: Returns a Filedtrip output structure containing:


`peakfeatures`

`sumfeatures`

`modefreqs`

`modepowers`

`xrecon` (optional)

### Community Guidelines

We welcome contributions, feedback, and questions from the community. Please be kind and constructive in all interactions

Contributing

* Contributions are welcome through pull requests

* Please fork the repository and create a new branch for your changes

* Ensure that your code is clear, documented, and consistent with the existing style

* Provide a clear description of the changes in your pull request

Reporting Issues

* Bugs and issues should be reported via the repository issue tracker

* Please include a clear description of the problem, relevant code or configuration
  
Support

* For questions or usage help, please open an issue in the repository. We will do our best to respond, but response times may vary


### License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0).
See the LICENSE file for details.