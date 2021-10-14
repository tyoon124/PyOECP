# PyOECP
PyOECP is a purely Python-based flexible open-source software for estimating and modeling the complex permittivity obtained from the open-ended coaxial probe (OECP) technique. The software library contains the dielectric spectra of common reference liquids, which can be used to transform the reflection coefficient into the dielectric spectra. Several Python routines that are commonly employed (e.g., SciPy and NumPy) in the field of science and engineering are only required so that the users can alter the software structure depending on their needs. The modeling algorithm exploits the Markov Chain Monte Carlo method for the data regression. The discrete relaxation models can be built by a proper combination of well-known relaxation models. In addition to these models, the electrode polarization, which is a common measurement artifact for interpreting the dielectric spectra, can be incorporated in the modeling algorithm. A continuous relaxation model, which solves the Fredholm integral equation of the first kind (a mathematically ill-posed problem) is also included. This open-source software enables users to freely adjust the physical parameters so that they can obtain physical insight into their materials under test and will be consistently updated for more accurate measurement and interpretation of dielectric spectra in an automated manner.

## Installation
This Package will be uploaded to the developer's webpage and pip.
## Usage
Some examples based on synthetic data and experimental data are available in the Examples folder. The analysis and interpretation is available in the relevant article.
## Citation Information
If you think this package was useful, please consider citing the following article.
@misc{yoon2021pyoecp,
      title={PyOECP: A flexible open-source software library for estimating and modeling the complex permittivity based on the open-ended coaxial probe (OECP) technique}, 
      author={Tae Jun Yoon and Katie A. Maerzke and Robert P. Currier and Alp T. Findikoglu},
      year={2021},
      eprint={2109.14889},
      archivePrefix={arXiv},
      primaryClass={physics.ins-det}
}
## OSS License Information
© 2021. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
others to do so.
