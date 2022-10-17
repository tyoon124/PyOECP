# PyOECP
PyOECP is a purely Python-based flexible open-source software for estimating and modeling the complex permittivity obtained from the open-ended coaxial probe (OECP) technique. The software library contains the dielectric spectra of common reference liquids, which can be used to transform the reflection coefficient into the dielectric spectra. Several Python routines that are commonly employed (e.g., SciPy and NumPy) in the field of science and engineering are only required so that the users can alter the software structure depending on their needs. The modeling algorithm exploits the Markov Chain Monte Carlo method for the data regression. The discrete relaxation models can be built by a proper combination of well-known relaxation models. In addition to these models, the electrode polarization, which is a common measurement artifact for interpreting the dielectric spectra, can be incorporated in the modeling algorithm. A continuous relaxation model, which solves the Fredholm integral equation of the first kind (a mathematically ill-posed problem) is also included. This open-source software enables users to freely adjust the physical parameters so that they can obtain physical insight into their materials under test and will be consistently updated for more accurate measurement and interpretation of dielectric spectra in an automated manner.

## Installation
**(Updated) This package can be installed by downloading the package and installing it. Specifically, the following is a procedure.**

1. Download the package by clicking the Code button and clikcing Download ZIP.
2. Unzip the package. Then, you will create a folder named (for instance) C:/Users/User/Download/PyOECP-main
3. In terminal (or Anaconda Prompt), go to the folder. Then, just type "pip install ." (Or, you can type "python setup.py install")
4. Enjoy! You can load the package in both Jupyter Notebook and your IDE.

## Usage
Some examples based on synthetic data and experimental data are available in the Examples folder. The analysis and interpretation is available in the relevant article.
## Citation Information
If you think this package was useful, please consider citing the following article.
@article{yoon2022pyoecp,
  title={PyOECP: A flexible open-source software library for estimating and modeling the complex permittivity based on the open-ended coaxial probe (OECP) technique},
  author={Yoon, Tae Jun and Maerzke, Katie A and Currier, Robert P and Findikoglu, Alp T},
  journal={Computer Physics Communications},
  pages={108517},
  year={2022},
  publisher={Elsevier}
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
