# The Fortran Asteroid Retrieval Toolkit

![License](https://img.shields.io/github/license/tylerjackoliver/ACROBAT)


<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/tylerjackoliver/Fortran_Asteroid_Retrieval_Toolkit">
    <img src="images/asteroid.png" alt="Logo" width="80" height="80">
  </a>
  <h3 align="center">FART: Fortran Asteroid Retrieval Toolkit</h3>
</p>

<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Project](#about-the-project)
  * [Built With](#built-with)
* [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#building)
* [Contact](#contact)



<!-- ABOUT THE PROJECT -->
## About The Project

As part of ongoing research into Asteroid Retrieval in the CR3BP, the Fortran Asteroid Retrieval Toolkit provides the basic computational tools to find Easily Retrievable Objects (EROs) for a given database of Near-Earth Objects, as in (Tyler, J., Wittig, A., 2020 \[Preprint\]). Several tools, written largely in Fortran, are included:

* The Fortran Asteroid Prefilter: Takes a given NEO database and filters it according to a Hohmann transfer approximation to obtain _retrieval candidates_.
* The Fortran Asteroid Retrieval Tool: Takes a given list of _retrieval candidates_ and optimises their transfers for minimum transfer velocity and whether the object can be considered an ERO.
* A script, written in Expect, to automatically extract epehemerides from the JPL HORIZONS system.
* MIDACO Wrappers, designed to interface C versions of MIDACO with Fortran objective function codes.

All of the codes are written in Fortran, with small interfaces to C++, and parallelised using the OpenMPI and MPI parallel programming paradigms. The problems are embarassingly parallel and scale nearly linearly with core count. **Basic Python versions of the Asteroid Retrieval _Tool_ are available on request.**

### Built With

* [BOOST](https://www.boost.org/)
* [OpenMP](https://www.openmp.org/)
* [Message-Passing Interface](https://www.mpich.org/) MPICH is officially supported, although any MPI implementation that supports Fortran will do.
* [CMake](https://cmake.org/)
* [SPICE](https://naif.jpl.nasa.gov/naif/)
* [Fortran Astrodynamics Toolkit](https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit)
* [Fortran B-Spline Library](https://github.com/jacobwilliams/bspline-fortran)
* [MIDACO](https://midaco-solver.com) - this is a paid optimiser for more than 4 design variables (as is the case here for the retrieval tool.) Either a license must be purchased, or the Tool must be rewritten to work for your optimiser. This is explained in detail in the local folder for the tool.

Both the Astrodynamics toolkit and B-Spline library are kindly provided by [Jacob Williams](https://github.com/jacobwilliams); his codes also depend on [FoBiS.py](https://github.com/szaghi/FoBiS), an automated build library for Fortran written in Python.

`gfortran`, `gcc` and `g++` are the officially-supported compilers. In theory, the Intel Compilers should also work, but use at your own risk.

<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

The dependencies given above are required prior to build. The program is currently built and tested with the GNU compiler systems, wrapped by the MPICH library.

### Building

1. Clone the repo
```sh
git clone https://github.com/tylerjackoliver/Fortran-Asteroid-Retrieval-Toolkit.git
```
2. Navigate to the tool you wish to use.

3. For Fortran codes, install the dependencies and run the CMake wrapper
```sh
sh build.sh
```
4. For Expect codes, use as per your system instructions.

<!-- CONTACT -->
## Contact

Jack Tyler - [@tylerjackoliver](https://twitter.com/tylerjackoliver) - jack.tyler@soton.ac.uk

Project Link: [https://github.com/tylerjackoliver/Fortran-Asteroid-Retrieval-Toolkit](https://github.com/tylerjackoliver/Fortran-Asteroid-Retrieval-Toolkit)
