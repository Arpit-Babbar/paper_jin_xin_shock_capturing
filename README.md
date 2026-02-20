# Jin-Xin relaxation as a shock capturing method for high order DG/FR schemes

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/10.5281/zenodo.18711723)](https://zenodo.org/doi/10.5281/zenodo.18711723)

This repository contains information and code to reproduce the results
presented in the article
```bibtex
@online{babbar2026jinxin,
  title={Jin-Xin relaxation as a shock capturing method for high order DG/FR schemes},
  author={Artiano, Marco and Babbar, Arpit and Schlottke-Lakemper, Michael and Gassner, Gregor and Ranocha, Hendrik},
  year={2026},
  month={3},
  eprint={?},
  eprinttype={arxiv},
  eprintclass={math.NA}
}
```

If you find these results useful, please cite the article mentioned above.
If you use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{babbar2026jinxinrepro,
  title={Reproducibility repository for
         "Jin-Xin relaxation as a shock capturing method for high order DG/FR schemes"},
  author={Artiano, Marco and Babbar, Arpit and Schlottke-Lakemper, Michael and Gassner, Gregor and Ranocha, Hendrik},
  year={2026},
  howpublished={\url{https://github.com/Arpit-Babbar/paper_jin_xin_shock_capturing}},
  doi={10.5281/zenodo.18711723}
}
```

## Abstract

Jin-Xin relaxation is a method for approximating hyperbolic conservation laws by a system of hyperbolic equations with an $\varepsilon$ dependent stiff source term.
The system formally relaxes to the original conservation law as $\varepsilon \to 0$.
An asymptotic analysis of the Jin-Xin relaxation system shows that it can be seen as a convection-diffusion equation with a diffusion coefficient that depends on the relaxation parameter $\varepsilon$.
This work makes use of this property to use the Jin-Xin relaxation system as a shock capturing method for high order discontinuous Galerkin (DG) or flux reconstruction (FR) schemes.
The idea is to use a smoothness indicator to choose the $\varepsilon$ value in each cell, so that we can use larger $\varepsilon$ values in non-smooth regions to add extra numerical dissipation.
We show how this can be done by using a single stage method by using the compact Runge-Kutta FR method that handles the stiff source term by using IMplicit-EXplicit Runge-Kutta (IMEX-RK) schemes.
Numerical results involving Burgers' equation and the Euler equations are shown to demonstrate the effectiveness of the proposed method.

## Numerical experiments

In order to generate the results from this repository, you need to install [Julia](https://julialang.org).
We recommend using `juliaup`, as detailed in the official website [https://julialang.org](https://julialang.org).

The results have been generated using Julia version 1.10.10, and we recommend installing the same.
Once you have installed Julia, you can clone this repository, enter this directory and start the executable `julia` with the following steps

```shell
git clone https://github.com/Arpit-Babbar/jin_xin_shock_capturing
cd jin_xin_shock_capturing
julia --project=. --threads=auto
```

Then enter the following commands to generate all the data, and plot the 1-D results

```julia
julia> include("run_all_and_plot1d.jl") # Generate all data
```

If you wish to visualize the 2D figures, you need [ParaView](https://www.paraview.org) and its command line version `pvpython`.
Then, in your shell, you can run

```shell
bash plot2d.sh
```

All the figures are now ready and available in the home directory of the repository

## Authors

- [Marco Artiano](https://github.com/MarcoArtiano) (Johannes Gutenberg University Mainz, Germany)
- [Arpit Babbar](https://babbar.dev) (Johannes Gutenberg University Mainz, Germany)
- [Michael Schlottke-Lakemper](https://lakemper.eu) (University of Augsburg, Germany)
- [Gregor Gassner](https://www.mi.uni-koeln.de/NumSim/gregor-gassner) (University of Cologne, Germany)
- [Hendrik Ranocha](https://ranocha.de) (Johannes Gutenberg University Mainz, Germany)

## License

The code in this repository is published under the MIT license, see the `LICENSE` file.

## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
