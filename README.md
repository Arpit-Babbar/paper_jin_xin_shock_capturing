# Jin-Xin relaxation as a shock capturing method for high order DG/FR schemes

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/?)](https://zenodo.org/doi/?)

This repository contains information and code to reproduce the results
presented in the article
```bibtex
@online{babbar2026jinxin,
  title={Jin-Xin relaxation as a shock capturing method for high order DG/FR schemes},
  author={Babbar, Arpit and Artiano, Marco and ...},
  year={2026},
  month={2},
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
  author={Babbar, Arpit and Artiano, Marco and ...},
  year={2026},
  howpublished={\url{https://github.com/Arpit-Babbar/paper_jin_xin_shock_capturing}},
  doi={?}
}
```

## Abstract

TODO

## Numerical experiments

In order to generate the results from this repository, you need to install [Julia](https://julialang.org).
We recommend using `juliaup`, as detailed in the official website [https://julialang.org](https://julialang.org).

The results have been generated using Julia version 1.10.?, and we recommend installing the same.
Once you have installed Julia, you can clone this repository, enter this directory and start the executable `julia` with the following steps

```shell
git clone https://github.com/Arpit-Babbar/jin_xin_shock_capturing
cd jin_xin_shock_capturing
julia --project=. --threads=auto
```

Then enter the following commands to generate all the data, and plot the 1-D results

```julia
julia> include("run_all.jl") # Generate all data
julia> include("plot_all_1d.jl") # Plot 1-D figures
```

If you wish to visualize the 2D figures, you need [ParaView](https://www.paraview.org) and its command line version `pvpython`.
Then, in your shell, you can run

```shell
bash plot_all_2d.sh
```

All the figures are now ready and available in the home directory of the repository

## Authors

- [Arpit Babbar](https://babbar.dev) (Johannes Gutenberg University Mainz, Germany)
- TODO


## License

The code in this repository is published under the MIT license, see the `LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
