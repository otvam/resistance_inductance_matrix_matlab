# MATLAB Code for Extracting Resistance/Inductance Matrices

![license - BSD](https://img.shields.io/badge/license-BSD-green)
![language - MATLAB](https://img.shields.io/badge/language-MATLAB-blue)
![category - power electronics](https://img.shields.io/badge/category-power%20electronics-lightgrey)
![status - unmaintained](https://img.shields.io/badge/status-unmaintained-red)

This **MATLAB code** extracts the **resistance/inductance matrix** of a system:
* Different **operating points** are provided (measured or simulated current excitations).
* The associated **losses/energies** are provided (measured or simulated values).
* The corresponding resistance/inductance matrix is **extracted**.

The following the **quadratic form** is linking the aforementionned variables:
* `U = 0.5*I_vec'*Q_mat*I_vec`
* `I_vec` is the vector containing the applied currents
* `Q_mat` is the resistance/inductance matrix
* `U` is the loss/energy value

This code is handling the following cases:
* `n_op` is the number of provided operating points
* `n_var` is the number of independent coefficients for the resistance/inductance matrix
* `n_op<n_var`: under-determined equation system => invalid problem
* `n_op==n_var`: determined equation system => linear equation system
* `n_op>n_var`: over-determined equation system => least-square fit

The condition number of the equation system and the residuum are computed.
This is required to ensure that the resistance/inductance matrix is numerically robust.
Finally, the code checks that the extracted matrix is physical (symmetric and positive definite).

## Compatibility

* Tested with MATLAB R2021a.
* No toolboxes are required.
* Compatibility with GNU Octave not tested but probably easy to achieve.

## Author

**Thomas Guillod** - [GitHub Profile](https://github.com/otvam)

## License

This project is licensed under the **BSD License**, see [LICENSE.md](LICENSE.md).
