# MATLAB Code for Extracting Resistance/Inductance Matrices

This **MATLAB code** extracts the **resistance/inductance matrix** of a system:
* Different **operating points** are provided (current excitations).
* The associated **losses/energies** are provided.
* The corresponding resistance/inductance matrix is extracted.

The following the **quadratic form** is linking the different variables:
* `U = 0.5*I'*Q_mat*I`
* `I` is the current excitation vector
* `U` is the loss/energy value
* `Q_mat` is the resistance/inductance matrix

This code is handling the following cases:
* `n_op` is the number of provided operating conditions
* `n_var` is the number of inpedent coefficient of the resistance/inductance matrix
* `n_op<n_var`: under-determined equation system => invalid problem
* `n_op==n_var`: determined equation system => linear equation system
* `n_op>n_var`: over-determined equation system => least-square fit

The condition number of the equation system and the residuum are provided.
Finally, the code ensures that the extracted matrix is physical (symmetric and positive definite).

## Compatibility

* Tested with MATLAB R2021a.
* No toolboxes are required.
* Compatibility with GNU Octave not tested but probably easy to achieve.

## Author

**Thomas Guillod** - [GitHub Profile](https://github.com/otvam)

## License

This project is licensed under the **BSD License**, see [LICENSE.md](LICENSE.md).
