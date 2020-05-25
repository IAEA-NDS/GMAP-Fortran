## About GMAP

GMAP developed by Wolfgang P. Poenitz is a code to obtain evaluations
of cross sections and their uncertainties based on the combined data
from various experiments. The code employs the Bayesian version
of the Generalized Least Squares method and is named after the
mathematicians Gauss, Markov, and Aitken, who all contributed
to the statistical theory around the linear least squares method.
Input to GMAP are experimental data preprocessed by the [DATP] code,
also developed by Poenitz. 

Notably, these codes have been employed for the evaluation
of [neutron cross section standards]. Among other available documents,
there is a [short user guide] for DATP and GMAP and a more detailed
[report] on the technicalities of the codes.

[DATP]: https://www-nds.iaea.org/IAES-NDS/DATP-Fortran
[neutron cross section standards]: https://www-nds.iaea.org/standards/
[short user guide]: https://www-nds.iaea.org/standards/Codes/GMA-User-Guide.pdf
[report]: https://www-nds.iaea.org/standards/Reports/ANL-NDM-139.pdf

## About this repository

The purpose of this repository is to track modifications to GMAP
and tro provide a reference for future developments related to this code.

### Modifications effected so far

- Converted single precision numeric literals to double precision, e.g.,
`0.5` to `0.5D0`, as otherwise the implicit conversion to double precision
in assignments introduces small numerical differences beyond the 
representational capabilities of single precision floats, e.g.,
`x=0.99` may get converted to `x=0.990000090...` if `x` is double precision.

- Modified a few format descriptors for printing table headings.
Those format descriptors which involved strings and spanned several lines
were not interpreted in the same way using the Intel and GNU Fortran
compiler, i.e., printed strings had a different number of whitespace
characters.

Modifications can also be traced by inspecting the 
[commit log](https://github.com/IAEA-NDS/GMAP-Fortran/commits/master).

### Additional `debug` branch

During the translation of this Fortran version to Python, statements have
been introducedd to write additional information to file `debug.out` for
debugging and ensuring the equivalence of the codes. The Fortran code
with these additional statements is stored in the `debug` branch and
those extra code blocks are enclosed by preprocessor directives `
#ifdef VERIFY ... #endif` to facilitate enabling and disabling them.
