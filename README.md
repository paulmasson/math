# Math

Nontrivial mathematics in the browser. 100% pure JavaScript. For modern mathematical physicists and their ilk.

All mathematical functions are exposed in the global context for ease of use. While this is contrary to most JavaScript development, it makes life much simpler for the end user. This library is not intended to be sliced and diced by Node users, so please do not complain to me about namespace pollution. The current behavior is intended: mathematics is a global phenomenon.

Complex numbers are supported as dictionaries with real and imaginary attributes. Fully complex functions include exponentials, logarithms, circular and hyperbolic functions, all implemented with built-in JavaScript functions. Also available on the entire complex plane are the gamma and beta functions, implemented with the Lanczos approximation.

Hypergeometric functions <sub>0</sub>*F*<sub>1</sub>, <sub>1</sub>*F*<sub>1</sub> and <sub>2</sub>*F*<sub>1</sub> of real or complex arguments are supported as power series. These make available Bessel *J*, *Y*, *I* and *K* functions on the complex plane, as well as Airy Ai and Bi functions and Hankel functions. Elliptic integrals and elliptic functions of real or complex arguments are available using Carlson symmetric integrals and Jacobi theta function power series. Weierstrass elliptic functions and related support functions are available. Hermite, Laguerre, Legendre and Chebyshev polynomials as well as spherical harmonics are supported.

Here is the complete list of mathematical [functions](https://paulmasson.github.io/math/docs/functions.html) available in this library.

Supported operations include numerical differentiation and integration of complex functions, integration of ordinary differential equations for real values, eigensystems of real symmetric matrices, real or complex root finding, multidimensional real root finding, real cubic splines and Fourier analysis.

This library and [MathCell](https://github.com/paulmasson/mathcell) power the interactive elements of [Analytic Physics](http://analyticphysics.com).

### Usage ###

Download the complete <a href="https://raw.githubusercontent.com/paulmasson/math/master/build/math.js">library</a> and include it in the body of a page:

```
<script src="math.js"></script>
```

See the [documentation](https://paulmasson.github.io/math/) for details of special functions and common mathematical operations.
