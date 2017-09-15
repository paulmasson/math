# Math

Nontrivial mathematics in the browser. 100% pure JavaScript. For modern mathematical physicists and their ilk.

All mathematical functions are exposed in the global context for ease of use. While this is contrary to most JavaScript development, it makes life much simpler for the end user. This library is not intended to be sliced and diced by Node users, so please do not complain to me about namespace pollution. The current behavior is intended: mathematics is a global phenomenon.

Complex numbers are supported as dictionaries with real and imaginary attributes. Fully complex functions include exponentials, logarithms, circular and hyperbolic functions, all implemented with built-in JavaScript functions. Also available in the entire complex plane are the gamma and beta functions, implemented with the Lanczos approximation.

Hypergeometric functions <sub>0</sub>*F*<sub>1</sub>, <sub>1</sub>*F*<sub>1</sub> and <sub>2</sub>*F*<sub>1</sub> of a real argument are currently supported as power series. These make available Bessel *J*, *Y*, *I* and *K* functions of a real argument, as well as Airy Ai and Bi functions. Elliptic integrals and elliptic functions of a real argument are also available using Carlson symmetric integrals and Jacobi theta function power series.

Integration of functions and ordinary differential equations is available for real values, as well as eigensystems of real symmetric matrices, root finding and cubic splines.

This library and [MathCell](https://github.com/paulmasson/mathcell) power the interactive elements of [Analytic Physics](http://analyticphysics.com).

### Usage ###

Download the complete <a href="https://raw.githubusercontent.com/paulmasson/math/master/build/math.js">library</a> and include it in the body of a page:

```
<script src="math.js"></script>
```

See the [documentation](https://paulmasson.github.io/math/) for details of special functions and common mathematical operations.
