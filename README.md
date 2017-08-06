# math.js

Nontrivial mathematics in the browser. 100% pure JavaScript.

Complex numbers are supported as a dictionary with `re` and `im` attributes. Fully complex functions include exponentials, logarithms and circular and hyperbolic trigonometric functions, all implemented with builtin JavaScript functions. Also available in the entire complex plane are the gamma and beta functions, implemented with the Lanczos approximation.

Hypergeometric functions 0F1, 1F1 and 2F1 of a real argument are currently supported as power series. These make available Bessel J, Y, I and K functions of a real argument, as well as Airy Ai and Bi functions. Elliptic integrals and functions of a real argument are coming soon.

Basic support for integration of functions and ordinary differential equations are available for real values, as well as eigensystems of real matrices.

### Usage ###

Download the complete <a href="https://raw.githubusercontent.com/paulmasson/math/master/build/math.js">library</a> and include it in the body of a page:

```
<script src="math.js"></script>
```

A free online content delivery network can also be used:

```
<script src="https://cdn.rawgit.com/paulmasson/math/master/build/math.js"></script>
```

