
function zeta( x ) {

  // functional equation
  if ( x < 0 ) return 2**x * pi**(x-1) * sin(pi*x/2) * gamma(1-x) * zeta(1-x); 

  // summation of Dirichlet eta function
  var s = 0;
  for ( var i = 1 ; i < 1e7 ; i++ ) {
    var last = s;
    s -= (-1)**i / i**x;
    if ( Math.abs(s-last) < 1e-8 ) break;
  }

  return s / ( 1 - 2**(1-x) );

}


