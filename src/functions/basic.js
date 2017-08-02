
function factorial( n ) {

  if ( Number.isInteger(n) && n >= 0 ) {

    var result = 1;
    for ( var i = 2 ; i <= n ; i++ ) result *= i;
    return result;

  } else return gamma( n+1 );

}

function binomial( n, m ) {

  return factorial(n) / factorial(n-m) / factorial(m);

}


// Rounding functions

function roundTo( x, n ) {

  var exponent = Math.floor( Math.log10( Math.abs(x) ) );
  n = n - exponent;
  return Math.round( 10**n * x ) / 10**n;

}

function ceilTo( x, n ) {

  var exponent = Math.floor( Math.log10( Math.abs(x) ) );
  n = n - exponent;
  return Math.ceil( 10**n * x ) / 10**n;

}

function floorTo( x, n ) {

  var exponent = Math.floor( Math.log10( Math.abs(x) ) );
  n = n - exponent;
  return Math.floor( 10**n * x ) / 10**n;

}

function chop( x, tolerance=1e-10 ) {

  if ( isComplex(x) ) return { re: chop(x.re), im: chop(x.im) };

  if ( Math.abs(x) < tolerance ) x = 0;
  return x;

}


