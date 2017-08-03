
// rounding functions

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

  if ( Array.isArray(x) ) {
    var v = vector( x.length );
    for ( var i = 0 ; i < x.length ; i++ ) v[i] = chop( x[i] );
    return v;
  }

  if ( isComplex(x) ) return { re: chop(x.re), im: chop(x.im) };

  if ( Math.abs(x) < tolerance ) x = 0;
  return x;

}


