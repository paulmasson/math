
// rounding functions

function roundTo( x, n ) {

  if ( x === 0 ) return x;

  if ( Array.isArray(x) ) {
    var v = vector( x.length );
    for ( var i = 0 ; i < x.length ; i++ ) v[i] = roundTo( x[i], n );
    return v;
  }

  var exponent = Math.floor( Math.log10( Math.abs(x) ) );
  n = n - exponent - 1;
  return Math.round( 10**n * x ) / 10**n;

}

function ceilTo( x, n ) {

  if ( x === 0 ) return x;

  if ( Array.isArray(x) ) {
    var v = vector( x.length );
    for ( var i = 0 ; i < x.length ; i++ ) v[i] = ceilTo( x[i], n );
    return v;
  }

  var exponent = Math.floor( Math.log10( Math.abs(x) ) );
  n = n - exponent - 1;
  return Math.ceil( 10**n * x ) / 10**n;

}

function floorTo( x, n ) {

  if ( x === 0 ) return x;

  if ( Array.isArray(x) ) {
    var v = vector( x.length );
    for ( var i = 0 ; i < x.length ; i++ ) v[i] = floorTo( x[i], n );
    return v;
  }

  var exponent = Math.floor( Math.log10( Math.abs(x) ) );
  n = n - exponent - 1;
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


function kronecker( i, j ) {

  return i === j ? 1 : 0;

}

