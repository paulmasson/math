
function chop( x, tolerance=1e-10 ) {

  if ( Array.isArray(x) ) {
    var v = vector( x.length );
    for ( var i = 0 ; i < x.length ; i++ ) v[i] = chop( x[i] );
    return v;
  }

  if ( isComplex(x) ) return complex( chop(x.re), chop(x.im) );

  if ( Math.abs(x) < tolerance ) x = 0;

  return x;

}


function kronecker( i, j ) {

  return i === j ? 1 : 0;

}


function piecewise() {

  var pieces = arguments;

  return function( x ) {

    for ( var i = 0 ; i < pieces.length ; i++ ) {
      var domain = pieces[i][1];
      if ( x >= domain[0] && x <= domain[1] ) return pieces[i][0](x);
    }

    return 0;

  }

}

