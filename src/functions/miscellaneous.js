
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

function round( x, y ) {

  if ( arguments.length === 2 ) return mul( y, round( div(x,y) ) );

  if ( isComplex(x) ) return complex( Math.round(x.re), Math.round(x.im) );

  return Math.round(x);

}

function ceiling( x ) {

  if ( isComplex(x) ) return complex( Math.ceil(x.re), Math.ceil(x.im) );

  return Math.ceil(x);

}

function floor( x ) {

  if ( isComplex(x) ) return complex( Math.floor(x.re), Math.floor(x.im) );

  return Math.floor(x);

}

function sign( x ) {

  if ( isZero(x) ) return x;

  return div( x, abs(x) );

}

function integerPart( x ) {

  if ( isComplex(x) ) return complex( Math.trunc(x.re), Math.trunc(x.im) );

  return Math.trunc(x);

}

function fractionalPart( x ) { return sub( x, integerPart(x) ); }


function random( x, y ) {

  if ( arguments.length === 0 ) return Math.random();

  if ( arguments.length === 1 )
    if ( isComplex(x) )
      return { re: x.re * Math.random(), im: x.im * Math.random() };
    else return x * Math.random();

  if ( isComplex(x) || isComplex(y) ) {
    if ( !isComplex(x) ) x = complex(x);
    if ( !isComplex(y) ) y = complex(y);
    return { re: random( x.re, y.re ), im: random( x.im, y.im ) };
  }

  return Math.abs(x-y) * Math.random() + ( x > y ? y : x );

}


function kronecker( i, j ) {

  if ( arguments.length === 2 ) {

    if ( isComplex(i) || isComplex(j) ) {

      if ( !isComplex(i) ) i = complex(i);
      if ( !isComplex(j) ) j = complex(j);

      return kronecker( i.re, j.re) * kronecker( i.im, j.im );

    }

    return i === j ? 1 : 0;

  }

  var result = kronecker( i, j );

  for ( var k = 2 ; k < arguments.length ; k++ )
    result *= kronecker( i, arguments[k] );

  return result;

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

