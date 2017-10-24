
// complex circular functions

function sin( x ) {

  if ( isComplex(x) )

    return { re: Math.sin(x.re) * Math.cosh(x.im),
             im: Math.cos(x.re) * Math.sinh(x.im) };

  else return Math.sin(x);

}

function cos( x ) {

  if ( isComplex(x) )

    return { re: Math.cos(x.re) * Math.cosh(x.im),
             im: -Math.sin(x.re) * Math.sinh(x.im) };

  else return Math.cos(x);

}

function tan( x ) {

  if ( isComplex(x) ) return div( sin(x), cos(x) );

  else return Math.tan(x);

 }

function cot( x ) {

 if ( isComplex(x) ) return div( cos(x), sin(x) );

  else return 1 / Math.tan(x);

}

function sec( x ) {

  if ( isComplex(x) ) return div( 1, cos(x) );

  else return 1 / Math.cos(x);

}

function csc( x ) {

  if ( isComplex(x) ) return div( 1, sin(x) );

  else return 1 / Math.sin(x);

}


// inverse circular functions

function arcsin( x ) {

  if ( isComplex(x) ) {

    var s = sqrt( sub( 1, mul( x, x ) ) );
    s = add( mul( complex(0,1), x ), s ); 
    return mul( complex(0,-1), log( s ) );

  } else if ( Math.abs(x) <= 1 ) return Math.asin(x);

  else return arcsin( complex(x) );

}

function arccos( x ) {

  if ( isComplex(x) ) {

    return sub( pi/2, arcsin(x) );

  } else if ( Math.abs(x) <= 1 ) return Math.acos(x);

  else return arccos( complex(x) );

}

function arctan( x ) {

  if ( isComplex(x) ) {

    var s = sub( log( sub( 1, mul( complex(0,1), x ) ) ),
                 log( add( 1, mul( complex(0,1), x ) ) ) );
    return mul( complex(0,1/2), s );

  } else return Math.atan(x);

}

function arccot( x ) {

  if ( isComplex(x) ) return arctan( div( 1, x ) );

  else return Math.atan( 1/x );

}

function arcsec( x ) {

  if ( isComplex(x) ) return arccos( div( 1, x ) );

  else if ( Math.abs(x) >= 1 ) return Math.acos( 1/x );

  else return arcsec( complex(x) );

}

function arccsc( x ) {

  if ( isComplex(x) ) return arcsin( div( 1, x ) );

  else if ( Math.abs(x) >= 1 ) return Math.asin( 1/x );

  else return arccsc( complex(x) );

}


// complex hyperbolic functions

function sinh( x ) {

  if ( isComplex(x) )

    return { re: Math.sinh(x.re) * Math.cos(x.im),
             im: Math.cosh(x.re) * Math.sin(x.im) };

  else return Math.sinh(x);

}

function cosh( x ) {

  if ( isComplex(x) )

    return { re: Math.cosh(x.re) * Math.cos(x.im),
             im: Math.sinh(x.re) * Math.sin(x.im) };

  else return Math.cosh(x);

}

function tanh( x ) {

  if ( isComplex(x) ) return div( sinh(x), cosh(x) );

  else return Math.tanh(x);

}

function coth( x ) {

  if ( isComplex(x) ) return div( cosh(x), sinh(x) );

  else return 1 / Math.tanh(x);

}

function sech( x ) {

  if ( isComplex(x) ) return div( 1, cosh(x) );

  else return 1 / Math.cosh(x);

}

function csch( x ) {

  if ( isComplex(x) ) return div( 1, sinh(x) );

  else return 1 / Math.sinh(x);

}


// inverse hyperbolic functions

function arcsinh( x ) {

  if ( isComplex(x) ) {

    var s = sqrt( add( mul( x, x ), 1 ) );
    s = add( x, s );
    return log( s );

  } else return Math.asinh(x);

}

function arccosh( x ) {

  if ( isComplex(x) ) {

    var s = mul( sqrt( add( x, 1 ) ), sqrt( sub( x, 1 ) ) );
    s = add( x, s ); 
    return log( s );

  } else if ( x >= 1 ) return Math.acosh(x);

  else return arccosh( complex(x) );

}

function arctanh( x ) {

  if ( isComplex(x) ) {

    var s = sub( log( add( 1, x ) ), log( sub( 1, x ) ) );
    return mul( 1/2, s );

  } else if ( Math.abs(x) <= 1 ) return Math.atanh(x);

  else return arctanh( complex(x) );

}

function arccoth( x ) {

  if ( isComplex(x) ) {

    if ( x.re === 0 && x.im === 0 ) throw 'Indeterminate value';

    return arctanh( div( 1, x ) );

  } else if ( Math.abs(x) > 1 ) return Math.atanh( 1/x );

  else return arccoth( complex(x) );

}

function arcsech( x ) {

  if ( isComplex(x) ) {

    if ( x.re === 0 && x.im === 0 ) throw 'Indeterminate value';

    // adjust for branch cut along negative axis
    if ( x.im === 0 ) x.im = -Number.MIN_VALUE;

    return arccosh( div( 1, x ) );

  } else if ( x > 0 && x < 1 ) return Math.acosh( 1/x );

  else return arcsech( complex(x) );

}

function arccsch( x ) {

  if ( isComplex(x) ) {

    return arcsinh( div( 1, x ) );

  } else return Math.asinh( 1/x );

}

