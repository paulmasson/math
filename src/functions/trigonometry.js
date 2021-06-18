
// complex circular functions

function sin( x ) {

  if ( isArbitrary(x) ) {

    if ( isComplex(x) )

      return { re: mul( sin(x.re), cosh(x.im) ),
               im: mul( cos(x.re), sinh(x.im) ) };

    x = x % twoPi;

    // reduce to [-pi/2,pi/2] with successive reductions
    if ( x > halfPi ) return sin( onePi - x );
    if ( x < -halfPi ) return sin( -onePi - x );

    var s = x;
    var p = x;
    var i = arb2;

    while ( p !== 0n ) {
      p = div( mul( p, -arb1, x, x ), mul( i, i + arb1 ) );
      s += p;
      i += arb2;
    }

    return s;

  }

  if ( isComplex(x) )

    return { re: Math.sin(x.re) * Math.cosh(x.im),
             im: Math.cos(x.re) * Math.sinh(x.im) };

  return Math.sin(x);

}

function cos( x ) {

  if ( isArbitrary(x) ) {

    if ( isComplex(x) )

      return { re: mul( cos(x.re), cosh(x.im) ),
               im: mul( arbitrary(-1), sin(x.re), sinh(x.im) ) };

    x = x % twoPi;

    // reduce to [-pi/2,pi/2] with successive reductions
    if ( x > halfPi ) return -cos( onePi - x );
    if ( x < -halfPi ) return -cos( -onePi - x );

    var s = arb1;
    var p = arb1;
    var i = arb1;

    while ( p !== 0n ) {
      p = div( mul( p, -arb1, x, x ), mul( i, i + arb1 ) );
      s += p;
      i += arb2;
    }

    return s;

  }

  if ( isComplex(x) )

    return { re: Math.cos(x.re) * Math.cosh(x.im),
             im: -Math.sin(x.re) * Math.sinh(x.im) };

  return Math.cos(x);

}

function tan( x ) {

  if ( isComplex(x) ) return div( sin(x), cos(x) );

  return Math.tan(x);

 }

function cot( x ) {

  if ( isComplex(x) ) return div( cos(x), sin(x) );

  return 1 / Math.tan(x);

}

function sec( x ) {

  if ( isComplex(x) ) return div( 1, cos(x) );

  return 1 / Math.cos(x);

}

function csc( x ) {

  if ( isComplex(x) ) return div( 1, sin(x) );

  return 1 / Math.sin(x);

}


// inverse circular functions

function arcsin( x ) {

  if ( isComplex(x) ) {

    var s = sqrt( sub( 1, mul( x, x ) ) );
    s = add( mul( complex(0,1), x ), s ); 
    return mul( complex(0,-1), log( s ) );

  }

  if ( Math.abs(x) <= 1 ) return Math.asin(x);

  return arcsin( complex(x) );

}

function arccos( x ) {

  if ( isComplex(x) ) {

    return sub( pi/2, arcsin(x) );

  }

  if ( Math.abs(x) <= 1 ) return Math.acos(x);

  return arccos( complex(x) );

}

function arctan( x ) {

  if ( isComplex(x) ) {

    var s = sub( log( sub( 1, mul( complex(0,1), x ) ) ),
                 log( add( 1, mul( complex(0,1), x ) ) ) );
    return mul( complex(0,.5), s );

  }

  return Math.atan(x);

}

function arccot( x ) {

  if ( isComplex(x) ) return arctan( div( 1, x ) );

  return Math.atan( 1/x );

}

function arcsec( x ) {

  if ( isComplex(x) ) return arccos( div( 1, x ) );

  if ( Math.abs(x) >= 1 ) return Math.acos( 1/x );

  return arcsec( complex(x) );

}

function arccsc( x ) {

  if ( isComplex(x) ) return arcsin( div( 1, x ) );

  if ( Math.abs(x) >= 1 ) return Math.asin( 1/x );

  return arccsc( complex(x) );

}


// complex hyperbolic functions

function sinh( x ) {

  if ( isArbitrary(x) ) return div( sub( exp(x), exp( mul(-arb1,x) ) ), arb2 );

  if ( isComplex(x) )

    return { re: Math.sinh(x.re) * Math.cos(x.im),
             im: Math.cosh(x.re) * Math.sin(x.im) };

  return Math.sinh(x);

}

function cosh( x ) {

  if ( isArbitrary(x) ) return div( add( exp(x), exp( mul(-arb1,x) ) ), arb2 );

  if ( isComplex(x) )

    return { re: Math.cosh(x.re) * Math.cos(x.im),
             im: Math.sinh(x.re) * Math.sin(x.im) };

  return Math.cosh(x);

}

function tanh( x ) {

  if ( isComplex(x) ) return div( sinh(x), cosh(x) );

  return Math.tanh(x);

}

function coth( x ) {

  if ( isComplex(x) ) return div( cosh(x), sinh(x) );

  return 1 / Math.tanh(x);

}

function sech( x ) {

  if ( isComplex(x) ) return div( 1, cosh(x) );

  return 1 / Math.cosh(x);

}

function csch( x ) {

  if ( isComplex(x) ) return div( 1, sinh(x) );

  return 1 / Math.sinh(x);

}


// inverse hyperbolic functions

function arcsinh( x ) {

  if ( isComplex(x) ) {

    var s = sqrt( add( mul( x, x ), 1 ) );
    s = add( x, s );
    return log( s );

  }

  return Math.asinh(x);

}

function arccosh( x ) {

  if ( isComplex(x) ) {

    var s = mul( sqrt( add( x, 1 ) ), sqrt( sub( x, 1 ) ) );
    s = add( x, s ); 
    return log( s );

  }

  if ( x >= 1 ) return Math.acosh(x);

  return arccosh( complex(x) );

}

function arctanh( x ) {

  if ( isComplex(x) ) {

    var s = sub( log( add( 1, x ) ), log( sub( 1, x ) ) );
    return mul( .5, s );

  }

  if ( Math.abs(x) <= 1 ) return Math.atanh(x);

  return arctanh( complex(x) );

}

function arccoth( x ) {

  if ( isComplex(x) ) {

    if ( x.re === 0 && x.im === 0 ) throw Error( 'Indeterminate arccoth value' );

    return arctanh( div( 1, x ) );

  }

  if ( Math.abs(x) > 1 ) return Math.atanh( 1/x );

  return arccoth( complex(x) );

}

function arcsech( x ) {

  if ( isComplex(x) ) {

    if ( x.re === 0 && x.im === 0 ) throw Error( 'Indeterminate arcsech value' );

    // adjust for branch cut along negative axis
    if ( x.im === 0 ) x.im = -Number.MIN_VALUE;

    return arccosh( div( 1, x ) );

  }

  if ( x > 0 && x < 1 ) return Math.acosh( 1/x );

  return arcsech( complex(x) );

}

function arccsch( x ) {

  if ( isComplex(x) ) {

    return arcsinh( div( 1, x ) );

  }

  return Math.asinh( 1/x );

}


// miscellaneous

function sinc( x ) {

  if ( isComplex(x) ) {

    if ( x.re === 0 && x.im === 0 ) return complex(1);

    return div( sin(x), x );

  }

  if ( x === 0 ) return 1;

  return Math.sin(x) / x;

}

function haversine( x ) { return div( sub( 1, cos(x) ), 2 ); }

function inverseHaversine( x ) { return arccos( sub( 1, mul(2,x) ) ); }


// analyticphysics.com / The Complex Gudermannian Function

function gudermannian( x ) { return mul( 2, arctan( tanh( div(x,2) ) ) ); }

function inverseGudermannian( x ) { return mul( 2, arctanh( tan( div(x,2) ) ) ); }

