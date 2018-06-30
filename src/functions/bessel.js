
function besselJ( n, x ) {

  if ( isComplex(n) || isComplex(x) ) {

    if ( !isComplex(n) ) n = complex(n,0);
    if ( !isComplex(x) ) x = complex(x,0);

    if ( Number.isInteger(n.re) && n.re < 0 && n.im === 0 )
      return mul( pow(-1,n), besselJ( mul(-1,n), x ) );

    var product = div( pow( div(x,2), n ), gamma( add(n,1) ) );
    return mul( product, hypergeometric0F1( add(n,1), mul(-.25, pow(x,2) ) ) );

  } 

  if ( Number.isInteger(n) && n < 0 ) return (-1)**n * besselJ( -n, x );

  if ( !Number.isInteger(n) && x < 0 ) return besselJ( n, complex(x) );

  return (x/2)**n * hypergeometric0F1( n+1, -.25*x**2 ) / gamma(n+1);

}

function besselY( n, x ) {

  // for averaging over integer orders until write code for limit
  var delta = 1e-5;

  if ( isComplex(n) || isComplex(x) ) {

    if ( !isComplex(n) ) n = complex(n,0);
    if ( !isComplex(x) ) x = complex(x,0);

    if ( Number.isInteger(n.re) && n.im === 0 )
      return div( add( besselY( n.re + delta, x ), besselY( n.re - delta, x ) ), 2 );

    var sum = sub( mul( besselJ(n,x), cos( mul(n,pi) ) ), besselJ( mul(-1,n), x ) );
    return div( sum, sin( mul(n,pi) ) );

  }

  if ( x < 0 ) return besselY( n, complex(x) );

  if ( Number.isInteger(n) )
    return ( besselY( n + delta, x ) + besselY( n - delta, x ) ) / 2;

  return ( besselJ(n,x) * cos(n*pi) - besselJ(-n,x) ) / sin(n*pi);

}

function besselI( n, x ) {

  if ( isComplex(n) || isComplex(x) ) {

    if ( !isComplex(n) ) n = complex(n,0);
    if ( !isComplex(x) ) x = complex(x,0);

    if ( Number.isInteger(n.re) && n.re < 0 && n.im === 0 )
      return besselI( mul(-1,n), x );

    var product = div( pow( div(x,2), n ), gamma( add(n,1) ) );
    return mul( product, hypergeometric0F1( add(n,1), mul(.25, pow(x,2) ) ) );

  }

  if ( Number.isInteger(n) && n < 0 ) return besselI( -n, x );

  if ( !Number.isInteger(n) && x < 0 ) return besselI( n, complex(x) );

  return (x/2)**n * hypergeometric0F1( n+1, .25*x**2 ) / gamma(n+1);

}

function besselK( n, x ) {

  var useAsymptotic = 5;

  // for averaging over integer orders until write code for limit
  var delta = 1e-5;

  if ( isComplex(n) || isComplex(x) ) {

    if ( !isComplex(n) ) n = complex(n,0);
    if ( !isComplex(x) ) x = complex(x,0);

    // asymptotic form as per Johansson
    if ( abs(x) > useAsymptotic ) {

      var t1 = mul( sqrt( div( pi/2, x ) ), exp( mul(-1,x) ) );
      var t2 = hypergeometric2F0( add(n,.5), sub(.5,n), div( -1, mul(2,x) ) );

      return mul( t1, t2 );

    }

    if ( Number.isInteger(n.re) && n.im === 0 )
      return div( add( besselK( n.re + delta, x ), besselK( n.re - delta, x ) ), 2 );

    var product = div( pi/2, sin( mul(n,pi) ) );
    return mul( product, sub( besselI( mul(-1,n), x ), besselI(n,x) ) );

  }

  if ( x > useAsymptotic )
    return sqrt(pi/2/x) * exp(-x) * hypergeometric2F0( n+.5, .5-n, -1/2/x );

  if ( x < 0 ) return besselK( n, complex(x) );

  if ( Number.isInteger(n) )
    return ( besselK( n + delta, x ) + besselK( n - delta, x ) ) / 2;

  return pi/2 * ( besselI(-n,x) - besselI(n,x) ) / sin(n*pi);

}

function hankel1( n, x ) {

  return add( besselJ(n,x), mul( complex(0,1), besselY(n,x) ) );

}

function hankel2( n, x ) {

  return sub( besselJ(n,x), mul( complex(0,1), besselY(n,x) ) );

}


function airyAi( x ) {

  if ( isComplex(x) ) {

    return mul( 1/pi, mul( sqrt( div( x, 3 ) ), besselK( 1/3, mul( 2/3, pow( x, 3/2 ) ) ) ) );

  }

  if ( x === 0 ) return 1 / 3**(2/3) / gamma(2/3);

  if ( x < 0 ) return sqrt(-x) / 2 * ( besselJ( 1/3, 2/3*(-x)**(3/2) )
                                       - besselY( 1/3, 2/3*(-x)**(3/2) ) / sqrt(3) );

  return 1/pi * sqrt(x/3) * besselK( 1/3, 2/3*x**(3/2) );

}

function airyBi( x ) {

  if ( isComplex(x) ) {

    return mul( sqrt( div( x, 3 ) ), add( besselI( 1/3, mul( 2/3, pow( x, 3/2 ) ) ),
                                          besselI( -1/3, mul( 2/3, pow( x, 3/2 ) ) ) ) );

  }

  if ( x === 0 ) return 1 / 3**(1/6) / gamma(2/3);

  if ( x < 0 ) return -sqrt(-x) / 2 * ( besselJ( 1/3, 2/3*(-x)**(3/2) ) / sqrt(3)
                                        + besselY( 1/3, 2/3*(-x)**(3/2) ) );

  return sqrt(x/3) * ( besselI( 1/3, 2/3*x**(3/2) ) + besselI( -1/3, 2/3*x**(3/2) ) );

}
