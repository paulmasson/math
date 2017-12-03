
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

  return (x/2)**n * hypergeometric0F1( n+1, -.25*x**2 ) / gamma(n+1);

}

function besselY( n, x ) {

  // mpmath displaces integer orders to avoid poles
  // displacement value is dependent on transition to asymptotic
  //   hypergeometric solution, where lose a few digits of precision

  var displacement = 1e-7;

  if ( isComplex(n) || isComplex(x) ) {

    if ( !isComplex(n) ) n = complex(n,0);
    if ( !isComplex(x) ) x = complex(x,0);

    if ( Number.isInteger(n.re) ) n = add( n, displacement );

    var sum = sub( mul( besselJ(n,x), cos( mul(n,pi) ) ), besselJ( mul(-1,n), x ) );
    return div( sum, sin( mul(n,pi) ) );

  } else {

    if ( Number.isInteger(n) ) n += displacement;

    return ( besselJ(n,x) * cos(n*pi) - besselJ(-n,x) ) / sin(n*pi);

  }

}

function besselI( n, x ) {

  if ( isComplex(n) || isComplex(x) ) {

    if ( !isComplex(n) ) n = complex(n,0);
    if ( !isComplex(x) ) x = complex(x,0);

    var product = div( pow( div(x,2), n ), gamma( add(n,1) ) );
    return mul( product, hypergeometric0F1( add(n,1), mul(.25, pow(x,2) ) ) );

  } else return (x/2)**n * hypergeometric0F1( n+1, .25*x**2 ) / gamma(n+1);

}

function besselK( n, x ) {

  var useAsymptotic = 10;

  if ( isComplex(n) || isComplex(x) ) {

    if ( !isComplex(n) ) n = complex(n,0);
    if ( !isComplex(x) ) x = complex(x,0);

    // asymptotic form as per Johansson
    if ( abs(x) > useAsymptotic ) {

      var t1 = mul( sqrt( div( pi/2, x ) ), exp( mul(-1,x) ) );
      var t2 = hypergeometric2F0( add(n,.5), sub(.5,n), div( -1, mul(2,x) ) );

      return mul( t1, t2 );

    }

    var product = div( pi/2, sin( mul(n,pi) ) );
    return mul( product, sub( besselI( mul(-1,n), x ), besselI(n,x) ) );

  } else {

    if ( x > useAsymptotic )
      return sqrt(pi/2/x) * exp(-x) * hypergeometric2F0( n+.5, .5-n, -1/2/x );

    return pi/2 * ( besselI(-n,x) - besselI(n,x) ) / sin(n*pi);

  }

}


function airyAi( x ) {

  if ( x > 0 ) return 1/pi * sqrt(x/3) * besselK(1/3, 2/3*x**(3/2));

  if ( x === 0 ) return 1 / 3**(2/3) / gamma(2/3);

  if ( x < 0 ) return sqrt(-x) / 2 * (  besselJ(1/3, 2/3*(-x)**(3/2))
                                       - besselY(1/3, 2/3*(-x)**(3/2)) / sqrt(3) );

}

function airyBi( x ) {

  if ( x > 0 ) return sqrt(x/3) * ( besselI(1/3, 2/3*x**(3/2))
                                    + besselI(-1/3, 2/3*x**(3/2)) );

  if ( x === 0 ) return 1 / 3**(1/6) / gamma(2/3);

  if ( x < 0 ) return -sqrt(-x) / 2 * ( besselJ(1/3, 2/3*(-x)**(3/2)) / sqrt(3)
                                        + besselY(1/3, 2/3*(-x)**(3/2)) );

}
