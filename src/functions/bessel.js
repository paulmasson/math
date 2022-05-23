
function besselJ( n, x ) {

  if ( isComplex(n) || isComplex(x) ) {

    if ( isNegativeInteger(n) ) return mul( pow(-1,n), besselJ( mul(-1,n), x ) );

    var product = div( pow( div(x,2), n ), gamma( add(n,1) ) );
    return mul( product, hypergeometric0F1( add(n,1), mul(-.25, pow(x,2) ) ) );

  } 

  if ( isNegativeInteger(n) ) return (-1)**n * besselJ( -n, x );

  if ( !Number.isInteger(n) && x < 0 ) return besselJ( n, complex(x) );

  return (x/2)**n * hypergeometric0F1( n+1, -.25*x**2 ) / gamma(n+1);

}

function besselJZero( n, m, derivative=false ) {

  if ( isComplex(n) || isComplex(m) )
    throw Error( 'Complex arguments not supported in Bessel zero' );

  if ( n < 0 ) throw Error( 'Negative order for Bessel zero' );
  if ( !isPositiveInteger(m) ) throw Error( 'Unsupported index for Bessel zero' );

  // approximations from dlmf.nist.gov/10.21#vi
  var delta = pi/4;

  if ( derivative ) {

    if ( n === 0 && m === 1 ) return 0;

    var b = ( m + n/2 - 3/4 ) * pi;
    var e = b - ( 4*n**2 + 3 ) / ( 8*b );

    // keep search evaluation real
    return findRoot( x => diff( x => besselJ(n,x), x ), [ e-delta < 0 ? 0 : e-delta, e+delta ] );

  } else {

    var a = ( m + n/2 - 1/4 ) * pi;
    var e = a - ( 4*n**2 - 1 ) / ( 8*a );

    return findRoot( x => besselJ(n,x), [ e-delta, e+delta ] );

  }

}

function besselY( n, x ) {

  if ( isComplex(n) || isComplex(x) ) {

    // dlmf.nist.gov/10.2.3
    if ( isInteger(n) )
      return div( add( diff( n => besselJ(n,x), n ),
                       mul( pow(-1,n), diff( n => besselJ(n,x), neg(n) ) ) ), pi );

    var sum = sub( mul( besselJ(n,x), cos( mul(n,pi) ) ), besselJ( mul(-1,n), x ) );
    return div( sum, sin( mul(n,pi) ) );

  }

  if ( x < 0 ) return besselY( n, complex(x) );

  // dlmf.nist.gov/10.2.3
  if ( Number.isInteger(n) )
    return ( diff( n => besselJ(n,x), n ) + (-1)**n * diff( n => besselJ(n,x), -n ) ) / pi;

  return ( besselJ(n,x) * cos(n*pi) - besselJ(-n,x) ) / sin(n*pi);

}

function besselYZero( n, m, derivative=false ) {

  if ( isComplex(n) || isComplex(m) )
    throw Error( 'Complex arguments not supported in Bessel zero' );

  if ( n < 0 ) throw Error( 'Negative order for Bessel zero' );
  if ( !isPositiveInteger(m) ) throw Error( 'Unsupported index for Bessel zero' );

  // approximations from dlmf.nist.gov/10.21#vi
  var delta = pi/4;

  if ( derivative ) {

    var b = ( m + n/2 - 1/4 ) * pi;
    var e = b - ( 4*n**2 + 3 ) / ( 8*b );

    return findRoot( x => diff( x => besselY(n,x), x ), [ e-delta, e+delta ] );

  } else {

    var a = ( m + n/2 - 3/4 ) * pi;
    var e = a - ( 4*n**2 - 1 ) / ( 8*a );

    return findRoot( x => besselY(n,x), [ e-delta, e+delta ] );

  }

}

function besselI( n, x ) {

  if ( isComplex(n) || isComplex(x) ) {

    if ( isNegativeInteger(n) ) return besselI( mul(-1,n), x );

    var product = div( pow( div(x,2), n ), gamma( add(n,1) ) );
    return mul( product, hypergeometric0F1( add(n,1), mul(.25, pow(x,2) ) ) );

  }

  if ( isNegativeInteger(n) ) return besselI( -n, x );

  if ( !Number.isInteger(n) && x < 0 ) return besselI( n, complex(x) );

  return (x/2)**n * hypergeometric0F1( n+1, .25*x**2 ) / gamma(n+1);

}

function besselK( n, x ) {

  var useAsymptotic = 10;

  if ( isComplex(n) || isComplex(x) ) {

    // asymptotic form as per Johansson arxiv.org/abs/1606.06977
    if ( abs(x) > useAsymptotic ) {

      var t1 = mul( sqrt( div( pi/2, x ) ), exp( neg(x) ) );
      var t2 = hypergeometric2F0( add(n,.5), sub(.5,n), div(-.5,x) );

      return mul( t1, t2 );

    }

    // dlmf.nist.gov/10.27.5
    if ( isInteger(n) )
      return mul( pow(-1,add(n,1)), .5,
                  add( diff( n => besselI(n,x), n ), diff( n => besselI(n,x), neg(n) ) ) );

    var product = div( pi/2, sin( mul(n,pi) ) );
    return mul( product, sub( besselI( mul(-1,n), x ), besselI(n,x) ) );

  }

  if ( x > useAsymptotic )
    return sqrt(pi/2/x) * exp(-x) * hypergeometric2F0( n+.5, .5-n, -.5/x );

  if ( x < 0 ) return besselK( n, complex(x) );

  // dlmf.nist.gov/10.27.5
  if ( Number.isInteger(n) )
    return (-1)**(n+1)/2 * ( diff( n => besselI(n,x), n ) + diff( n => besselI(n,x), -n ) );

  return pi/2 * ( besselI(-n,x) - besselI(n,x) ) / sin(n*pi);

}

function hankel1( n, x ) {

  return add( besselJ(n,x), mul( complex(0,1), besselY(n,x) ) );

}

function hankel2( n, x ) {

  return sub( besselJ(n,x), mul( complex(0,1), besselY(n,x) ) );

}


// dlmf.nist.gov/9.2.ii and dlmf.nist.gov/9.6.i

function airyAi( x ) {

  if ( isComplex(x) ) {

    if ( isZero(x) ) return complex( 1 / 3**(2/3) / gamma(2/3) );

    if ( x.re < 0 ) {

      var z = mul( 2/3, pow( neg(x), 3/2 ) );
      return mul( 1/3, sqrt(neg(x)), add( besselJ( 1/3, z ), besselJ( -1/3, z ) ) );

    }

    var z = mul( 2/3, pow( x, 3/2 ) );
    return mul( 1/pi, sqrt( div( x, 3 ) ), besselK( 1/3, z ) );

  }

  if ( x === 0 ) return 1 / 3**(2/3) / gamma(2/3);

  if ( x < 0 ) {

    var z = 2/3 * (-x)**(3/2);
    return sqrt(-x) / 3 * ( besselJ( 1/3, z ) + besselJ( -1/3, z ) );

  }

  var z = 2/3 * x**(3/2);
  return 1/pi * sqrt(x/3) * besselK( 1/3, z );

}

function airyAiPrime( x ) {

  if ( isComplex(x) ) {

    if ( isZero(x) ) return complex( -1 / 3**(1/3) / gamma(1/3) );

    if ( x.re < 0 ) {

      var z = mul( 2/3, pow( neg(x), 3/2 ) );
      return mul( 1/3, x, sub( besselJ( -2/3, z ), besselJ( 2/3, z ) ) );

    }

    var z = mul( 2/3, pow( x, 3/2 ) );
    return mul( -1/pi/sqrt(3), x, besselK( 2/3, z ) );

  }

  if ( x === 0 ) return -1 / 3**(1/3) / gamma(1/3);

  if ( x < 0 ) {

    var z = 2/3 * (-x)**(3/2);
    return x/3 * ( besselJ( -2/3, z ) - besselJ( 2/3, z ) );

  }

  var z = 2/3 * x**(3/2);
  return -1/pi/sqrt(3) * x * besselK( 2/3, z );

}

function airyBi( x ) {

  if ( isComplex(x) ) {

    if ( isZero(x) ) return complex( 1 / 3**(1/6) / gamma(2/3) );

    if ( x.re < 0 ) {

      var z = mul( 2/3, pow( neg(x), 3/2 ) );
      return mul( sqrt( div(neg(x),3) ), sub( besselJ( -1/3, z ), besselJ( 1/3, z ) ) );

    }

    var z = mul( 2/3, pow( x, 3/2 ) );
    return mul( sqrt( div( x, 3 ) ), add( besselI( 1/3, z ), besselI( -1/3, z ) ) );

  }

  if ( x === 0 ) return 1 / 3**(1/6) / gamma(2/3);

  if ( x < 0 ) {

    var z = 2/3 * (-x)**(3/2);
    return sqrt(-x/3) * ( besselJ( -1/3, z ) - besselJ( 1/3, z ) );

  }

  var z = 2/3 * x**(3/2);
  return sqrt(x/3) * ( besselI( 1/3, z ) + besselI( -1/3, z ) );

}

function airyBiPrime( x ) {

  if ( isComplex(x) ) {

    if ( isZero(x) ) return complex( 3**(1/6) / gamma(1/3) );

    if ( x.re < 0 ) {

      var z = mul( 2/3, pow( neg(x), 3/2 ) );
      return mul( 1/sqrt(3), neg(x), add( besselJ( 2/3, z ), besselJ( -2/3, z ) ) );

    }

    var z = mul( 2/3, pow( x, 3/2 ) );
    return mul( 1/sqrt(3), x, add( besselI( 2/3, z ), besselI( -2/3, z ) ) );

  }

  if ( x === 0 ) return 3**(1/6) / gamma(1/3);

  if ( x < 0 ) {

    var z = 2/3 * (-x)**(3/2);
    return -x/sqrt(3) * ( besselJ( 2/3, z ) + besselJ( -2/3, z ) );

  }

  var z = 2/3 * x**(3/2);
  return x/sqrt(3) * ( besselI( 2/3, z ) + besselI( -2/3, z ) );

}


function sphericalBesselJ( n, x ) {

  return mul( div( sqrt(pi/2), sqrt(x) ), besselJ( add( n, .5 ), x ) );

}

function sphericalBesselY( n, x ) {

  return mul( div( sqrt(pi/2), sqrt(x) ), besselY( add( n, .5 ), x ) );

}

function sphericalHankel1( n, x ) {

  return add( sphericalBesselJ(n,x), mul( complex(0,1), sphericalBesselY(n,x) ) );

}

function sphericalHankel2( n, x ) {

  return sub( sphericalBesselJ(n,x), mul( complex(0,1), sphericalBesselY(n,x) ) );

}


function struveH( n, x ) {

  return mul( pow( x, add(n,1) ), inv( mul( pow(2,n), sqrt(pi), gamma( add(n,3/2) ) ) ),
              hypergeometric1F2( 1, 3/2, add(n,3/2), mul( -1/4, pow(x,2) ) ) );

}

function struveL( n, x ) {

  // one sign different from struveH

  return mul( pow( x, add(n,1) ), inv( mul( pow(2,n), sqrt(pi), gamma( add(n,3/2) ) ) ),
              hypergeometric1F2( 1, 3/2, add(n,3/2), mul( 1/4, pow(x,2) ) ) );

}

