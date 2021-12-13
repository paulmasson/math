
function factorial( n ) {

  if ( isComplex(n) ) {

    if ( n.im === 0 && isPositiveIntegerOrZero(n.re) ) return complex( factorial(n.re) );

    return gamma( add(n,1) );

  }

  if ( isPositiveIntegerOrZero(n) ) {

    if ( factorialCache[n] ) return factorialCache[n];

    var last = factorialCache.length - 1;
    var result = factorialCache[last];

    for ( var i = last + 1 ; i <= n ; i++ ) {
      result *= i;
      factorialCache[i] = result;
    }

    return result;

  }

  return gamma( n+1 );

}

function factorial2( n ) {

  if ( isZero(n) ) return 1;

  if ( isPositiveInteger(n) ) {

    // bitwise test for odd integer, upward recursion for possible caching
    var result = n & 1 ? 1 : 2;
    for ( var i = result + 2 ; i <= n ; i += 2 ) result *= i;
    return result;

  }

  var f1 = pow( 2, div(n,2) );
  var f2 = pow( pi/2, div( sub( cos(mul(pi,n)), 1 ), 4 ) );
  var f3 = gamma( add( div(n,2) , 1 ) );

  return mul( f1, f2, f3 );

}

function binomial( n, m ) {

  if ( Number.isInteger(m) && m < 0 && n >= 0 ) return 0;

  if ( Number.isInteger(n) && Number.isInteger(m) && n >= 0 && m > n ) return 0;

  if ( isComplex(n) || isComplex(m) )
    return div( factorial(n), mul( factorial( sub(n,m) ), factorial(m) ) );

  return factorial(n) / factorial(n-m) / factorial(m);

}

function pochhammer( x, n ) {

  var one = isArbitrary(x) ? arb1 : 1;

  if ( isPositiveInteger(n) ) {

    var result = x, current = one, count = 1;

    while ( count < n ) {
      result = mul( result, add( x, current ) );
      current = add( current, one );
      count++;
    }

    return result;

  }

  if ( isZero(n) ) return one;

  return div( gamma( add(x,n) ), gamma(x) );

}


// log of gamma less likely to overflow than gamma
// Lanczos approximation as evaluated by Paul Godfrey

function logGamma( x ) {

  var c = [ 57.1562356658629235, -59.5979603554754912, 14.1360979747417471,
            -.491913816097620199, .339946499848118887e-4, .465236289270485756e-4,
            -.983744753048795646e-4, .158088703224912494e-3, -.210264441724104883e-3,
            .217439618115212643e-3, -.164318106536763890e-3, .844182239838527433e-4,
            -.261908384015814087e-4, .368991826595316234e-5 ];

  if ( isComplex(x) ) {

    if ( isNegativeIntegerOrZero(x) ) throw Error( 'Gamma function pole' );

    // reflection formula with modified Hare correction to imaginary part
    if ( x.re < 0 ) {
      var t = sub( log( div( pi, sin( mul(pi,x) ) ) ), logGamma( sub(1,x) ) );
      var s = x.im < 0 ? -1 : 1;
      var d = x.im === 0 ? 1/4 : 0;
      var k = Math.ceil( x.re/2 - 3/4 + d );
      return add( t, complex( 0, 2*s*k*pi ) );
    }

    var t = add( x, 5.24218750000000000 );
    t = sub( mul( add( x, .5 ), log(t)), t );
    var s = .999999999999997092;
    for ( var j = 0 ; j < 14 ; j++ ) s = add( s, div( c[j], add( x, j+1 ) ) );
    var u = add( t, log( mul( 2.5066282746310005, div( s, x ) ) ) );

    // adjustment to keep imaginary part on same sheet
    if ( s.re < 0 ) {
      if( x.im < 0 && div(s,x).im < 0 ) u = add( u, complex(0,2*pi) );
      if( x.im > 0 && div(s,x).im > 0 ) u = add( u, complex(0,-2*pi) );
    }

    return u;

  } else {

    if ( isNegativeIntegerOrZero(x) ) throw Error( 'Gamma function pole' ); 

    var t = x + 5.24218750000000000;
    t = ( x + .5 ) * log(t) - t;
    var s = .999999999999997092;
    for ( var j = 0 ; j < 14 ; j++ ) s += c[j] / (x+j+1);
    return t + log( 2.5066282746310005 * s / x );

  }

}

function gamma( x, y, z ) {

  if ( arguments.length === 2 ) {

    if ( isZero(x) ) {

      if ( isZero(y) ) throw Error( 'Gamma function pole' );

      // combination of logarithms adds/subtracts complex(0,pi)
      var sign = y.im > 0 ? -1 : y.im < 0 ? 1 : 0;

      var result = add( neg(expIntegralEi(neg(y))), complex(0,sign*pi) );

      if ( !isComplex(y) && y > 0 ) return result.re;

      return result;

    }

    // dlmf.nist.gov/8.4.15
    if ( isNegativeInteger(x) ) {

      var n = isComplex(x) ? -x.re : -x;
      var t = mul( exp(neg(y)), summation( k => div( (-1)**k*factorial(k), pow(y,k+1) ), [0,n-1] ) );

      // dlmf.nist.gov/8.4.4
      var result = mul( (-1)**n/factorial(n), sub( gamma(0,y), t ) );

      if ( isComplex(x) && !isComplex(result) ) return complex(result); // complex in, complex out

      return result;

    }

    return sub( gamma(x), gamma(x,0,y) );

  }

  if ( arguments.length === 3 ) {

    if ( !isZero(y) ) return sub( gamma(x,0,z), gamma(x,0,y) );

    return mul( pow(z,x), inv(x), hypergeometric1F1( x, add(x,1), neg(z) ) );

  }

  if ( isPositiveInteger(x) ) return factorial( sub(x,1) );

  // logGamma complex on negative axis
  if ( !isComplex(x) && x < 0 ) return exp( logGamma( complex(x) ) ).re;

  return exp( logGamma(x) );

}

function gammaRegularized( x, y, z ) {

  if ( arguments.length === 3 ) return div( gamma(x,y,z), gamma(x) );

  return div( gamma(x,y), gamma(x) );

}

function beta( x, y, z, w ) {

  if ( arguments.length === 4 )

    return sub( beta(y,z,w), beta(x,z,w) );

  if ( arguments.length === 3 )

    return mul( pow(x,y), inv(y), hypergeometric2F1( y, sub(1,z), add(y,1), x ) );

  return div( mul( gamma(x), gamma(y) ), gamma( add(x,y) ) ); 

}

function betaRegularized( x, y, z, w ) {

  if ( arguments.length === 4 )

    return div( beta(x,y,z,w), beta(z,w) );

  return div( beta(x,y,z), beta(y,z) );

}

function polygamma( n, x ) {

  if ( arguments.length === 1 ) return digamma(x);

  if ( !isPositiveInteger(n) ) throw Error( 'Unsupported polygamma index' );

  return mul( (-1)**(n+1) * factorial(n), hurwitzZeta( n+1, x ) );

}


function digamma( x ) {

  return diff( x => logGamma(x), x );

}


function erf( x ) {

  var useAsymptotic = 5;

  var absArg = Math.abs( arg(x) );

  if ( abs(x) > useAsymptotic && ( absArg < pi/4 || absArg > 3*pi/4 ) )

    return sub( 1, erfc(x) );

  return mul( 2/sqrt(pi), x, hypergeometric1F1( .5, 1.5, neg(mul(x,x)) ) );

}

function erfc( x ) {

  var useAsymptotic = 5;

  var absArg = Math.abs( arg(x) );

  if ( abs(x) > useAsymptotic && ( absArg < pi/4 || absArg > 3*pi/4 ) ) {

    // as per dlmf.nist.gov/7.12.1 this could be an independent sum for minor improvement
    // these numbers are tiny and need to stay in this function even though
    //   there is some code duplication with erf

    var t = mul( 1/sqrt(pi), exp( neg(mul(x,x)) ), inv(x),
                 hypergeometric2F0( .5, 1, neg(inv(mul(x,x))) ) );

    if ( x.re < 0 || x < 0 ) return add( 2, t );

    return t;

  }

  return sub( 1, erf(x) );

}

function erfi( x ) {

  return mul( complex(0,-1), erf( mul( complex(0,1), x ) ) );

}


function fresnelS( x ) {

  // can also be evaluated with hypergeometric1F2

  var m1 = hypergeometric1F1( .5, 1.5, mul( complex(0,pi/2), pow(x,2) ) );
  var m2 = hypergeometric1F1( .5, 1.5, mul( complex(0,-pi/2), pow(x,2) ) );

  var result = mul( x, sub( m1, m2 ), complex(0,-.5) );

  if ( isComplex(x) ) return result;

  return result.re;

}

function fresnelC( x ) {

  // can also be evaluated with hypergeometric1F2

  var m1 = hypergeometric1F1( .5, 1.5, mul( complex(0,pi/2), pow(x,2) ) );
  var m2 = hypergeometric1F1( .5, 1.5, mul( complex(0,-pi/2), pow(x,2) ) );

  var result = mul( x, add( m1, m2 ), .5 );

  if ( isComplex(x) ) return result;

  return result.re;

}


function expIntegralEi( x, tolerance=1e-10 ) {

  var useAsymptotic = 26;

  if ( isComplex(x) ) {

    if ( abs(x) > useAsymptotic ) {

      var s = complex(1);
      var p = complex(1);
      var i = 1;

      while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {
        p = mul( p, i, inv(x) );
        s = add( s, p );
        i++;
      }

      // combination of logarithms adds/subtracts complex(0,pi)
      var sign = x.im > 0 ? 1 : x.im < 0 ? -1 : 0;

      return add( mul( s, exp(x), inv(x) ), complex(0,sign*pi) );

    }

    // determined from pattern on test page
    var distanceScale = abs( sub(x,useAsymptotic) ) / useAsymptotic;
    var useArbitrary = distanceScale > 1;

    if ( useArbitrary ) {

      // use only decimals needed
      var n = 17 + Math.round( 10 * ( distanceScale - 1 ) );
      setPrecisionScale( n );

      var y = arbitrary( x );

      var s = arbitrary( complex(0) );
      var p = arbitrary( complex(1) );
      var i = arb1;

      while ( div(p.re,i) !== 0n || div(p.im,i) !== 0n ) {
        p = div( mul(p,y), i );
        s = add( s, div(p,i) );
        i = add( i, arb1 );
      }

      s = add( s, getConstant( 'eulerGamma' ), ln(y) );

      s = arbitrary( s );

      setPrecisionScale( defaultDecimals );

    } else {

      var s = complex(0);
      var p = complex(1);
      var i = 1;

      while ( Math.abs(p.re/i) > tolerance || Math.abs(p.im/i) > tolerance ) {
        p = mul( p, x, 1/i );
        s = add( s, div(p,i) );
        i++;
      }

      s = add( s, eulerGamma, log(x) );

    }

    // real on negative real axis, set phase explicitly rather than log combo
    if ( x.re < 0 && x.im === 0 ) s.im = 0;

    return s;

  } else {

    if ( x < 0 ) return expIntegralEi( complex(x) ).re;

    if ( Math.abs(x) > useAsymptotic ) {

      var s = 1;
      var p = 1;
      var i = 1;

      while ( Math.abs(p) > tolerance ) {
        p *= i / x;
        s += p;
        i++;
      }

      return s * Math.exp(x) / x;

    }

    var s = 0;
    var p = 1;
    var i = 1;

    while ( Math.abs(p/i) > tolerance ) {
      p *= x / i;
      s += p / i;
      i++;
    }

    return s + eulerGamma + Math.log(x);

  }

}

function logIntegral( x ) {

  return expIntegralEi( log(x) );

}

function sinIntegral( x ) {

  if ( isZero(x) ) return isComplex(x) ? complex(0) : 0;

  var ix = mul( complex(0,1), x );

  var result = mul( complex(0,.5), add( gamma(0,neg(ix)), neg(gamma(0,ix)),
                                        log(neg(ix)), neg(log(ix)) ) );

  if ( isComplex(x) ) return result;

  return result.re;

}

function cosIntegral( x ) {

  // complex for negative real argument

  var ix = mul( complex(0,1), x );

  return sub( log(x), mul( .5, add( gamma(0,neg(ix)), gamma(0,ix),
                                    log(neg(ix)), log(ix) ) ) );

}

function sinhIntegral( x ) {

  if ( isZero(x) ) return isComplex(x) ? complex(0) : 0;

  var result = mul( .5, add( gamma(0,x), neg(gamma(0,neg(x))), log(x), neg(log(neg(x))) ) );

  if ( isComplex(x) ) return result;

  return result.re;

}

function coshIntegral( x ) {

  // complex for negative real argument

  return mul( -.5, add( gamma(0,x), gamma(0,neg(x)), neg(log(x)), log(neg(x)) ) );

}

function expIntegralE( n, x ) {

  if ( isZero(n) ) return div( exp(neg(x)), x );

  if ( isZero(x) && ( n > 1 || n.re > 1 ) ) return inv( sub(n,1) );

  return mul( pow( x, sub(n,1) ), gamma( sub(1,n), x ) );

}

