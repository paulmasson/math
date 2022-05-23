
function factorial( n ) {

  if ( isComplex(n) ) {

    if ( isPositiveIntegerOrZero(n) ) return complex( factorial(n.re) );

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

  if ( isComplex(n) ) {

    if ( isPositiveIntegerOrZero(n) ) return complex( factorial2(n.re) );

    var p1 = pow( 2, div(n,2) );
    var p2 = pow( pi/2, div( sub( cos(mul(pi,n)), 1 ), 4 ) );
    var p3 = gamma( add( div(n,2) , 1 ) );

    return mul( p1, p2, p3 );

  }

  if ( n === 0 ) return 1;

  if ( isPositiveInteger(n) ) {

    // bitwise test for odd integer, upward recursion for possible caching
    var result = n & 1 ? 1 : 2;
    for ( var i = result + 2 ; i <= n ; i += 2 ) result *= i;
    return result;

  }

  return 2**(n/2) * (pi/2)**((cos(pi*n)-1)/4) * gamma( n/2+1 );

}

function binomial( n, m ) {

  if ( isComplex(n) || isComplex(m) )
    return div( factorial(n), mul( factorial( sub(n,m) ), factorial(m) ) );

  if ( Number.isInteger(m) && m < 0 && n >= 0 ) return 0;

  if ( Number.isInteger(n) && Number.isInteger(m) && n >= 0 && m > n ) return 0;

  return factorial(n) / factorial(n-m) / factorial(m);

}

function multinomial( n ) {

  if ( arguments.length === 1 ) return isComplex(n) ? complex(1) : 1;

  var s = 0, p = 1, nNeg = 0;

  for ( var i = 0 ; i < arguments.length ; i++ ) {

    var ni = arguments[i];
    s = add( s, ni );

    // catch negative integers for cancellation in ratio
    // reduce all to factorial(-1) divided by factor
    // extra overhead with pow() keeps complex in complex out

    if ( isNegativeInteger( ni ) ) {
      var k = sub( -1, ni );
      p = div( p, mul( chop(pow(-1,k)), factorial(k) ) );
      nNeg++;
    } else
      // standard factorial
      p = mul( p, factorial( arguments[i] ) );

  }

  if ( isNegativeInteger(s) ) {
    nNeg--;
    if ( nNeg > 0 ) return isComplex(s) ? complex(0) : 0;
    var k = sub( -1, s );
    return inv( mul( p, chop(pow(-1,k)), factorial(k) ) );
  }

  if ( nNeg > 0 ) return isComplex(s) ? complex(0) : 0;

  return div( factorial(s), p );

}

function pochhammer( x, n ) {

  var one = isArbitrary(x) ? arb1 : 1;

  if ( isZero(n) )
    if ( isComplex(x) || isComplex(n) ) return complex(one);
    else return one;

  if ( isComplex(n) ) {

    if ( isPositiveInteger(n) ) {
      var result = pochhammer( x, n.re );
      if ( isComplex(result) ) return result;
      else return complex(result);
    }

    if ( isArbitrary(x) && !isArbitrary(n) ) n = arbitrary(n);

    return div( gamma( add(x,n) ), gamma(x) );

  }

  if ( isPositiveInteger(n) ) {

    var result = x, current = one, count = 1;

    while ( count < n ) {
      result = mul( result, add( x, current ) );
      current = add( current, one );
      count++;
    }

    return result;

  }

  if ( isArbitrary(x) && !isArbitrary(n) ) n = arbitrary(n);

  return div( gamma( add(x,n) ), gamma(x) );

}


function logGamma( x ) {

  if ( isArbitrary(x) ) {

    var useAsymptotic = 10n * arb1;

    if ( isNegativeIntegerOrZero( arbitrary(x) ) ) throw Error( 'Gamma function pole' );

    if ( abs(x) < useAsymptotic ) {

      // slightly faster to get smaller integer near transition circle
      var xRe = isComplex(x) ? x.re : x;
      var n = Math.ceil( arbitrary(useAsymptotic) - arbitrary(xRe) );

      // ln(pochhammer(x,n))
      var lnPochhammer = ln(x), current = arb1, count = 1;

      while ( count < n ) {
        lnPochhammer = add( lnPochhammer, ln( add( x, current ) ) );
        current = add( current, arb1 );
        count++;
      }

      return sub( logGamma(add(x,arbitrary(n))), lnPochhammer );

    }

    if ( x < 0n ) x = complex(x);

    // reflection formula with non-Hare correction to imaginary part
    if ( x.re < 0n ) {

      // expand sine as exponentials for more accurate result
      var t, e, p, k = arb1;

      if ( abs(x.im) === 0 ) {

        t = neg( ln( sin( mul( x, onePi ) ) ) );

      } else if ( x.im < 0n ) {

        e = exp( mul( x, complex(0n,-twoPi) ) );
        p = complex( e.re, e.im );
        t = complex( e.re, e.im );
        while ( p.re !== 0n || p.im !== 0n ) {
          k += arb1;
          p = mul( p, e );
          t = add( t, div( p, k ) );
        }
        t = add( t, mul( x, complex(0n,-onePi) ), ln(arb2) );

      } else {

        e = exp( mul( x, complex(0n,twoPi) ) );
        p = complex( e.re, e.im );
        t = complex( e.re, e.im );
        while ( p.re !== 0n || p.im !== 0n ) {
          k += arb1;
          p = mul( p, e );
          t = add( t, div( p, k ) );
        }
        t = add( t, mul( x, complex(0n,onePi) ), complex(0n,-onePi), ln(arb2) );

      }

      return add( t, ln(onePi), neg( logGamma( sub(arb1,x) ) ), complex(0n,halfPi) );

    }

    // Johansson arxiv.org/abs/2109.08392

    var k = 1, K = arb1, y, p = 1n, s = 0n;

    if ( isComplex(x) ) {
      function compare() { return p.re !== 0n || p.im !== 0n; };
      y = complex( x.re, x.im );
    } else {
      function compare() { return p !== 0n; };
      y = x;
    }

    while ( compare() ) {

      if ( k === bernoulli2nN.length ) {
        console.log( 'Not enough Bernoulli numbers for logGamma' );
        break;
      }

      p = div( div( bernoulli2nN[k], bernoulli2nD[k] ),
               mul( K+arb1, K, y ) );
      s = add( s, p );

      y = mul( y, x, x );
      K += arb2;
      k++;

    }

    return add( mul( sub(x,arb1/2n), ln(x) ), neg(x), div(ln(twoPi),arb2), s );

  }

  // log of gamma less likely to overflow than gamma
  // Lanczos approximation as evaluated by Paul Godfrey

  var c = [ 57.1562356658629235, -59.5979603554754912, 14.1360979747417471,
            -.491913816097620199, .339946499848118887e-4, .465236289270485756e-4,
            -.983744753048795646e-4, .158088703224912494e-3, -.210264441724104883e-3,
            .217439618115212643e-3, -.164318106536763890e-3, .844182239838527433e-4,
            -.261908384015814087e-4, .368991826595316234e-5 ];

  if ( isNegativeIntegerOrZero(x) ) throw Error( 'Gamma function pole' );

  if ( isComplex(x) ) {

    // reflection formula with modified Hare correction to imaginary part
    if ( x.re < 0 ) {

      var logRatio = log( div( pi, sin( mul(pi,x) ) ) );
      // rounding errors can lead to wrong side of branch point
      if ( isNegativeIntegerOrZero( ( x.re + .5 ) / 2 ) )
        logRatio.im = pi * ( x.im > 0 ? 1 : -1 ); // avoid Math.sign(0) = 0

      var t = sub( logRatio, logGamma( sub(1,x) ) );
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

    if ( x < 0 ) return logGamma( complex(x) ); 

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

      var result = neg( expIntegralEi( neg(y), true ) );

      // complex on negative real axis
      if ( y < 0 || y.re < 0 && y.im === 0 ) result = sub( result, complex(0,pi) );

      return result;

    }

    // dlmf.nist.gov/8.4.15
    if ( isNegativeInteger(x) ) {

      var n = isComplex(x) ? -x.re : -x;
      var s = inv(y), p = s; // mul returns new object

      for ( var k = 1 ; k < n ; k++ ) {
        p = mul( p, -k, inv(y) );
        s = add( s, p );
      }

      var t = mul( exp(neg(y)), s );

      // dlmf.nist.gov/8.4.4
      var result = sub( neg( expIntegralEi( neg(y), true, 1e-14 ) ), t )

      // complex on negative real axis
      if ( y < 0 || y.re < 0 && y.im === 0 ) result = sub( result, complex(0,pi) );

      result = mul( (-1)**n/factorial(n), result );

      if ( isComplex(x) && !isComplex(result) ) return complex(result); // complex in, complex out

      return result;

    }

    return mul( exp(neg(y)), hypergeometricU( sub(1,x), sub(1,x), y ) );

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


function expIntegralEi( x, adjustImForGamma=false, tolerance=1e-10 ) {

  // direct series hangs at -1 and +i/-i
  // complex average until investigate further

  if ( isUnity(neg(x)) ) return complexAverage( x => expIntegralEi(x,adjustImForGamma,tolerance), x );

  var ix = mul( x, complex(0,1) );
  if ( isUnity(ix) || isUnity(neg(ix)) )
    return complexAverage( x => expIntegralEi(x,adjustImForGamma,tolerance), x, complex(0,1e-5) );

  var useAsymptotic = 26;

  if ( isComplex(x) ) {

    if ( abs(x) > useAsymptotic ) {

      var s = complex(1);
      var p = complex(1), pLast = p; // mul returns new object
      var i = 1;

      while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {

        p = mul( p, i, inv(x) );
        if ( abs(p) > abs(pLast) ) break;

        s = add( s, p );
        i++;
        pLast = p;

      }

      if ( adjustImForGamma )

        return mul( s, exp(x), inv(x) );

      else {

        // combination of logarithms adds/subtracts complex(0,pi)
        var sign = x.im > 0 ? 1 : x.im < 0 ? -1 : 0;

        return add( mul( s, exp(x), inv(x) ), complex(0,sign*pi) );

      }

    }

    // determined from pattern on test page
    var distanceScale = abs( sub(x,useAsymptotic) ) / useAsymptotic;
    var useArbitrary = distanceScale > 1;

    // arbitrary precision series is unstable around the origin but not needed

    if ( useArbitrary && abs(x) > 1.5 ) {

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

      if ( adjustImForGamma ) {
        var sign = x.im > 0n ? -1n : x.im < 0n ? 1n : 0n;
        s = add( s, complex(0n,sign*onePi) );
      }

      s = arbitrary( s );

      resetPrecisionScale();

      // adjust phase on positive imaginary axis according to test
      if ( x.im > 0 && x.re === 0 ) s.im += pi;

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

      if ( adjustImForGamma ) {
        var sign = x.im > 0 ? -1 : x.im < 0 ? 1 : 0;
        s = add( s, complex(0,sign*pi) );
      }

    }

    // real on negative real axis, set phase explicitly rather than log combo
    if ( x.re < 0 && x.im === 0 ) s.im = 0;

    return s;

  } else {

    if ( x < 0 ) return expIntegralEi( complex(x) ).re;

    if ( Math.abs(x) > useAsymptotic ) {

      var s = 1;
      var p = 1, pLast = p;
      var i = 1;

      while ( Math.abs(p) > tolerance ) {

        p *= i / x;
        if ( Math.abs(p) > Math.abs(pLast) ) break;

        s += p;
        i++;
        pLast = p;

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

  var result = sub( log(x), mul( .5, add( gamma(0,neg(ix)), gamma(0,ix),
                                     log(neg(ix)), log(ix) ) ) );

  // real for positive real argument
  if ( x > 0 ) return result.re;
  if ( x.re > 0 && x.im === 0 ) result.im = chop( result.im );

  return result;

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

  var p = pow( x, sub(n,1) );

  // real on negative real axis for integer powers
  if ( isInteger(n) && x.re < 0 ) p.im = chop( p.im );

  return mul( p, gamma( sub(1,n), x ) );

}

