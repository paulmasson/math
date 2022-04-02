
function hypergeometric0F1( a, x, tolerance=1e-10 ) {

  var useAsymptotic = 100;

  if ( isNegativeIntegerOrZero(a) ) throw Error( 'Hypergeometric function pole' );

  if ( isComplex(a) || isComplex(x) ) {

    // asymptotic form as per Johansson arxiv.org/abs/1606.06977
    if ( abs(x) > useAsymptotic ) {

      // transform variables for convenience
      var b = sub( mul(2,a), 1 );
      a = sub( a, .5 );
      x = mul( 4, sqrt(x) );

      // copied from hypergeometric1F1
      var t1 = mul( gamma(b), pow( neg(x), neg(a) ), inv( gamma(sub(b,a)) ) );
      t1 = mul( t1, hypergeometric2F0( a, add(a,neg(b),1), div(-1,x) ) );

      var t2 = mul( gamma(b), pow( x, sub(a,b) ), exp(x), inv( gamma(a) ) );
      t2 = mul( t2, hypergeometric2F0( sub(b,a), sub(1,a), div(1,x) ) );

      return mul( exp( div(x,-2) ), add( t1, t2 ) );

    }

    var s = complex(1);
    var p = complex(1);
    var i = 1;

    while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {
      p = mul( p, x, inv(a), 1/i );
      s = add( s, p );
      a = add( a, 1 );
      i++;
    }

    return s;

  } else {

    // asymptotic form is complex
    if ( Math.abs(x) > useAsymptotic ) return hypergeometric0F1( a, complex(x) ).re;

    var s = 1;
    var p = 1;
    var i = 1;

    while ( Math.abs(p) > tolerance ) {
      p *= x / a / i;
      s += p;
      a++;
      i++;
    }

    return s;

  }

}


function hypergeometric1F1( a, b, x, tolerance=1e-10 ) {

  if ( isEqualTo(a,b) ) return exp(x);

  if ( isArbitrary(x) ) {

    if ( !isArbitrary(a) ) a = arbitrary(a);
    if ( !isArbitrary(b) ) b = arbitrary(b);

    if ( isNegativeIntegerOrZero(arbitrary(b)) ) throw Error( 'Hypergeometric function pole' );

    // asymptotic forms not yet implemented

    if ( isComplex(x) ) {

      // Kummer transformation
      if ( x.re < 0n ) return mul( exp(x), hypergeometric1F1( sub(b,a), b, neg(x) ) );

      var s = complex(arb1);
      var p = complex(arb1);
      var i = arb1;

      while ( abs(p.re) !== 0n || abs(p.im) !== 0n ) {
        p = mul( p, x, a, inv(b), inv(i) );
        s = add( s, p );
        a = add( a, arb1 );
        b = add( b, arb1 );
        i += arb1;
      }

      return s;

    } else {

      // Kummer transformation
      if ( x < 0n ) return mul( exp(x), hypergeometric1F1( b-a, b, -x ) );

      var s = arb1;
      var p = arb1;
      var i = arb1;

      while ( abs(p) !== 0n ) {
        p = mul( p, x, a, inv(b), inv(i) );
        s += p;
        a += arb1;
        b += arb1;
        i += arb1;
      }

      return s;

    }

  }

  var useAsymptotic = 30;

  if ( isNegativeIntegerOrZero(b) ) throw Error( 'Hypergeometric function pole' );

  if ( isComplex(a) || isComplex(b) || isComplex(x) ) {

    if ( isZero(a) ) return complex(1);

    if ( !isComplex(x) ) x = complex(x);

    // Kummer transformation
    if ( x.re < 0 ) return mul( exp(x), hypergeometric1F1( sub(b,a), b, neg(x) ) );

    // asymptotic form as per Johansson arxiv.org/abs/1606.06977
    if ( abs(x) > useAsymptotic ) {

      if ( isZero(a) || isNegativeIntegerOrZero(sub(b,a)) )
        return complexAverage( a => hypergeometric1F1(a,b,x), a );

      var t1 = mul( gamma(b), pow( neg(x), neg(a) ), inv( gamma(sub(b,a)) ),
                    hypergeometric2F0( a, add(a,neg(b),1), div(-1,x) ) );

      var t2 = mul( gamma(b), pow( x, sub(a,b) ), exp(x), inv( gamma(a) ),
                    hypergeometric2F0( sub(b,a), sub(1,a), div(1,x) ) );

      return add( t1, t2 );

    }

    var s = complex(1);
    var p = complex(1);
    var i = 1;

    while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {
      p = mul( p, x, a, inv(b), 1/i );
      s = add( s, p );
      a = add( a, 1 );
      b = add( b, 1 );
      i++;
    }

    return s;

  } else {

    if ( isZero(a) ) return 1;

    // Kummer transformation
    if ( x < 0 ) return exp(x) * hypergeometric1F1( b-a, b, -x );

    // asymptotic form is complex
    if ( Math.abs(x) > useAsymptotic ) return hypergeometric1F1( a, b, complex(x) ).re;

    var s = 1;
    var p = 1;
    var i = 1;

    while ( Math.abs(p) > tolerance ) {
      p *= x * a / b / i;
      s += p;
      a++;
      b++;
      i++;
    }

    return s;

  }

}


function hypergeometricU( a, b, x ) {

  if ( isEqualTo( add(a,1), b ) ) return pow( x, neg(a) );

  if ( isZero(a) )
    if ( isComplex(a) || isComplex(b) || isComplex(x) ) return complex(1);
    else return 1;

  var useAsymptotic = 25;

  // asymptotic form as per Johansson arxiv.org/abs/1606.06977
  if ( abs(x) > useAsymptotic ) {

    return mul( pow( x, neg(a) ), hypergeometric2F0( a, add(a,neg(b),1), neg(inv(x)) ) );

  }

  // special case through independent expIntegralEi for speed
  if ( isUnity(a) && isUnity(b) ) return mul( exp(x), gamma(0,x) );

  if ( isUnity(b) ) return complexAverage( b => hypergeometricU(a,b,x), b );

  if ( isPositiveInteger(a) && isEqualTo(a,b) ) return mul( exp(x), gamma(sub(1,a),x) );

  var t1 = mul( gamma(sub(b,1)), inv( gamma(a) ), pow( x, sub(1,b) ),
                hypergeometric1F1( add(a,neg(b),1), sub(2,b), x ) );

  var t2 = mul( gamma(sub(1,b)), inv( gamma(add(a,neg(b),1)) ), hypergeometric1F1( a, b, x ) );

  var max = Math.max( abs(t1), abs(t2) );
  if ( max < 100 ) return add( t1, t2 );

  var extraDecimals = +max.toExponential().split( 'e' )[1];
  setPrecisionScale( defaultDecimals + extraDecimals );

  a = arbitrary(a);
  b = arbitrary(b);
  x = arbitrary(x);

  t1 = mul( gamma(sub(b,arb1)), inv( gamma(a) ), pow( x, sub(arb1,b) ),
            hypergeometric1F1( add(a,neg(b),arb1), sub(arb2,b), x ) );

  t2 = mul( gamma(sub(arb1,b)), inv( gamma(add(a,neg(b),arb1)) ), hypergeometric1F1( a, b, x ) );

  var result = arbitrary( add( t1, t2 ) );

  resetPrecisionScale();

  return result;

}

function whittakerM( k, m, x ) {

  return mul( exp( mul(-.5,x) ), pow( x, add(m,.5) ),
              hypergeometric1F1( add(m,neg(k),.5), add(mul(2,m),1), x ) );

}

function whittakerW( k, m, x ) {

  return mul( exp( mul(-.5,x) ), pow( x, add(m,.5) ),
              hypergeometricU( add(m,neg(k),.5), add(mul(2,m),1), x ) );

}


function hypergeometric2F0( a, b, x, tolerance=1e-10 ) {

  var terms = 50;

  if ( isComplex(a) || isComplex(b) || isComplex(x) ) {

    var s = complex(1);
    var p = complex(1), pLast = p; // mul returns new object
    var converging = false; // first few terms can be larger than unity
    var i = 1;

    while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {

      p = mul( p, x, a, b, 1/i );

      if ( abs(p) > abs(pLast) && converging ) break; // prevent runaway sum
      if ( abs(p) < abs(pLast) ) converging = true;
      if ( i > terms ) throw Error( 'Not converging after ' + terms + ' terms' );

      s = add( s, p );
      a = add( a, 1 );
      b = add( b, 1 );
      i++;
      pLast = p;

    }

    return s;

  } else {

    var s = 1;
    var p = 1, pLast = p;
    var converging = false; // first few terms can be larger than unity
    var i = 1;

    while ( Math.abs(p) > tolerance ) {

      p *= x * a * b / i;

      if ( Math.abs(p) > Math.abs(pLast) && converging ) break; // prevent runaway sum
      if ( Math.abs(p) < Math.abs(pLast) ) converging = true;
      if ( i > terms ) throw Error( 'Not converging after ' + terms + ' terms' );

      s += p;
      a++;
      b++;
      i++;
      pLast = p;

    }

    return s;

  }

}


function hypergeometric2F1( a, b, c, x, tolerance=1e-10 ) {

  if ( isEqualTo(a,c) ) return pow( sub(1,x), neg(b) );

  if ( isEqualTo(b,c) ) return pow( sub(1,x), neg(a) );

  if ( isNegativeIntegerOrZero(c) ) throw Error( 'Hypergeometric function pole' );

  if ( isComplex(a) || isComplex(b) || isComplex(c) || isComplex(x) ) {

    if ( isZero(a) || isZero(b) ) return complex(1);

    // choose smallest absolute value of transformed argument
    // transformations from Abramowitz & Stegun p.559
    // fewer operations compared to dlmf.nist.gov/15.8

    var absArray = [ abs(x), abs(div(x,sub(x,1))), abs(sub(1,x)),
                     abs(inv(x)), abs(inv(sub(1,x))), abs(sub(1,inv(x))) ];

    var index = absArray.indexOf( Math.min.apply( null, absArray ) );

    switch( index ) {

      case 0:

        break;

      case 1:

        return mul( pow( sub(1,x), neg(a) ), hypergeometric2F1( a, sub(c,b), c, div(x,sub(x,1)) ) );

      case 2:

        if ( isInteger(sub(c,add(a,b))) || isNegativeIntegerOrZero(sub(c,a)) )
          return complexAverage( a => hypergeometric2F1(a,b,c,x), a );

        if ( isNegativeIntegerOrZero(sub(c,b)) )
          return complexAverage( b => hypergeometric2F1(a,b,c,x), b );

        var t1 = mul( gamma(c), gamma( sub(c,add(a,b)) ), 
                      inv( gamma(sub(c,a)) ), inv( gamma(sub(c,b)) ),
                      hypergeometric2F1( a, b, add(a,b,neg(c),1), sub(1,x) ) );

        var t2 = mul( pow( sub(1,x), sub(c,add(a,b)) ),
                      gamma(c), gamma( sub(add(a,b),c) ), inv( gamma(a) ), inv( gamma(b) ),
                      hypergeometric2F1( sub(c,a), sub(c,b), add(c,neg(a),neg(b),1), sub(1,x) ) );

        return add( t1, t2 );

      case 3:

        if ( isInteger(sub(a,b)) || isNegativeIntegerOrZero(sub(c,a)) )
          return complexAverage( a => hypergeometric2F1(a,b,c,x), a );

        if ( isNegativeIntegerOrZero(sub(c,b)) )
          return complexAverage( b => hypergeometric2F1(a,b,c,x), b );

        var t1 = mul( gamma(c), gamma(sub(b,a)), inv( gamma(b) ),
                      inv( gamma(sub(c,a)) ), pow( neg(x), neg(a) ),
                      hypergeometric2F1( a, add(1,neg(c),a), add(1,neg(b),a), inv(x) ) );

        var t2 = mul( gamma(c), gamma(sub(a,b)), inv( gamma(a) ),
                      inv( gamma(sub(c,b)) ), pow( neg(x), neg(b) ),
                      hypergeometric2F1( b, add(1,neg(c),b), add(1,neg(a),b), inv(x) ) );

        return add( t1, t2 );

      case 4:

        if ( isInteger(sub(a,b)) || isNegativeIntegerOrZero(sub(c,a)) )
          return complexAverage( a => hypergeometric2F1(a,b,c,x), a );

        if ( isNegativeIntegerOrZero(sub(c,b)) )
          return complexAverage( b => hypergeometric2F1(a,b,c,x), b );

        var t1 = mul( pow( sub(1,x), neg(a) ), gamma(c), gamma(sub(b,a)),
                      inv( gamma(b) ), inv( gamma(sub(c,a)) ),
                      hypergeometric2F1( a, sub(c,b), add(a,neg(b),1), inv(sub(1,x)) ) );

        var t2 = mul( pow( sub(1,x), neg(b) ), gamma(c), gamma(sub(a,b)),
                      inv( gamma(a) ), inv( gamma(sub(c,b)) ),
                      hypergeometric2F1( b, sub(c,a), add(b,neg(a),1), inv(sub(1,x)) ) );

        return add( t1, t2 );

      case 5:

        if ( isInteger(sub(c,add(a,b))) || isNegativeIntegerOrZero(sub(c,a)) )
          return complexAverage( a => hypergeometric2F1(a,b,c,x), a );

        if ( isNegativeIntegerOrZero(sub(c,b)) )
          return complexAverage( b => hypergeometric2F1(a,b,c,x), b );

        var t1 = mul( gamma(c), gamma( sub(c,add(a,b)) ), inv( gamma(sub(c,a)) ),
                      inv( gamma(sub(c,b)) ), pow( x, neg(a) ),
                      hypergeometric2F1( a, add(a,neg(c),1), add(a,b,neg(c),1), sub(1,inv(x)) ) );

        var t2 = mul( gamma(c), gamma( sub(add(a,b),c) ), inv( gamma(a) ), inv( gamma(b) ),
                      pow( sub(1,x), sub(c,add(a,b)) ), pow( x, sub(a,c) ),
                      hypergeometric2F1( sub(c,a), sub(1,a), add(c,neg(a),neg(b),1), sub(1,inv(x)) ) );

        return add( t1, t2 );

    }

    var s = complex(1);
    var p = complex(1);
    var i = 1;

    while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {
      p = mul( p, x, a, b, inv(c), 1/i );
      s = add( s, p );
      a = add( a, 1 );
      b = add( b, 1 );
      c = add( c, 1 );
      i++;
    }

    return s;

  } else {

    if ( isZero(a) || isZero(b) ) return 1;

    // transformation from Abramowitz & Stegun p.559
    if ( x < -1 ) {

      var t1 = gamma(c) * gamma(b-a) / gamma(b) / gamma(c-a)
                 * (-x)**(-a) * hypergeometric2F1( a, 1-c+a, 1-b+a, 1/x );
      var t2 = gamma(c) * gamma(a-b) / gamma(a) / gamma(c-b)
                 * (-x)**(-b) * hypergeometric2F1( b, 1-c+b, 1-a+b, 1/x );

      return t1 + t2;

    }

    if ( x === -1 ) return hypergeometric2F1( a, b, c, complex(x) ).re;

    if ( x === 1 )
      if ( c - a - b > 0 ) return gamma(c) * gamma(c-a-b) / gamma(c-a) / gamma(c-b);
      else throw Error( 'Divergent Gauss hypergeometric function' );

    if ( x > 1 ) return hypergeometric2F1( a, b, c, complex(x) );

    var s = 1;
    var p = 1;
    var i = 1;

    while ( Math.abs(p) > tolerance ) {
      p *= x * a * b / c / i;
      s += p;
      a++;
      b++;
      c++;
      i++;
    }

    return s;

  }

}


function hypergeometric1F2( a, b, c, x, tolerance=1e-10 ) {

  var useAsymptotic = 150;

  if ( isComplex(a) || isComplex(b) || isComplex(c) || isComplex(x) ) {

    // functions.wolfram.com/HypergeometricFunctions/Hypergeometric1F2/06/02/03/0002/

    if ( abs(x) > useAsymptotic ) {

      var aNegbNegc = add( a, neg(b), neg(c) );

      var p = div( add( aNegbNegc, 1/2 ), 2 );

      var ck = [ 1 ];

      ck[1] = add( mul( add(mul(3,a),b,c,-2), aNegbNegc, 1/2 ), mul(2,b,c), -3/8 );

      var toSquare = add( mul(b,c), mul( aNegbNegc, add(mul(3,a),b,c,-2), 1/4 ), -3/16 );

      ck[2] = add( mul( 2, toSquare, toSquare ),
                   neg( mul( sub(mul(2,a),3), b, c ) ),
                   mul( aNegbNegc, add( mul(-8,mul(a,a)), mul(11,a), b, c, -2 ), 1/4 ),
                   -3/16 );

      var plusI = complex(0,1), minusI = complex(0,-1);

      function w( k ) { return mul( 1/2**k, ck[k], pow(neg(x),-k/2) ); }

      var w1 = w(1), w2 = w(2);

      var s1 = add( 1, mul(minusI,w1), neg(w2) );
      var s2 = add( 1, mul(plusI,w1), neg(w2) );

      var k = 2, wk = w(k), powPlusI = powMinusI = -1;

      // wk not always complex
      while ( abs(wk) > tolerance ) {

        k++;
        powPlusI = mul( powPlusI, plusI );
        powMinusI = mul( powMinusI, minusI );

        var t1 = mul( add( 3*k**2, mul( 2*k, add(mul(-3,a),b,c,-2) ),
                               mul(3,mul(a,a)), neg(mul(sub(b,c),sub(b,c))),
                               neg(mul(2,a,add(b,c,-2))), 1/4 ),
                      ck[k-1] );

        var t2 = mul( add(k,neg(a),b,neg(c),-1/2), add(k,neg(a),neg(b),c,-1/2),
                      add(k,neg(a),b,c,-5/2), ck[k-2] );

        ck[k] = div( sub( t1, t2 ), 2*k );

        wk = w(k);

        s1 = add( s1, mul( powMinusI, wk ) );
        s2 = add( s2, mul( powPlusI, wk ) );

      }

      var exponent = add( mul(pi,p), mul(2,sqrt(neg(x))) );

      var u1 = exp( mul( plusI, exponent ) );
      var u2 = exp( mul( minusI, exponent ) );

      var s = add( mul(u1,s1), mul(u2,s2) );

      var t1 = mul( 1/(2*sqrt(pi)), inv(gamma(a)), pow(neg(x),p), s );

      var t2 = mul( inv(gamma(sub(b,a))), inv(gamma(sub(c,a))), pow(neg(x),neg(a)),
                    hypergeometricSeries( [ a, add(a,neg(b),1), add(a,neg(c),1) ], [], inv(x), true ) );

      return mul( gamma(b), gamma(c), add( t1, t2 ) );

    }

    return hypergeometricSeries( [a], [b,c], x, true );

  } else {

    // asymptotic form is complex
    if ( Math.abs(x) > useAsymptotic ) return hypergeometric1F2( a, b, c, complex(x) ).re;

    return hypergeometricSeries( [a], [b,c], x );

  }

}


// convenience function for less-used hypergeometrics
// accessing array slower than local variables
// for loops faster than forEach or reduce

function hypergeometricSeries( A, B, x, complexArguments=false, tolerance=1e-10 ) {

  if ( complexArguments ) {

    var s = complex(1);
    var p = complex(1);
    var i = 1;

    while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {

      for ( var j = 0 ; j < A.length ; j++ ) {
        p = mul( p, A[j] );
        A[j] = add( A[j], 1 );
      }

      for ( var j = 0 ; j < B.length ; j++ ) {
        p = div( p, B[j] );
        B[j] = add( B[j], 1 );
      }

      p = mul( p, x, 1/i );
      s = add( s, p );
      i++;

    }

    return s;

  } else {

    var s = 1;
    var p = 1;
    var i = 1;

    while ( Math.abs(p) > tolerance ) {

      for ( var j = 0 ; j < A.length ; j++ ) {
        p *= A[j];
        A[j]++;
      }

      for ( var j = 0 ; j < B.length ; j++ ) {
        p /= B[j];
        B[j]++;
      }

      p *= x / i;
      s += p;
      i++;

    }

    return s;

  }

}


function hypergeometricPFQ( A, B, x ) {

  // dlmf.nist.gov/16.11 for general transformations

  // for B.length > A.length terms can get very large
  // roundoff errors even for formally convergent series

  if ( abs(x) > 1 ) throw Error( 'General hypergeometric argument currently restricted' );

  // check for complex parameters
  var cp = false;
  A.forEach( a => cp = cp || isComplex(a) );
  B.forEach( b => cp = cp || isComplex(b) );

  if ( cp || isComplex(x) ) return hypergeometricSeries( A, B, x, true );

  return hypergeometricSeries( A, B, x );

}

