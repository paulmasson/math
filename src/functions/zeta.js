
function zeta( x, tolerance=1e-10 ) {

  if ( isUnity(x) ) throw Error( 'Riemann zeta pole' );

  // functional equation dlmf.nist.gov/25.4.2 - connects both half planes
  if ( x < 0 || x.re < 0 )
    return mul( pow(2,x), pow(pi,sub(x,1)), sin( mul(pi/2,x) ), gamma( sub(1,x) ), zeta( sub(1,x) ) );

  // direct summation more accurate in right-hand plane
  var directSummation = 5;

  if ( x > directSummation || x.re > directSummation ) {

    if ( isComplex(x) ) {

      var s = complex(1);
      var p = complex(1);
      var i = 2;

      while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance*tolerance ) {
        p = pow( i, neg(x) );
        s = add( s, p );
        i++;
      }

      return s;

    } else {

      var s = 1;
      var p = 1;
      var i = 2;

      while ( Math.abs(p) > tolerance ) {
        p = 1 / i**x;
        s += p;
        i++;
      }

      return s;

    }

  }

  // Borwein, Efficient Algorithm

  var n = 14; // from error bound for tolerance

  if ( isComplex(x) && x.im !== 0 )
    n = Math.max( n, Math.ceil( log( 2 / abs(gamma(x)) / tolerance ) / log( 3 + sqrt(8) ) ) );

  var d = [ 1 ], prod = n;
  for ( var i = 1 ; i <= n ; i++ ) {
    d.push( d[i-1] + n * prod / factorial( 2*i ) * 4**i );
    // manually evaluate factorial( n+i-1 ) / factorial( n-i ) to avoid overflow
    prod *= (n+i)*(n-i);
  }

  if ( isComplex(x) ) {

    var s = summation( k => div( (-1)**k * ( d[k] - d[n] ), pow( k+1, x ) ), [0,n-1] );

    return div( div( s, -d[n] ), sub( 1, pow( 2, sub(1,x) ) ) );

  } else {

    var s = summation( k => (-1)**k * ( d[k] - d[n] ) / (k+1)**x, [0,n-1] );

    return -s / d[n] / ( 1 - 2**(1-x) );

  }

}

function dirichletEta( x ) { return mul( zeta(x), sub( 1, pow( 2, sub(1,x) ) ) ); }

function riemannXi( x ) {

  if ( isZero(x) || isUnity(x) ) return isComplex(x) ? complex(.5) : .5;

  var half = div( x, 2 );

  return mul( half, sub(x,1), pow( pi, neg(half) ), gamma(half), zeta(x) );

}


function bernoulli( n, x ) {

  if ( arguments.length === 2 ) {

    if ( isZero(n) ) return isComplex(n) || isComplex(x) ? complex(1) : 1;

    // avoid Hurwitz zeta parameter poles
    if (  isNegativeIntegerOrZero(x) ) return complexAverage( x => bernoulli(n,x), x );

    return mul( neg(n), hurwitzZeta( sub(1,n), x ) );

  }

  if ( Number.isInteger(n) && n >= 0 ) {

    if ( n === 0 ) return 1;

    if ( n === 1 ) return -.5;

    if ( n & 1 ) return 0;

    var m = n/2;
    if ( m <= bernoulli2nN.length )
      return arbitrary( div( bernoulli2nN[m], bernoulli2nD[m] ) );

    return -n * zeta(1-n);

  }

  if ( isPositiveIntegerOrZero(n) ) return complex( bernoulli( n.re ) );

  // generalized Bernoulli number
  console.log( 'Returning generalized Bernoulli number' );
  return bernoulli( n, 0 );

}

function harmonic( n ) {

  if ( !Number.isInteger(n) ) throw Error( 'Noninteger index for harmonic number' );

  if ( n > 1e3 ) return log(n) + eulerGamma + 1/2/n - 1/12/n**2;

  return summation( i => 1/i, [1,n] );

}


function hurwitzZeta( x, a, tolerance=1e-10 ) {

  if ( isUnity(x) ) throw Error( 'Hurwitz zeta pole' );

  if ( isNegativeIntegerOrZero(a) ) throw Error( 'Hurwitz zeta parameter pole' );

  if ( isComplex(x) || isComplex(a) ) {

    if ( !isComplex(x) ) x = complex(x);
    if ( !isComplex(a) ) a = complex(a);

    // direct summation more accurate than dlmf.nist.gov/25.11.4 for positive a

    if ( a.re < 0 ) {
      var m = -Math.floor(a.re);
      return add( hurwitzZeta(x,add(a,m)), summation( i => pow( add(a,i), neg(x) ), [0,m-1] ) );
    }

    // Johansson arxiv.org/abs/1309.2877

    var n = 15; // recommendation of Vepstas, Efficient Algorithm, p.12

    // Euler-Maclaurin has differences of large values in left-hand plane
    var useArbitrary = x.re < 0;

    if ( useArbitrary ) {

      setPrecisionScale( 20 - Math.round(x.re) );

      x = arbitrary(x), a = arbitrary(a);
      var arbN = arbitrary(n), arb3 = arbitrary(3);

      var S = 0n;
      for ( var i = 0 ; i < n ; i++ )
        S = add( S, pow( div( add(a,arbN), add(a,arbitrary(i)) ), x ) );

      var I = div( add(a,arbN), sub(x,arb1) );

      var p = div( mul( arb1/2n, x ), add(a,arbN) );
      var b = div( bernoulli2nN[1], bernoulli2nD[1] );
      var t = mul( b, p );
      var i = arb2;
      var j = 2;

      while ( p.re !== 0n || p.im !== 0n ) {
        if ( j === bernoulli2nN.length ) break;
        p = div( mul( p, add( x, 2n*i - arb2 ), add( x, 2n*i - arb3 ) ),
                 mul( 2n*i, 2n*i - arb1, add(a,arbN), add(a,arbN) ) );
        b = div( bernoulli2nN[j], bernoulli2nD[j] );
        t = add( t, mul( b, p ) );
        i += arb1;
        j++;
      }

      var T = add( arb1/2n, t );

      var result = arbitrary( mul( add( S, I, T ), pow( add(a,arbN), mul(-arb1,x) ) ) );

      resetPrecisionScale();

      return result;

    }

    var S = summation( i => pow( add(a,i), neg(x) ), [0,n-1] );

    var I = div( pow( add(a,n), sub(1,x) ), sub(x,1) );

    var p = mul( .5, x, inv(add(a,n)) );
    var t = mul( bernoulli(2), p );
    var i = 2;

    // converges rather quickly
    while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {
      p = div( mul( p, add( x, 2*i - 2 ), add( x, 2*i - 3 ) ),
               mul( 2*i * (2*i-1), add(a,n), add(a,n) ) );
      t = add( t, mul( bernoulli(2*i), p ) );
      i++;
    }

    var T = div( add( .5, t ), pow( add(a,n), x ) );

    return add( S, I, T );

  } else {

    if ( a < 0 ) return hurwitzZeta( x, complex(a) );

    // direct summation more accurate than dlmf.nist.gov/25.11.4

    // Euler-Maclaurin has differences of large values in left-hand plane
    // switch to different summation: dlmf.nist.gov/25.11.9

    var switchForms = -5;

    if ( x < switchForms ) {

      x = 1 - x;
      var t = Math.cos( pi*x/2 - 2*pi*a );
      var s = t;
      var i = 1;

      while ( Math.abs(t) > tolerance ) {
        i++;
        t = Math.cos( pi*x/2 - 2*i*pi*a ) / i**x;
        s += t;
      }

      return 2 * gamma(x) / (2*pi)**x * s;

    }

    // Johansson arxiv.org/abs/1309.2877

    var n = 15; // recommendation of Vepstas, Efficient Algorithm, p.12

    var S = summation( i => 1 / (a+i)**x, [0,n-1] );

    var I = (a+n)**(1-x) / (x-1);

    var p = x / 2 / (a+n);
    var t = bernoulli(2) * p;
    var i = 2;

    // converges rather quickly
    while ( Math.abs(p) > tolerance ) {
      p *= ( x + 2*i - 2 ) * ( x + 2*i - 3 ) / ( 2*i * (2*i-1) * (a+n)**2 );
      t += bernoulli(2*i) * p;
      i++;
    }

    var T = ( .5 + t ) / (a+n)**x;

    return S + I + T;

  }

}


function polylog( n, x, tolerance=1e-10 ) {

  if ( isEqualTo(x,1) ) return zeta(n);

  if ( isEqualTo(x,-1) ) return neg( dirichletEta(n) );

  if ( isEqualTo(n,1) ) return neg( log( sub(1,x) ) );

  if ( isEqualTo(n,0) ) return div( x, sub(1,x) );

  if ( isEqualTo(n,-1) ) return div( x, mul( sub(1,x), sub(1,x) ) );

  if ( abs(x) >= 1 ) {

    var twoPiI = complex(0,2*pi);

    if ( isPositiveInteger(n) ) {

      // Crandall, Note on Fast Polylogarithm Computation

      var t1 = mul( (-1)**n, polylog( n, inv(x) ) );

      var t2 = mul( div( pow(twoPiI,n), factorial(n) ), bernoulli( n, div(log(x),twoPiI) ) );

      var y = isComplex(x) ? x : complex(x); // just for test
      var t3 = y.im < 0 || ( y.im === 0 && y.re >= 1 ) ?
               mul( twoPiI, div( pow(log(x),n-1), factorial(n-1) ) ) : 0;      

      var result = neg( add( t1, t2, t3 ) );

      // real on negative real axis
      if ( !isComplex(x) && x < 0 ) return result.re;

      return result;

    }

    var v = sub(1,n);
    var I = complex(0,1);
    var L = div( log(neg(x)), twoPiI );

    var z1 = mul( pow(I,v), hurwitzZeta( v, add(.5,L) ) );
    var z2 = mul( pow(I,neg(v)), hurwitzZeta( v, sub(.5,L) ) );

    return mul( gamma(v), pow(2*pi,neg(v)), add(z1,z2) )

  }

  if ( isComplex(n) || isComplex(x) ) {

    var s = x;
    var p = complex(1);
    var i = 2;

    while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {
      p = div( pow(x,i), pow(i,n) );
      s = add( s, p );
      i++;
    }

    return s;

  } else {

    var s = x;
    var p = 1;
    var i = 2;

    while ( Math.abs(p) > tolerance ) {
      p = x**i / i**n;
      s += p;
      i++;
    }

    return s;

  }

}

