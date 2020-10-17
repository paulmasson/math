
function zeta( x, tolerance=1e-10 ) {

  // direct summation fast in right-hand plane
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

  // Borwein p.3 simplified
  if ( isComplex(x) && x.im !== 0 )
    n = Math.max( n, Math.ceil( log( 2 / abs(gamma(x)) / tolerance ) / 1.5 ) );

  var d = [ 1 ];
  for ( var i = 1 ; i <= n ; i++ )
    // order of multiplication reduces overflow, but factorial overflows at 171
    d.push( d[i-1] + n * factorial( n+i-1 ) / factorial( n-i ) / factorial( 2*i ) * 4**i );

  if ( isComplex(x) ) {

    // functional equation dlmf.nist.gov/25.4.2
    if ( x.re < 0 )
      return mul( pow(2,x), pow(pi,sub(x,1)), sin( mul(pi/2,x) ), gamma( sub(1,x) ), zeta( sub(1,x) ) );

    var s = summation( k => div( (-1)**k * ( d[k] - d[n] ), pow( k+1, x ) ), [0,n-1] );

    return div( div( s, -d[n] ), sub( 1, pow( 2, sub(1,x) ) ) );

  } else {

    // functional equation dlmf.nist.gov/25.4.2
    if ( x < 0 ) return 2**x * pi**(x-1) * sin(pi*x/2) * gamma(1-x) * zeta(1-x);

    var s = summation( k => (-1)**k * ( d[k] - d[n] ) / (k+1)**x, [0,n-1] );

    return -s / d[n] / ( 1 - 2**(1-x) );

  }

}

function dirichletEta( x ) { return mul( zeta(x), sub( 1, pow( 2, sub(1,x) ) ) ); }


function bernoulli( n ) {

  if ( !Number.isInteger(n) ) throw Error( 'Noninteger argument for Bernoulli number' );

  if ( n < 0 ) throw Error( 'Unsupported argument for Bernoulli number' );

  if ( n === 0 ) return 1;

  if ( n === 1 ) return -.5;

  if ( n & 1 ) return 0;

  return (-1)**(n+1) * n * zeta(-n+1);

}

function harmonic( n ) {

  if ( !Number.isInteger(n) ) throw Error( 'Noninteger argument for harmonic number' );

  if ( n > 1e3 ) return log(n) + eulerGamma + 1/2/n - 1/12/n**2;

  return summation( i => 1/i, [1,n] );

}


function hurwitzZeta( x, a, tolerance=1e-10 ) {

  if ( isComplex(x) || isComplex(a) ) {

    if ( !isComplex(x) ) x = complex(x);
    if ( !isComplex(a) ) a = complex(a);

    if ( x.re === 1 && x.im === 0 ) throw Error( 'Hurwitz zeta pole' );

    // dlmf.nist.gov/25.11.4

    if ( a.re > 1 ) {
      var m = Math.floor(a.re);
      a = sub( a, m );
      return sub( hurwitzZeta(x,a), summation( i => pow( add(a,i), neg(x) ), [0,m-1] ) );
    }

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

      x = arbitrary(x), a = arbitrary(a), arbN = arbitrary(n), arb3 = arbitrary(3);

      var S = 0n;
      for ( var i = 0 ; i < n ; i++ )
        S = add( S, pow( div( add(a,arbN), add(a,arbitrary(i)) ), x ) );

      var I = div( add(a,arbN), sub(x,arb1) );

      var p = div( mul( arb1/2n, x ), add(a,arbN) );
      var b = div( arbitrary(bernoulli2nN[1]), arbitrary(bernoulli2nD[1]) );
      var t = mul( b, p );
      var i = arb2;
      var j = 2;

      while ( p.re !== 0n || p.im !== 0n ) {
        if ( j === bernoulli2nN.length ) break;
        p = div( mul( p, add( x, 2n*i - arb2 ), add( x, 2n*i - arb3 ) ),
                 mul( 2n*i, 2n*i - arb1, add(a,arbN), add(a,arbN) ) );
        b = div( arbitrary(bernoulli2nN[j]), arbitrary(bernoulli2nD[j]) );
        t = add( t, mul( b, p ) );
        i += arb1;
        j++;
      }

      var T = add( arb1/2n, t );

      return arbitrary( mul( add( S, I, T ), pow( add(a,arbN), mul(-arb1,x) ) ) );

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

    if ( x === 1 ) throw Error( 'Hurwitz zeta pole' );

    // dlmf.nist.gov/25.11.4

    if ( a > 1 ) {
      var m = Math.floor(a);
      a -= m;
      return hurwitzZeta(x,a) - summation( i => 1 / (a+i)**x, [0,m-1] );
    }

    if ( a < 0 ) return hurwitzZeta( x, complex(a) );

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

