
function zeta( x, tolerance=1e-10 ) {

  // Borwein algorithm

  var n = 14; // from error bound for tolerance

  if ( isComplex(x) && x.im !== 0 )
    n = Math.max( n, Math.ceil( log( 2 / abs(gamma(x)) / tolerance ) / log( 3 + sqrt(8) ) ) );

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

  return summation( i => 1/i, [1,n] );

}


function hurwitzZeta( x, a, tolerance=1e-10 ) {

  if ( isComplex(x) || isComplex(a) ) {


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
    // swith to difference summation: dlmf.nist.gov/25.11.9

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
    var i = 1;

    // converges rather quickly
    while ( Math.abs(p) > tolerance ) {
      i++;
      p *= ( x + 2*i - 2 ) * ( x + 2*i - 3 ) / ( 2*i * (2*i-1) * (a+n)**2 )
      t += bernoulli(2*i) * p;
    }

    var T = ( .5 + t ) / (a+n)**x;

    return S + I + T;

  }

}

