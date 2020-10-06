
function exp( x ) {

  if ( isComplex(x) )

    return { re: Math.exp(x.re) * Math.cos(x.im),
             im: Math.exp(x.re) * Math.sin(x.im) };

  return Math.exp(x);

}


function log( x, base ) {

  if ( isComplex(x) ) {

    if ( isComplex(base) ) return div( log(x), log(base) );

    return { re: log( abs(x), base ), im: log( Math.E, base ) * arg(x) };

  }

  if ( x < 0 ) return log( complex(x), base );

  if ( base === undefined ) return Math.log(x);

  return Math.log(x) / Math.log(base);

}

function ln( x ) {

  // Brent, Modern Computer Arithmetic, second AGM algorithm

  function arbitraryAGM( x, y ) {

    var t, u, arb2 = arbitrary(2);

    while( x !== y ) {
      t = x, u = y;
      x = div( t + u, arb2 );
      y = sqrt( mul(t,u) );
    }

    return x;

  }

  function arbitraryTheta2( x ) {

    var p = mul( arbitrary(2), x );
    var s = p;
    var i = 1;

    while ( p !== 0n ) {
      for ( var j = 0 ; j < 8*i ; j++ ) p = mul( p, x );
      s = s + p;
      i++;
    }

    return s;

  }

  function arbitraryTheta3( x ) {

    var p = arbitrary(2);
    var s = arbitrary(1);
    var i = 1;

    while ( p !== 0n ) {
      for ( var j = 0 ; j < 4*(2*i-1) ; j++ ) p = mul( p, x );
      s = s + p;
      i++;
    }

    return s;

  }

  if ( isArbitrary(x) ) {

    if ( x < 0n ) throw Error( 'Complex arbitrary logarithm not yet supported' );

    var arb1 = arbitrary(1);

    if ( x === arb1 ) return 0n;

    if ( x < arb1 ) return -ln( div( arb1, x ) );

    x = div( arb1, x );

    var t2 = arbitraryTheta2(x);
    var t3 = arbitraryTheta3(x);

    return div( getConstant('pi'), mul( arbitrary(4), arbitraryAGM( mul(t2,t2), mul(t3,t3) ) ) );

  }

  return log(x);

}


function lambertW( k, x, tolerance=1e-10 ) {

  if ( arguments.length === 1 ) {
    x = k;
    k = 0;
  }

  if ( Math.abs( x + Math.exp(-1) ) < tolerance ) return -1;

  // inversion by root finding

  switch ( k ) {

    case 0:

      if ( x < -Math.exp(-1) ) throw Error( 'Unsupported lambertW argument' );

      return findRoot( w => w * Math.exp(w) - x, [-1,1000], { tolerance: tolerance } );

    case -1:

      if ( x < -Math.exp(-1) || x > 0 ) throw Error( 'Unsupported lambertW argument' );

      return findRoot( w => w * Math.exp(w) - x, [-1000,-1], { tolerance: tolerance } );

    default:

      throw Error( 'Unsupported lambertW index' );

  }

}

