
function exp( x ) {

  if ( isArbitrary(x) ) {

    var ln10 = ln(arbitrary(10));
    var m = Math.trunc( arbitrary( div( x, ln10 ) ) );
    x = x - mul( arbitrary(m), ln10 );

    // direct sum faster than function inversion
    var arb1 = arbitrary(1);
    var s = arb1;
    var p = arb1;
    var i = arb1;

    while ( p !== 0n ) {
      p = div( mul( p, x ), i );
      s += p;
      i += arb1;
    }

    return mul( s, arbitrary( Number('1e'+m) ) );

  }

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

    if ( isComplex(x) ) {

      var maxIter = 10, i = 0;

      while( x.re !== y.re || x.im !== y.im ) {
        t = x, u = y;
        x = div( add( t, u ), arb2 );
        y = sqrt( mul( t, u ) );
        i++;
        if ( i > maxIter ) break; // convergence on complex plane not assured...
      }

    } else {

      while( x !== y ) {
        t = x, u = y;
        x = div( t + u, arb2 );
        y = sqrt( mul( t, u ) );
      }

    }

    return x;

  }

  function arbitraryTheta2( x ) {

    var p = mul( arbitrary(2), x );
    var s = p;
    var i = 1;

    if ( isComplex(x) ) {

      while ( p.re !== 0n || p.im !== 0n ) {
        for ( var j = 0 ; j < 8*i ; j++ ) p = mul( p, x );
        s = add( s, p );
        i++;
      }

    } else {

      while ( p !== 0n ) {
        for ( var j = 0 ; j < 8*i ; j++ ) p = mul( p, x );
        s = s + p;
        i++;
      }

    }

    return s;

  }

  function arbitraryTheta3( x ) {

    var p = arbitrary(2);
    var s = arbitrary(1);
    var i = 1;

    if ( isComplex(x) ) {

      while ( p.re !== 0n || p.im !== 0n ) {
        for ( var j = 0 ; j < 4*(2*i-1) ; j++ ) p = mul( p, x );
        s = add( s, p );
        i++;
      }

    } else {

      while ( p !== 0n ) {
        for ( var j = 0 ; j < 4*(2*i-1) ; j++ ) p = mul( p, x );
        s = s + p;
        i++;
      }

    }

    return s;

  }

  if ( isArbitrary(x) ) {

    var arb1 = arbitrary(1);

    if ( !isComplex(x) ) {

      if ( x < 0n ) return { re: ln( -x ), im: getConstant( 'pi' ) };

      if ( x === arb1 ) return 0n;

      if ( x < arb1 ) return -ln( div( arb1, x ) );

    }

    if ( abs(x) < arb1 ) return mul( -arb1, ln( div( arb1, x ) ) );

    x = div( arb1, x );

    var t2 = arbitraryTheta2(x);
    var t3 = arbitraryTheta3(x);
    var Pi = getConstant( 'pi' );

    var result = div( Pi, mul( arbitrary(4), arbitraryAGM( mul(t2,t2), mul(t3,t3) ) ) );

    // adjust imaginary part
    if ( x.re < 0n ) {
      if ( result.im > 0n ) result.im -= Pi;
      else result.im += Pi;
    }

    return result;

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

