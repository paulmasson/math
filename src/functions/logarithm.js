
function exp( x ) {

  if ( isArbitrary(x) ) {

    if ( isComplex(x) )

      return { re: mul( exp(x.re), cos(x.im) ),
               im: mul( exp(x.re), sin(x.im) ) };

    var m = Math.trunc( arbitrary( div( x, ln10 ) ) );
    x = x - mul( arbitrary(m), ln10 );

    // direct sum faster than function inversion
    var s = arb1;
    var p = arb1;
    var i = arb1;

    while ( p !== 0n ) {
      p = div( mul( p, x ), i );
      s += p;
      i += arb1;
    }

    // could also return as mantissa/exponent
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

    if ( !isComplex(x) ) {

      if ( x < 0n ) return { re: ln( -x ), im: getConstant( 'pi' ) };

      if ( x === arb1 ) return 0n;

      if ( x < arb1 ) return -ln( div( arb1, x ) );

    }

    if ( abs(x) < arb1 ) return mul( -arb1, ln( div( arb1, x ) ) );

    x = div( arb1, x );

    var t2 = arbitraryTheta2(x);
    var t3 = arbitraryTheta3(x);

    var result = div( onePi, mul( arbitrary(4), arbitraryAGM( mul(t2,t2), mul(t3,t3) ) ) );

    // adjust imaginary part
    if ( x.re < 0n ) {
      if ( result.im > 0n ) result.im -= onePi;
      else result.im += onePi;
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

  // restrict to real integers for convenience
  if ( !Number.isInteger(k) ) throw Error( 'Unsupported Lambert W index' );

  if ( isZero(x) )
    if ( k === 0 ) return x;
    else throw Error( 'Lambert W pole' );

  var expMinusOne = Math.exp(-1);

  if ( abs( add( x, expMinusOne ) ) < tolerance*tolerance && ( k === 0 || k === -1 ) )
    if ( isComplex(x) ) return complex(-1);
    else return -1;

  // inversion by complex root finding

  if ( abs(x) <= 1 && k === 0 ) var start = complex(0);
  else {
    var L = add( log(x), complex(0,2*pi*k) );
    var start = add( L, neg(log(L)) );
  }

  var result = findRoot( w => sub( mul(w,exp(w)), x ), start, { tolerance: tolerance } );

  if ( !isComplex(x) ) {
    if ( x > -expMinusOne && k === 0 ) return result.re;
    if ( x > -expMinusOne && x < 0 && k === -1 ) return result.re;
  }

  return result;

}

