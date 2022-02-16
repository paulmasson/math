
function exp( x ) {

  if ( isArbitrary(x) ) {

    if ( isComplex(x) ) {

      var expXre = exp(x.re);

      return { re: mul( expXre, cos(x.im) ),
               im: mul( expXre, sin(x.im) ) };

    }

    var m = Math.trunc( arbitrary( div( x, ln10 ) ) );

    if ( m > 0 ) {
      if ( defaultDecimals + m > constants.decimals )
        console.log( 'Not enough decimals in constants for arbitrary exponential' );
      else {
        setPrecisionScale( defaultDecimals + m );
        x *= BigInt( 10**m ); // pad to match new precision
      }
    }

    x -= BigInt(m) * ln10;

    // direct sum faster than function inversion
    var s = arb1;
    var p = arb1;
    var i = arb1;

    while ( p !== 0n ) {
      p = div( mul( p, x ), i );
      s += p;
      i += arb1;
    }

    if ( m > 0 )
      if ( defaultDecimals + m > constants.decimals )
        s *= BigInt( 10**m );
      else
        setPrecisionScale( defaultDecimals );
    else s /= BigInt( 10**-m ); // value approaches zero for fixed decimals

    // could also return as mantissa/exponent
    return s;

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

    var t, u;

    if ( isComplex(x) ) {

      var maxIter = 20, i = 0;

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

    var p = mul( arb2, x );
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

    var p = arb2;
    var s = arb1;
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

    if ( abs(x) < arb1 ) return neg( ln( div( arb1, x ) ) );

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

  if ( ( k === 0 || k === -1 ) && abs( add( x, expMinusOne ) ) < tolerance*tolerance )
    if ( isComplex(x) ) return complex(-1);
    else return -1;

  // handle real cases separately, complex otherwise

  if ( !isComplex(x) ) {
    if ( k === 0 && x > -expMinusOne )
      return findRoot( w => w * Math.exp(w) - x, [-1,1000], { tolerance: tolerance } );
    if ( k === -1 && x > -expMinusOne && x < 0 )
      return findRoot( w => w * Math.exp(w) - x, [-1000,-1], { tolerance: tolerance } );
    x = complex(x);
  }

  // inversion by complex root finding with custom Newton's method
  var maxIter = 100;

  if ( k === 0 && abs(x) < 1.25 ) {
    // unstable along branch cut, add small complex part to avoid error
    if ( x.im === 0 &&  x.re < -expMinusOne ) x.im = tolerance;
    // based on test page: unstable region jumps between sheets
    var w = x.re < .5 ? complex( 0, .5*Math.sign(x.im) ) : complex(0);
  } else {
    var L = add( log(x), complex(0,2*pi*k) );
    var w = sub( L, log(L) );
  }

  for ( var i = 0; i < maxIter ; i++ ) {
     var delta = div( sub( mul(w,exp(w)), x ), mul( exp(w), add(w,1) ) );
     w = sub( w, delta );
     if ( abs(delta) < tolerance ) return w;
  }

  throw Error( 'No Lambert W root found at ' + JSON.stringify(x) );

}

function inverseLambertW( x ) { return mul( x, exp(x) ); }

