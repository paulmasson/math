
function complex( x, y=0 ) {

  if ( y === 0 && isArbitrary(x) ) y = 0n;

  return { re: x, im: y };

}

var C = complex;

function isComplex( x ) { return typeof x === 'object' && 're' in x; }


var decimals, precisionScale, arb1, arb2, onePi, twoPi, halfPi, ln10;

function setPrecisionScale( n ) {

  decimals = n;
  precisionScale = 10n**BigInt(decimals);

  // set some commonly used constants
  arb1 = arbitrary(1);
  arb2 = 2n * arb1;
  onePi = getConstant( 'pi' );
  twoPi = 2n * onePi;
  halfPi = onePi / 2n;
  ln10 = getConstant( 'ln10' );

}

setPrecisionScale( defaultDecimals );

function arbitrary( x ) {

  if ( isComplex(x) ) return { re: arbitrary( x.re ), im: arbitrary( x.im ) };

  if ( isArbitrary(x) ) return Number(x) / 10**decimals;

  // BigInt from exponential form includes wrong digits
  // manual construction from string more accurate

  var parts = x.toExponential().split( 'e' );
  var mantissa = parts[0].replace( '.', '' );
  var digits = mantissa.length - ( mantissa[0] === '-' ? 2 : 1 )
  var padding = +parts[1] + decimals - digits;

  if ( padding < 0 ) return BigInt( Math.round( x * 10**decimals ) );

  return BigInt( mantissa + '0'.repeat(padding) );

}

var A = arbitrary;

function isArbitrary( x ) { return typeof x === 'bigint' || typeof x.re === 'bigint'; }


function isZero( x ) {

  if ( isComplex(x) ) return x.re === 0 && x.im === 0;
  return x === 0;

}

function isInteger( x ) {

  if ( isComplex(x) ) return Number.isInteger(x.re) && x.im === 0;
  return Number.isInteger(x);

}

function isPositiveInteger( x ) {

  if ( isComplex(x) ) return Number.isInteger(x.re) && x.re > 0 && x.im === 0;
  return Number.isInteger(x) && x > 0;

}

function isPositiveIntegerOrZero( x ) {

  if ( isComplex(x) ) return Number.isInteger(x.re) && x.re >= 0 && x.im === 0;
  return Number.isInteger(x) && x >= 0;

}

function isNegativeInteger( x ) {

  if ( isComplex(x) ) return Number.isInteger(x.re) && x.re < 0 && x.im === 0;
  return Number.isInteger(x) && x < 0;

}

function isNegativeIntegerOrZero( x ) {

  if ( isComplex(x) ) return Number.isInteger(x.re) && x.re <= 0 && x.im === 0;
  return Number.isInteger(x) && x <= 0;

}

function isEqualTo( x, y ) {

  if ( isComplex(x) || isComplex(y) ) {
    if ( !isComplex(x) ) x = complex(x);
    if ( !isComplex(y) ) y = complex(y);
    return x.re === y.re && x.im === y.im;
  }

  return x === y;

}


function re( x ) {

  if ( isComplex(x) ) return x.re;
  return x;

}

var real = re;

function im( x ) {

  if ( isComplex(x) ) return x.im;
  return 0;

}

var imag = im;

function abs( x ) {

  if ( isComplex(x) ) {

    if ( x.re === 0 && x.im === 0 ) return complex(0);

    if ( isArbitrary(x) ) return sqrt( mul(x.re,x.re) + mul(x.im,x.im) );

    if ( Math.abs(x.re) < Math.abs(x.im) )

      return Math.abs(x.im) * Math.sqrt( 1 + ( x.re / x.im )**2 );

    else

      return Math.abs(x.re) * Math.sqrt( 1 + ( x.im / x.re )**2 );

  }

  if ( isArbitrary(x) )
    if ( x < 0n ) return -x;
    else return x;

  return Math.abs(x);

}

function arg( x ) {

  if ( isComplex(x) ) return Math.atan2( x.im, x.re );

  return Math.atan2( 0, x );

}


// JavaScript does not support operator overloading

function add( x, y ) {

  if ( arguments.length > 2 ) {

    var z = add( x, y );
    for ( var i = 2 ; i < arguments.length ; i++ ) z = add( z, arguments[i] );
    return z; 

  }

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x);
    if ( !isComplex(y) ) y = complex(y);

    return { re: x.re + y.re, im: x.im + y.im };

  }

  return x + y;

}

function sub( x, y ) {

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x);
    if ( !isComplex(y) ) y = complex(y);

    return { re: x.re - y.re, im: x.im - y.im };

  }

  return x - y;

}

function mul( x, y ) {

  if ( arguments.length > 2 ) {

    var z = mul( x, y );
    for ( var i = 2 ; i < arguments.length ; i++ ) z = mul( z, arguments[i] );
    return z; 

  }

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x);
    if ( !isComplex(y) ) y = complex(y);

    if ( isArbitrary(x) )

      return { re: ( x.re * y.re - x.im * y.im ) / precisionScale,
               im: ( x.im * y.re + x.re * y.im ) / precisionScale };

    return { re: x.re * y.re - x.im * y.im,
             im: x.im * y.re + x.re * y.im };

  }

  if ( isArbitrary(x) ) return x * y / precisionScale;

  return x * y;

}

function neg( x ) {

  if ( isComplex(x) ) return { re: -x.re, im: -x.im };

  return -x;

}

function div( x, y ) {

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x);
    if ( !isComplex(y) ) y = complex(y);

    if ( y.re === 0 && y.im === 0 || y.re === 0n && y.im === 0n ) // operator precedence
      throw Error( 'Division by zero' );

    if ( isArbitrary(x) ) {

      var N = { re: x.re * y.re + x.im * y.im,
                im: x.im * y.re - x.re * y.im };
      var D = y.re * y.re + y.im * y.im;

      return { re: precisionScale * N.re / D,
               im: precisionScale * N.im / D };

    }

    if ( Math.abs(y.re) < Math.abs(y.im) ) {

      var f = y.re / y.im;
      return { re: ( x.re * f + x.im ) / ( y.re * f + y.im ),
               im: ( x.im * f - x.re ) / ( y.re * f + y.im ) };

    } else {

      var f = y.im / y.re;
      return { re: ( x.re + x.im * f ) / ( y.re + y.im * f ),
               im: ( x.im - x.re * f ) / ( y.re + y.im * f ) };

    }

  }

  if ( y === 0 || y === 0n ) throw Error( 'Division by zero' );

  if ( isArbitrary(x) ) return precisionScale * x / y;

  return x / y;

}

function inv( x ) { return div( 1, x ); }

function pow( x, y ) {

  if ( isArbitrary(x) || isArbitrary(y) ) {

    if ( !isArbitrary(x) ) x = arbitrary(x);
    if ( !isArbitrary(y) ) y = arbitrary(y);

    return exp( mul( y, ln(x) ) );

  }

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x);
    if ( !isComplex(y) ) y = complex(y);

    if ( x.re === 0 && x.im === 0 && y.re > 0 )
      return complex(0);
    if ( x.re === 0 && x.im === 0 && y.re === 0 && y.im === 0 )
      return complex(1);
    if ( x.re === 0 && x.im === 0 && y.re < 0 )
      throw Error( 'Power singularity' );

    return exp( mul( y, log(x) ) );

  }

  if ( x === 0 && y < 0 ) throw Error( 'Power singularity' );

  if ( x < 0 && !Number.isInteger(y) ) return pow( complex(x), y );

  return x**y;

}

function root( x, y ) { return pow( x, div( 1, y ) ); }

function surd( x, n ) {

  if ( isComplex(x) || isComplex(n) ) throw Error( 'Surd requires real inputs' );

  if ( !isInteger(n) ) throw Error( 'Second parameter of surd must be integer' );

  if ( n & 1 ) {
    var sign = Math.sign(x); // zero at origin anyway
    return sign * root( sign*x, n );
  }

  if ( x < 0 ) throw Error( 'First parameter of surd must be positive for even integers' );

  return root( x, n );

}

function sqrt( x ) {

  if ( isComplex(x) ) {

    if ( isArbitrary(x) ) {

      if ( x.im === 0n )
        if ( x.re < 0n ) return { re: 0n, im: sqrt(-x.re) };
        else return { re: sqrt(x.re), im: 0n };

      // need evaluation independent of natural logarithm

      var c = abs(x);
      var sign = x.im < 0n ? -1n : 1n;

      return { re: sqrt( div( c + x.re, arb2 ) ), im: sign * sqrt( div( c - x.re, arb2 ) ) }

    }

    if ( x.im === 0 )
      if ( x.re < 0 ) return { re: 0, im: Math.sqrt(-x.re) };
      else return { re: Math.sqrt(x.re), im: 0 };

    // expression above suffers from rounding errors when not arbitrary,
    //   especially affecting elliptic integral of third kind

    return exp( mul( .5, log(x) ) );

  }

  if ( isArbitrary(x) ) {

    if ( x === 0n ) return 0n;

    if ( x < 0n ) throw Error( 'Cannot evaluate real square root of ' + x );

    // Brent, Modern Computer Arithmetic, SqrtInt algorithm

    var u = x, s, t;

    while ( u !== s ) {
      s = u;
      t = s + div( x, s );
      u = div( t, arb2 );
    }

    return s;

  }

  if ( x < 0 ) return { re: 0, im: Math.sqrt(-x) };

  return Math.sqrt(x);

}


function complexAverage( f, x, offset=1e-5 ) {

  return div( add( f(add(x,offset)), f(sub(x,offset)) ), 2 );

}


function complexFromString( s, returnAsString=false ) {

  var lead = '', real, imag;

  if ( s[0] === '+' || s[0] === '-' ) {
    lead = s[0];
    s = s.slice(1);
  }

  if ( s.includes('+') || s.includes('-') ) {
    if ( s.includes('+') ) {
      real = lead + s.slice( 0, s.indexOf('+') );
      imag = s.slice( s.indexOf('+') + 1, s.length - 1 );
    } else {
      real = lead + s.slice( 0, s.indexOf('-') );
      imag = s.slice( s.indexOf('-'), s.length - 1 );
    }
  } else {
    if ( s.includes('i') ) {
      real = '0';
      imag = lead + s.slice( 0, s.length - 1 );
    } else {
      real = lead + s;
      imag = '0';
    }
  }

  if ( imag === '' || imag === '-' ) imag += '1';

  if ( returnAsString ) return `{ re: ${real}, im: ${imag} }`;

  return { re: +real, im: +imag };

}

