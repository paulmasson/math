
function complex( x, y ) {

  var y = y || 0;
  return { re: x, im: y };

}

var C = complex;

function isComplex( x ) { return typeof x === 'object' && 're' in x; }


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

    if ( x.re === 0 && x.im === 0 ) return 0;

    if ( Math.abs(x.re) < Math.abs(x.im) )

      return Math.abs(x.im) * Math.sqrt( 1 + ( x.re / x.im )**2 );

    else

      return Math.abs(x.re) * Math.sqrt( 1 + ( x.im / x.re )**2 );

  }

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

    return { re: x.re * y.re - x.im * y.im,
             im: x.im * y.re + x.re * y.im };

  }

  return x * y;

}

function neg( x ) { return mul( -1, x ); }

function div( x, y ) {

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x);
    if ( !isComplex(y) ) y = complex(y);

    if ( y.re === 0 && y.im === 0 ) throw Error( 'Division by zero' );

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

  if ( y === 0 ) throw Error( 'Division by zero' );

  return x / y;

}

function inv( x ) { return div( 1, x ); }

function pow( x, y ) {

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x);
    if ( !isComplex(y) ) y = complex(y);

    if ( x.re === 0 && x.im === 0 && y.re > 0 )
      return complex(0);
    if ( x.re === 0 && x.im === 0 && y.re === 0 && y.im === 0 )
      return complex(1);
    if ( x.re === 0 && x.im === 0 && y.re < 0 )
      throw Error( 'Power singularity' );

    var r = Math.sqrt( x.re * x.re + x.im * x.im );
    var phi = Math.atan2( x.im, x.re );

    var R = r**y.re * Math.exp( -phi * y.im );
    var Phi = phi * y.re + y.im * Math.log(r);

    return { re: R * Math.cos(Phi), im: R * Math.sin(Phi) };

  }

  if ( x === 0 && y < 0 ) throw Error( 'Power singularity' );

  if ( x < 0 && !Number.isInteger(y) ) return pow( complex(x), y );

  return x**y;

}

function root( x, y ) { return pow( x, div( 1, y ) ); }

function sqrt( x ) {

  if ( isComplex(x) ) {

    var R = ( x.re * x.re + x.im * x.im )**(1/4);
    var Phi = Math.atan2( x.im, x.re ) / 2;

    return { re: R * Math.cos(Phi), im: R * Math.sin(Phi) };

  }

  if ( x < 0 ) return sqrt( complex(x) );

  return Math.sqrt(x);

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

