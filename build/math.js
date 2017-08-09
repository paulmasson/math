
var pi = Math.PI;

var e = Math.exp(1);


function complex( x, y ) {

  var y = y || 0;
  return { re: x, im: y };

}

var C = complex;

// arrays of length greater than two are true for this test
// need a fast way to distinguish between array and dictionary

var isComplex = isNaN;


// JavaScript does not yet support operator overloading

function add( x, y ) {

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x,0);
    if ( !isComplex(y) ) y = complex(y,0);

    return { re: x.re + y.re, im: x.im + y.im };

  } else return x + y;

}

function sub( x, y ) {

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x,0);
    if ( !isComplex(y) ) y = complex(y,0);

    return { re: x.re - y.re, im: x.im - y.im };

  } else return x - y;

}

function mul( x, y ) {

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x,0);
    if ( !isComplex(y) ) y = complex(y,0);

    return { re: x.re * y.re - x.im * y.im,
             im: x.im * y.re + x.re * y.im };

  } else return x * y;

}

function div( x, y ) {

  // need to handle 0/0...

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x,0);
    if ( !isComplex(y) ) y = complex(y,0);

    var ySq = y.re * y.re + y.im * y.im;
    return { re: ( x.re * y.re + x.im * y.im ) / ySq,
             im: ( x.im * y.re - x.re * y.im ) / ySq };

  } else return x / y;

}

function pow( x, y ) {

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x,0);
    if ( !isComplex(y) ) y = complex(y,0);

    var r = Math.sqrt( x.re * x.re + x.im * x.im );
    var phi = Math.atan2( x.im, x.re );

    var R = r**y.re * Math.exp( -phi * y.im );
    var Phi = phi * y.re + y.im * Math.log(r); // NaN at origin...

    return { re: R * Math.cos(Phi), im: R * Math.sin(Phi) };

  } else return x**y;


}

function root( x, y ) { return pow( x, div( 1, y ) ); }

function sqrt( x ) {

  if ( isComplex(x) ) {

    var R = ( x.re * x.re + x.im * x.im )**(1/4);
    var Phi = Math.atan2( x.im, x.re ) / 2;

    return { re: R * Math.cos(Phi), im: R * Math.sin(Phi) };

  } else return Math.sqrt(x);

}



function besselJ( n, x ) {

  return (x/2)**n * hypergeometric0F1( n+1, -.25*x**2 ) / gamma(n+1);

}

function besselY( n, x ) {

  return ( besselJ(n,x) * cos(n*pi) - besselJ(-n,x) ) / sin(n*pi);


}

function besselI( n, x ) {

  return (x/2)**n * hypergeometric0F1( n+1, .25*x**2 ) / gamma(n+1);

}

function besselK( n, x ) {

  return pi/2 * ( besselI(-n,x) - besselI(n,x) ) / sin(n*pi);


}


function airyAi( x ) {

  if ( x > 0 ) return 1/pi * sqrt(x/3) * besselK(1/3, 2/3*x**(3/2));

  if ( x === 0 ) return 1 / 3**(2/3) / gamma(2/3);

  if ( x < 0 ) return sqrt(-x) / 2 * (  besselJ(1/3, 2/3*(-x)**(3/2))
                                       - besselY(1/3, 2/3*(-x)**(3/2)) / sqrt(3) );

}

function airyBi( x ) {

  if ( x > 0 ) return sqrt(x/3) * ( besselI(1/3, 2/3*x**(3/2))
                                    + besselI(-1/3, 2/3*x**(3/2)) );

  if ( x === 0 ) return 1 / 3**(1/6) / gamma(2/3);

  if ( x < 0 ) return -sqrt(-x) / 2 * ( besselJ(1/3, 2/3*(-x)**(3/2)) / sqrt(3)
                                        + besselY(1/3, 2/3*(-x)**(3/2)) );

}

function jacobiTheta( n, x, q ) {

  switch( n ) {

    case 1:

      var s = 0;
      for ( var i = 0 ; i < 100 ; i++ ) s += (-1)**i * q**(i*i+i) * sin( (2*i+1) * x );
      return 2 * q**(1/4) * s;

    case 2:

      var s = 0;
      for ( var i = 0 ; i < 100 ; i++ ) s += q**(i*i+i) * cos( (2*i+1) * x );
      return 2 * q**(1/4) * s;

    case 3:

      var s = 0;
      for ( var i = 1 ; i < 100 ; i++ ) s += q**(i*i) * cos( 2*i * x );
      return 1 + 2 * s;

    case 4:

      var s = 0;
      for ( var i = 1 ; i < 100 ; i++ ) s += (-q)**(i*i) * cos( 2*i * x );
      return 1 + 2 * s;

    default:

      throw( 'Undefined Jacobi theta index' );

  }

}


function sn( x, m ) {

  var q = exp( -pi * ellipticK(1-m) / ellipticK(m) );

  var t = x / jacobiTheta(3,0,q)**2;

  return jacobiTheta(3,0,q) / jacobiTheta(2,0,q)
         * jacobiTheta(1,t,q) / jacobiTheta(4,t,q);

}

function cn( x, m ) {

  var q = exp( -pi * ellipticK(1-m) / ellipticK(m) );

  var t = x / jacobiTheta(3,0,q)**2;

  return jacobiTheta(4,0,q) / jacobiTheta(2,0,q)
         * jacobiTheta(2,t,q) / jacobiTheta(4,t,q);

}
function dn( x, m ) {

  var q = exp( -pi * ellipticK(1-m) / ellipticK(m) );

  var t = x / jacobiTheta(3,0,q)**2;

  return jacobiTheta(4,0,q) / jacobiTheta(3,0,q)
         * jacobiTheta(3,t,q) / jacobiTheta(4,t,q);

}


// Carlson symmetric integrals

function carlsonRC( x, y, z, tolerance=1e-10 ) {

  return 1;

}

function carlsonRD( x, y, z, tolerance=1e-10 ) {

  return 1;

}

function carlsonRF( x, y, z, tolerance=1e-10 ) {

  if ( y === z ) return carlsonRC( x, y );
  if ( x === z ) return carlsonRC( y, x );
  if ( x === y ) return carlsonRC( z, x );

  var xm = x;
  var ym = y;
  var zm = z;

  var A0 = (x+y+z)/3;
  var Am = A0;
  var Q = Math.pow( 3*tolerance, -1/6 )
          * Math.max( Math.abs(A0-x), Math.abs(A0-y), Math.abs(A0-z) );
  var g = 0.25;
  var pow4 = 1;
  var m = 0;
  while ( true ) {
    var xs = Math.sqrt(xm);
    var ys = Math.sqrt(ym);
    var zs = Math.sqrt(zm);
    var lm = xs*ys + xs*zs + ys*zs;
    var Am1 = ( Am + lm ) * g;
    xm = ( xm + lm ) * g;
    ym = ( ym + lm ) * g;
    zm = ( zm + lm ) * g;
    if ( pow4 * Q < Math.abs(Am) ) break;
    Am = Am1;
    m += 1;
    pow4 *= g;
  }
  var t = pow4/Am;
  var X = (A0-x)*t;
  var Y = (A0-y)*t;
  var Z = -X-Y;
  var E2 = X*Y-Z**2;
  var E3 = X*Y*Z;

  return Math.pow(Am,-0.5) * (9240-924*E2+385*E2**2+660*E3-630*E2*E3)/9240;

}

function carlsonRG( x, y, z, tolerance=1e-10 ) {

  return 1;

}

function carlsonRJ( x, y, z, tolerance=1e-10 ) {

  return 1;

}


// elliptic integrals

function ellipticF( x, m ) {

  if ( arguments.length === 1 ) {
    m = x;
    x = pi / 2;
  }

  var period = 0;
  if ( Math.abs(x) > pi / 2 ) {
    var p = Math.round( x / pi );
    x = x - p * pi;
    period = 2 * p * ellipticK( m );
  }

  return sin(x) * carlsonRF( cos(x)**2, 1 - m * sin(x)**2, 1 ) + period;

}

function ellipticK( m ) {

  return ellipticF( m );

}

function ellipticE( x, m ) {

  if ( arguments.length === 1 ) {
    m = x;
    x = pi / 2;
  }

  var period = 0;
  if ( Math.abs(x) > pi / 2 ) {
    var p = Math.round( x / pi );
    x = x - p * pi;
    period = 2 * p * ellipticE( m );
  }

  return sin(x) * carlsonRF( cos(x)**2, 1 - m * sin(x)**2, 1 )
         - m / 3 * sin(x)**3 * carlsonRD( cos(x)**2, 1 - m * sin(x)**2, 1 )
         + period;

}

function ellipticPi( n, x, m ) {

  if ( arguments.length === 2 ) {
    m = x;
    x = pi / 2;
  }

  var period = 0;
  if ( Math.abs(x) > pi / 2 ) {
    var p = Math.round( x / pi );
    x = x - p * pi;
    period = 2 * p * ellipticPi( n, m );
  }

  return sin(x) * carlsonRF( cos(x)**2, 1 - m * sin(x)**2, 1 )
         + n / 3 * sin(x)**3
           * carlsonRJ( cos(x)**2, 1 - m * sin(x)**2, 1, 1 - n * sin(x)**2 )
         + period;

}


function factorial( n ) {

  if ( Number.isInteger(n) && n >= 0 ) {

    var result = 1;
    for ( var i = 2 ; i <= n ; i++ ) result *= i;
    return result;

  } else return gamma( n+1 );

}

function binomial( n, m ) {

  return factorial(n) / factorial(n-m) / factorial(m);

}


// log of gamma less likely to overflow than gamma
// Lanczos approximation as evaluated by Paul Godfrey

function logGamma( x ) {

  var c = [ 57.1562356658629235, -59.5979603554754912, 14.1360979747417471,
            -0.491913816097620199, .339946499848118887e-4, .465236289270485756e-4,
            -.983744753048795646e-4, .158088703224912494e-3, -.210264441724104883e-3,
            .217439618115212643e-3, -.164318106536763890e-3, .844182239838527433e-4,
            -.261908384015814087e-4, .368991826595316234e-5 ];

  if ( isComplex(x) ) {

    if ( Number.isInteger(x.re) && x.re <= 0 && x.im === 0 )
      throw( 'Gamma function pole' );

    // reflection formula
    if ( x.re < 0 ) {
      var y = mul( -1, x );
      return sub( log( div( -pi, mul( y, sin( mul(pi,y) ) ) ) ), logGamma(y) );
    }

    var t = add( x, 5.24218750000000000 );
    t = sub( mul( add( x, 0.5 ), log(t)), t );
    var s = 0.999999999999997092;
    for ( var j = 0 ; j < 14 ; j++ ) s = add( s, div( c[j], add( x, j+1 ) ) );
    return add( t, log( mul( 2.5066282746310005, div( s, x ) ) ) );

  } else {

    if ( Number.isInteger(x) && x <= 0 )
      throw( 'Gamma function pole' ); 

    var t = x + 5.24218750000000000;
    t = ( x + 0.5 ) * log(t) - t;
    var s = 0.999999999999997092;
    for ( var j = 0 ; j < 14 ; j++ ) s += c[j] / (x+j+1);
    return t + log( 2.5066282746310005 * s / x );

  }

}

function gamma( x ) {

  // logGamma complex on negative axis
  if ( !isComplex(x) && x < 0 )
    return exp( logGamma( complex(x) ) ).re;
  else return exp( logGamma(x) );

}

function beta( x, y ) {

  return div( mul( gamma(x), gamma(y) ), gamma( add(x,y) ) ); 

}


function hypergeometric0F1( a, x ) {

  if ( Number.isInteger(a) && a <= 0 ) throw( 'Hypergeometric function pole' );

  var s = 1;
  var p = 1;

  for ( var i = 1 ; i < 100 ; i++ ) {
    p *= x / a / i;
    s += p;
    a++;
  }

  return s;

}


function hypergeometric1F1( a, b, x ) {

  if ( Number.isInteger(b) && b <= 0 ) throw( 'Hypergeometric function pole' );

  var s = 1;
  var p = 1;

  for ( var i = 1 ; i < 100 ; i++ ) {
    p *= x * a / b / i;
    s += p;
    a++;
    b++;
  }

  return s;

}


function hypergeometric2F1( a, b, c, x ) {

  if ( Number.isInteger(c) && c <= 0 ) throw( 'Hypergeometric function pole' );

  if ( x === 1 ) return gamma(c) * gamma(c-a-b) / gamma(c-a) / gamma(c-b);

  var s = 1;
  var p = 1;

  for ( var i = 1 ; i < 1000 ; i++ ) {
    p *= x * a * b / c / i;
    s += p;
    a++;
    b++;
    c++;
  }

  return s;

}


function exp( x ) {

  if ( isComplex(x) )

    return { re: Math.exp(x.re) * Math.cos(x.im),
             im: Math.exp(x.re) * Math.sin(x.im) };

  else return Math.exp(x);

}


function log( x, base ) {

  if ( isComplex(x) ) {

    var r = sqrt( x.re * x.re + x.im * x.im );
    var phi = Math.atan2( x.im, x.re );

    return { re: log(r,base), im: log(e,base) * phi };

  } else if ( x < 0 ) return log( complex(x), base );

  else if ( base === undefined ) return Math.log(x);

  else return Math.log(x) / Math.log(base);

}

var ln = log;


// rounding functions

function roundTo( x, n ) {

  if ( x === 0 ) return x;

  if ( Array.isArray(x) ) {
    var v = vector( x.length );
    for ( var i = 0 ; i < x.length ; i++ ) v[i] = roundTo( x[i], n );
    return v;
  }

  var exponent = Math.floor( Math.log10( Math.abs(x) ) );
  n = n - exponent - 1;
  return Math.round( 10**n * x ) / 10**n;

}

function ceilTo( x, n ) {

  if ( x === 0 ) return x;

  if ( Array.isArray(x) ) {
    var v = vector( x.length );
    for ( var i = 0 ; i < x.length ; i++ ) v[i] = ceilTo( x[i], n );
    return v;
  }

  var exponent = Math.floor( Math.log10( Math.abs(x) ) );
  n = n - exponent - 1;
  return Math.ceil( 10**n * x ) / 10**n;

}

function floorTo( x, n ) {

  if ( x === 0 ) return x;

  if ( Array.isArray(x) ) {
    var v = vector( x.length );
    for ( var i = 0 ; i < x.length ; i++ ) v[i] = floorTo( x[i], n );
    return v;
  }

  var exponent = Math.floor( Math.log10( Math.abs(x) ) );
  n = n - exponent - 1;
  return Math.floor( 10**n * x ) / 10**n;

}

function chop( x, tolerance=1e-10 ) {

  if ( Array.isArray(x) ) {
    var v = vector( x.length );
    for ( var i = 0 ; i < x.length ; i++ ) v[i] = chop( x[i] );
    return v;
  }

  if ( isComplex(x) ) return { re: chop(x.re), im: chop(x.im) };

  if ( Math.abs(x) < tolerance ) x = 0;
  return x;

}


function kronecker( i, j ) {

  return i === j ? 1 : 0;

}


function hermite( n, x ) {

  function coefficients( n ) {

    var minus2 = [ 1 ];
    var minus1 = [ 2, 0 ];
    var t, current;

    if ( n === 0 ) return minus2;
    if ( n === 1 ) return minus1;

    for ( var i = 2 ; i <= n ; i++ ) {
      current = [];
      t = minus1.slice();
      t.push( 0 );
      minus2.unshift( 0, 0 );
      for ( var k = 0 ; k < t.length ; k++ )
        current.push( 2*t[k] - 2*(i-1)*minus2[k] );
      minus2 = minus1;
      minus1 = current;
    }

    return current;

  }

  if ( Number.isInteger(n) && n >= 0 )
    return polynomial( x, coefficients(n) );

  else {

    var s = hypergeometric1F1( -n/2, 1/2, x**2 ) / gamma( (1-n)/2 )
            - 2 * x * hypergeometric1F1( (1-n)/2, 3/2, x**2 ) / gamma( -n/2 );
    return 2**n * sqrt(pi) * s;

  }

}


// complex circular functions

function sin( x ) {

  if ( isComplex(x) )

    return { re: Math.sin(x.re) * Math.cosh(x.im),
             im: Math.cos(x.re) * Math.sinh(x.im) };

  else return Math.sin(x);

}

function cos( x ) {

  if ( isComplex(x) )

    return { re: Math.cos(x.re) * Math.cosh(x.im),
             im: -Math.sin(x.re) * Math.sinh(x.im) };

  else return Math.cos(x);

}

function tan( x ) {

  if ( isComplex(x) ) return div( sin(x), cos(x) );

  else return Math.tan(x);

 }

function cot( x ) {

 if ( isComplex(x) ) return div( cos(x), sin(x) );

  else return 1 / Math.tan(x);

}

function sec( x ) {

  if ( isComplex(x) ) return div( 1, cos(x) );

  else return 1 / Math.cos(x);

}

function csc( x ) {

  if ( isComplex(x) ) return div( 1, sin(x) );

  else return 1 / Math.sin(x);

}


// inverse circular functions

function arcsin( x ) {

  if ( isComplex(x) ) {

    var s = sqrt( sub( 1, mul( x, x ) ) );
    s = add( mul( complex(0,1), x ), s ); 
    return mul( complex(0,-1), log( s ) );

  } else if ( Math.abs(x) <= 1 ) return Math.asin(x);

  else return arcsin( complex(x) );

}

function arccos( x ) {

  if ( isComplex(x) ) {

    return sub( pi/2, arcsin(x) );

  } else if ( Math.abs(x) <= 1 ) return Math.acos(x);

  else return arccos( complex(x) );

}

function arctan( x ) {

  if ( isComplex(x) ) {

    var s = sub( log( sub( 1, mul( complex(0,1), x ) ) ),
                 log( add( 1, mul( complex(0,1), x ) ) ) );
    return mul( complex(0,1/2), s );

  } else return Math.atan(x);

}

function arccot( x ) {

  if ( isComplex(x) ) return arctan( div( 1, x ) );

  else return Math.atan( 1/x );

}

function arcsec( x ) {

  if ( isComplex(x) ) return arccos( div( 1, x ) );

  else if ( Math.abs(x) >= 1 ) return Math.acos( 1/x );

  else return arcsec( complex(x) );

}

function arccsc( x ) {

  if ( isComplex(x) ) return arcsin( div( 1, x ) );

  else if ( Math.abs(x) >= 1 ) return Math.asin( 1/x );

  else return arccsc( complex(x) );

}


// complex hyperbolic functions

function sinh( x ) {

  if ( isComplex(x) )

    return { re: Math.sinh(x.re) * Math.cos(x.im),
             im: Math.cosh(x.re) * Math.sin(x.im) };

  else return Math.sinh(x);

}

function cosh( x ) {

  if ( isComplex(x) )

    return { re: Math.cos(x.re) * Math.cosh(x.im),
             im: Math.sin(x.re) * Math.sinh(x.im) };

  else return Math.cosh(x);

}

function tanh( x ) {

  if ( isComplex(x) ) return div( sinh(x), cosh(x) );

  else return Math.tanh(x);

}

function coth( x ) {

  if ( isComplex(x) ) return div( cosh(x), sinh(x) );

  else return 1 / Math.tanh(x);

}

function sech( x ) {

  if ( isComplex(x) ) return div( 1, cosh(x) );

  else return 1 / Math.cosh(x);

}

function csch( x ) {

  if ( isComplex(x) ) return div( 1, sinh(x) );

  else return 1 / Math.sinh(x);

}


// inverse hyperbolic functions

function arcsinh( x ) {

  if ( isComplex(x) ) {

    var s = sqrt( add( mul( x, x ), 1 ) );
    s = add( x, s );
    return log( s );

  } else return Math.asinh(x);

}

function arccosh( x ) {

  if ( isComplex(x) ) {

    var s = mul( sqrt( add( x, 1 ) ), sqrt( sub( x, 1 ) ) );
    s = add( x, s ); 
    return log( s );

  } else if ( x >= 1 ) return Math.acosh(x);

  else return arccosh( complex(x) );

}

function arctanh( x ) {

  if ( isComplex(x) ) {

    var s = sub( log( add( 1, x ) ), log( sub( 1, x ) ) );
    return mul( 1/2, s );

  } else if ( Math.abs(x) <= 1 ) return Math.atanh(x);

  else return arctanh( complex(x) );

}

function arccoth( x ) {

  if ( isComplex(x) ) {

    if ( x.re === 0 && x.im === 0 ) throw( 'Indeterminate value' );

    return arctanh( div( 1, x ) );

  } else if ( Math.abs(x) > 1 ) return Math.atanh( 1/x );

  else return arccoth( complex(x) );

}

function arcsech( x ) {

  if ( isComplex(x) ) {

    if ( x.re === 0 && x.im === 0 ) throw( 'Indeterminate value' );

    // adjust for branch cut along negative axis
    if ( x.im === 0 ) x.im = -Number.MIN_VALUE;

    return arccosh( div( 1, x ) );

  } else if ( x > 0 && x < 1 ) return Math.acosh( 1/x );

  else return arcsech( complex(x) );

}

function arccsch( x ) {

  if ( isComplex(x) ) {

    return arcsinh( div( 1, x ) );

  } else return Math.asinh( 1/x );

}


function zeta( x ) {

  // functional equation
  if ( x < 0 ) return 2**x * pi**(x-1) * sin(pi*x/2) * gamma(1-x) * zeta(1-x); 

  // summation of Dirichlet eta function
  var s = 0;
  for ( var i = 1 ; i < 1e7 ; i++ ) {
    var last = s;
    s -= (-1)**i / i**x;
    if ( Math.abs(s-last) < 1e-8 ) break;
  }

  return s / ( 1 - 2**(1-x) );

}



function ode( f, y, [x0, x1], step=.001, method='euler' ) {

  if ( f(x0,y)[0] === undefined ) {
    g = f;
    f = function(x) { return [ g(x) ]; };
    y = [ y ];
  }

  var points = [ [x0].concat(y) ];

  switch( method ) {

    case 'euler':

      for ( var x = x0+step ; x <= x1 ; x += step ) {

        var k = [];
        for ( var i = 0 ; i < y.length ; i++ ) k.push( f(x,y)[i] * step );

        for ( var i = 0 ; i < y.length ; i++ ) y[i] += k[i];
        points.push( [x].concat(y) );

      }

      return points;

    case 'runge-kutta':

      for ( var x = x0+step ; x <= x1 ; x += step ) {

        var k1 = [], k2 = [], k3 = [], k4 = [];
        var y1 = [], y2 = [], y3 = [];

        for ( var i = 0 ; i < y.length ; i++ ) k1.push( f(x,y)[i] );
        for ( var i = 0 ; i < y.length ; i++ ) y1.push( y[i] + k1[i]*step/2 );
        for ( var i = 0 ; i < y.length ; i++ ) k2.push( f( x+step/2, y1 )[i] );
        for ( var i = 0 ; i < y.length ; i++ ) y2.push( y[i] + k2[i]*step/2 );
        for ( var i = 0 ; i < y.length ; i++ ) k3.push( f( x+step/2, y2 )[i] );
        for ( var i = 0 ; i < y.length ; i++ ) y3.push( y[i] + k3[i]*step );
        for ( var i = 0 ; i < y.length ; i++ ) k4.push( f( x+step, y3 )[i] );

        for ( var i = 0 ; i < y.length ; i++ )
          y[i] += ( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] ) * step / 6;

        points.push( [x].concat(y) );

      }

      return points;

    default:

      throw( 'Unsupported method' );

  }

}





function diff( f, x, n=1, method='ridders' ) {

  // central differences have h**2 error but division
  //   by h**n increases roundoff error
  // step sizes chosen as epsilon**(1/(n+2)) to minimize error

  function difference() {

    var s = 0;
    for ( var i = 0 ; i <= n ; i++ )
      s += (-1)**i * binomial(n,i) * f( x + (n-2*i)*h );

    return s / (2*h)**n

  }

  switch( method ) {

    case 'naive':

      // only accurate for first couple derivatives
      var h = (1e-8)**(1/(n+2));
      return difference();

    case 'ridders':

      var h = (1e-5)**(1/(n+2));
      var error = Number.MAX_VALUE;
      var maxIter = 10;
      var result;

      var d = [];
      for ( var i = 0 ; i < maxIter ; i++ ) d.push( [] );

      // Richardson extrapolation as per C. Ridders
      d[0][0] = difference();

      for ( var i = 1 ; i < maxIter ; i++ ) {

        h /= 2;
        d[0][i] = difference();

        for ( var j = 1 ; j <= i ; j++ ) {

          d[j][i] = ( 4**j * d[j-1][i] - d[j-1][i-1] ) / ( 4**j - 1 );

          var delta = Math.max( Math.abs( d[j][i] - d[j-1][i] ),
                                Math.abs( d[j][i] - d[j-1][i-1] ) );
          if ( delta <= error ) {
            error = delta;
            result = d[j][i];
          }

        }

        if ( Math.abs( d[i][i] - d[i-1][i-1] ) > error ) break;

      }

      return result;

    default:

      throw( 'Unsupported method' );

  }

}

var D = diff;


function integrate( f, a, b, method='adaptive-simpson') {

  var h = ( b - a ) / 2;
  var s = ( f(a) + f(b) ) / 2 + f( (a+b)/2 );

  function nextIter() {

      h /= 2;
      var x = a + h;
      while ( x < b ) {
        // only add new function evaluations
        s += f(x);
        x += 2*h;
      }

  }

  switch( method ) {

    case 'euler-maclaurin':

      // Euler-Maclaurin summation formula

      var tolerance = 1e-10;
      var maxIter = 50;

      var result = h * s;
      var previous = result;

      for ( var i = 0 ; i < maxIter ; i++ ) {

        nextIter();
        result = h * s;
        if ( Math.abs( result - previous ) < tolerance * Math.abs(previous) ) return result;
        previous = result;

      }

      throw( 'Maximum interations reached' );

    case 'romberg':

      var error = Number.MAX_VALUE;
      var maxIter = 30;
      var result = h * s;

      var d = [];
      for ( var i = 0 ; i < maxIter ; i++ ) d.push( [] );

      // Richardson extrapolation of Euler-Maclaurin trapezoids
      d[0][0] = result;

      for ( var i = 1 ; i < maxIter ; i++ ) {

        nextIter();
        d[0][i] = h * s;

        for ( var j = 1 ; j <= i ; j++ ) {

          d[j][i] = ( 4**j * d[j-1][i] - d[j-1][i-1] ) / ( 4**j - 1 );

          var delta = Math.max( Math.abs( d[j][i] - d[j-1][i] ),
                                Math.abs( d[j][i] - d[j-1][i-1] ) );
          if ( delta <= error ) {
            error = delta;
            result = d[j][i];
          }

        }

        if ( Math.abs( d[i][i] - d[i-1][i-1] ) > error ) break;

      }

      return result;

    case 'adaptive-simpson':

      // algorithm by Charles Collins

      var tolerance = 1e-8;
      var maxIter = 50;

      function adaptiveSimpson( a, b, fa, fm, fb, s, tolerance, depth ) {

        var h = b - a;
        var f1 = f( a + h/4 );
        var f2 = f( b - h/4 )

        var s1 = ( fa + 4*f1 + fm ) * h / 12;
        var s2 = ( fm + 4*f2 + fb ) * h / 12;
        var ss = s1 + s2;
        var error = ( ss - s ) / 15;

        if ( Math.abs(error) < tolerance  || depth > maxIter ) return ss + error;
        else {
          var m = a + h/2;
          return adaptiveSimpson( a, m, fa, f1, fm, s1, tolerance/2, depth+1 )
                 + adaptiveSimpson( m, b, fm, f2, fb, s2, tolerance/2, depth+1 );
        }

      }

      var fa = f(a);
      var fm = f( (a+b)/2 );
      var fb = f(a);
      var s = ( fa + 4*fm + fb ) * (b-a) / 6;
      var depth = 0;

      return adaptiveSimpson( a, b, fa, fm, fb, s, tolerance, depth );

    default:

      throw( 'Unsupported method' );

  }

}


function discreteIntegral( values, step ) {

  // Euler-Maclaurin summation over fixed intervals

  var s = ( values[0] + values[ values.length - 1 ] ) / 2;

  for ( var i = 1 ; i < values.length - 1 ; i++ ) s += values[i];

  return s * step;

}


function polynomial( x, coefficients, derivative=false ) {

  // Horner's method with highest power coefficient first

  var p = coefficients[0];
  var q = 0;

  for ( var i = 1 ; i < coefficients.length ; i++ ) {
    if ( derivative ) q = p + q * x;
    p = coefficients[i] + p * x;
  }

  if ( derivative ) return { polynomial:p, derivative:q };
  else return p;

}


function findRoot( f, a, b, tolerance=1e-10, method='bisect' ) {

  switch( method ) {

    case 'bisect':

      var fa = f(a);
      var fb = f(b);

      if ( fa * f(b) >= 0 ) throw( 'Change of sign necessary for bisection' );

      var root, h;
      if ( fa < 0 ) {
        root = a;
        h = b - a;
      } else {
        root = b;
        h = a - b;
      }

      var maxIter = 100;
      for ( var i = 0; i < maxIter ; i++ ) {
        h /= 2;
        var mid = root + h;
        fmid = f(mid);
        if ( fmid <= 0 ) root = mid;
        if ( fmid === 0 || Math.abs(h) < tolerance ) return root;
      }

      throw( 'No root found for tolerance ' + tolerance );

    default:

      throw( 'Unsupported method' );

  }

}


function eigensystem( A, symmetric=true ) {

  if ( symmetric ) return tridiagonalQL( tridiagonalForm(A) );
  else throw( 'Unsupported system' );

}

// sourced from http://math.nist.gov/javanumerics/jama/
// no need to reinvent this wheel...

function tridiagonalForm( A ) {

  var n = A.length;
  var V = [];
  for ( var i = 0 ; i < n ; i++ ) V[i] = A[i].slice(); // deeper copy

  var d = vector( n );
  var e = vector( n );

  for ( var j = 0 ; j < n ; j++ ) d[j] = V[n-1][j];

  // Householder reduction to tridiagonal form
   
  for ( var i = n - 1 ; i > 0 ; i-- ) {
   
    // scale to avoid under/overflow
   
    var scale = 0;
    var h = 0;
    for ( var k = 0 ; k < i ; k++ ) scale += Math.abs(d[k]);

    if ( scale === 0 ) {
      e[i] = d[i-1];
      for ( var j = 0 ; j < i ; j++ ) {
        d[j] = V[i-1][j];
        V[i][j] = 0;
        V[j][i] = 0;
      }
    } else {

      // generate Householder vector

      for ( var k = 0 ; k < i ; k++ ) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      var f = d[i-1];
      var g = Math.sqrt(h);
      if ( f > 0 ) g = -g;
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for ( var j = 0; j < i; j++ ) e[j] = 0;

      // apply similarity transformation to remaining columns

      for ( var j = 0 ; j < i ; j++ ) {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for ( var k = j + 1 ; k <= i - 1 ; k++ ) {
          g += V[k][j] * d[k];
          e[k] += V[k][j] * f;
        }
        e[j] = g;
      }
      f = 0;
      for ( var j = 0 ; j < i ; j++ ) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      var hh = f / ( h + h );
      for ( var j = 0 ; j < i ; j++ ) {
        e[j] -= hh * d[j];
      }
      for ( var j = 0 ; j < i ; j++ ) {
        f = d[j];
        g = e[j];
        for ( var k = j ; k <= i - 1 ; k++ ) {
          V[k][j] -= f * e[k] + g * d[k];
        }
        d[j] = V[i-1][j];
        V[i][j] = 0;
      }

    }

    d[i] = h;

  }
   
  // accumulate transformations
   
  for ( var i = 0 ; i < n - 1 ; i++ ) {
    V[n-1][i] = V[i][i];
    V[i][i] = 1;
    var h = d[i+1];
    if ( h !== 0 ) {
      for ( var k = 0 ; k <= i ; k++ ) d[k] = V[k][i+1] / h;
      for ( var j = 0 ; j <= i ; j++ ) {
        var g = 0;
        for ( var k = 0 ; k <= i ; k++ ) g += V[k][i+1] * V[k][j];
        for ( var k = 0 ; k <= i ; k++ ) V[k][j] -= g * d[k];
      }
    }
    for ( var k = 0; k <= i; k++) V[k][i+1] = 0;
  }
  for ( var j = 0 ; j < n ; j++ ) {
    d[j] = V[n-1][j];
    V[n-1][j] = 0;
  }
  V[n-1][n-1] = 1;
  e[0] = 0;

  return { diagonal:d, offDiagonal:e, eigenvectors:V };

}


function tridiagonalQL( tridiagonalForm ) {

  var d = tridiagonalForm.diagonal;
  var n = d.length;
  var e = tridiagonalForm.offDiagonal;
  var V = tridiagonalForm.eigenvectors;

  function hypot( a, b) {
    var r;
    if ( Math.abs(a) > Math.abs(b) ) {
      r = b/a;
      r = Math.abs(a) * Math.sqrt( 1 + r*r );
    } else if (b != 0) {
      r = a/b;
      r = Math.abs(b) * Math.sqrt( 1 + r*r );
    } else r = 0;

    return r;

  }

  for ( var i = 1 ; i < n ; i++ ) e[i-1] = e[i];
  e[n-1] = 0;

  var f = 0;
  var tst1 = 0;
  var eps = Math.pow( 2, -52 );

  for ( var l = 0 ; l < n ; l++ ) {

    // find small subdiagonal element

    tst1 = Math.max( tst1, Math.abs(d[l]) + Math.abs(e[l]) );
    var m = l;
    while ( m < n ) {
      if ( Math.abs(e[m]) <= eps*tst1 ) break;
      m++;
    }

    // if m === l, d[l] is an eigenvalue, otherwise iterate

    if ( m > l ) {

      var iter = 0;
      do {

        iter = iter + 1;
        if ( iter > 1000 ) throw( 'Eigenvalues not converging...' );

        // compute implicit shift

        var g = d[l];
        var p = ( d[l+1] - g ) / ( 2 * e[l] );
        var r = hypot(p,1);
        if ( p < 0 ) r = -r;
        d[l] = e[l] / ( p + r );
        d[l+1] = e[l] * ( p + r );
        var dl1 = d[l+1];
        var h = g - d[l];
        for ( var i = l + 2 ; i < n ; i++ ) d[i] -= h;
        f = f + h;

        // implicit QL transformation

        p = d[m];
        var c = 1;
        var c2 = c;
        var c3 = c;
        var el1 = e[l+1];
        var s = 0;
        var s2 = 0;
        for ( var i = m - 1 ; i >= l ; i-- ) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = hypot(p,e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * ( c * g + s * d[i] );

          // accumulate transformation

          for ( var k = 0 ; k < n ; k++ ) {
            h = V[k][i+1];
            V[k][i+1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
   
        // check for convergence
   
      } while ( Math.abs(e[l]) > eps*tst1 );
    }
    d[l] = d[l] + f;
    e[l] = 0;
  }
   
  // sort eigenvalues and corresponding vectors

  for ( var i = 0 ; i < n - 1 ; i++ ) {
    var k = i;
    var p = d[i];
    for ( var j = i + 1 ; j < n ; j++ ) {
      if ( d[j] < p ) {
        k = j;
        p = d[j];
      }
    }
    if ( k != i ) {
      d[k] = d[i];
      d[i] = p;
      for ( var j = 0 ; j < n ; j++ ) {
        p = V[j][i];
        V[j][i] = V[j][k];
        V[j][k] = p;
      }
    }
  }

  return { eigenvalues:d, eigenvectors:V };

}


function hessenbergForm( A ) {

}


function luDecomposition( A, tolerance=1e-7 ) {

  var size = A.length;
  var LU = [];
  for ( var i = 0 ; i < size ; i++ ) LU[i] = A[i].slice(); // deeper copy

  var P = identity( size );
  pivots = 0;

  for ( var i = 0 ; i < size ; i++ ) {

    var maxValue = 0;
    var maxIndex = i;

    for ( var j = i ; j < size ; j++ ) {
      var element = Math.abs( LU[j][i] );
      if ( element > maxValue ) {
        maxValue = element;
        maxIndex = j;
      }
    }

    if ( maxValue < tolerance ) throw( 'Matrix is degenerate' );

    if ( maxIndex !== i ) {

      // pivot matrix rows
      var t = LU[i];
      LU[i] = LU[maxIndex];
      LU[maxIndex] = t;

      // pivot permutation rows
      var t = P[i];
      P[i] = P[maxIndex];
      P[maxIndex] = t;

      pivots++;

    }

    for ( var j = i + 1 ; j < size ; j++ ) {
      LU[j][i] /= LU[i][i];
      for ( var k = i + 1; k < size ; k++ )
        LU[j][k] -= LU[j][i] * LU[i][k];
    }
  }

  var L = identity( size );
  for ( var i = 1 ; i < size ; i++ )
    for ( var j = 0 ; j < i ; j++ ) L[i][j] = LU[i][j];

  var U = matrix( size );
  for ( var i = 0 ; i < size ; i++ )
    for ( var j = i ; j < size ; j++ ) U[i][j] = LU[i][j];

  return { L:L, U:U, P:P, pivots:pivots };

}

function luSolve( A, b ) {

  var size = A.length;
  var lu = luDecomposition(A);

  var x = vector( size );
  var y = vector( size );
  var pb = vector( size );

  for ( var i = 0 ; i < size ; i++ )
    for ( var j = 0 ; j < size ; j++ )
      pb[i] += lu.P[i][j] * b[j];

  // forward solve
  for ( var i = 0 ; i < size ; i++ ) {
    y[i] = pb[i];
    for ( var j = 0 ; j < i ; j++ ) y[i] -= lu.L[i][j] * y[j];
    y[i] /= lu.L[i][i];
  }

  // backward solve
  for ( var i = size - 1 ; i >= 0 ; i-- ) {
    x[i] = y[i];
    for ( var j = i + 1 ; j < size ; j++ ) x[i] -= lu.U[i][j] * x[j];
    x[i] /= lu.U[i][i];
  }

  return x;

}

function determinant( A ) {

  var lu = luDecomposition(A);

  var product = 1;
  for ( var i = 0 ; i < A.length; i++ ) product *= lu.U[i][i];

  return (-1)**lu.pivots * product;

}

function inverse( A ) {

  // calling luSolve for each column is not efficient
  //   but avoids code duplication

  var I = matrix( A.length );

  for ( var i = 0 ; i < A.length ; i++ ) {

    var b = vector( A.length );
    b[i] = 1;

    var x = luSolve( A, b );
    for ( var j = 0 ; j < A.length ; j++ ) I[j][i] = x[j];

  }

  return I;

}



function vector( size, value=0 ) {

  var v = [];
  for ( var i = 0 ; i < size ; i++ ) v.push( value );

  return v;

}

function matrix( rows, columns, value=0 ) {

  var columns = columns || rows;

  var m = [];
  for ( var i = 0 ; i < rows ; i++ ) {
    m.push( [] );
    for ( var j = 0 ; j < columns ; j++ ) m[i].push( value );
  }

  return m;

}

function identity( rows, value=1 ) {

  var m = matrix( rows );
  for ( var i = 0 ; i < rows ; i++ ) m[i][i] = value;

  return m;

}

function transpose( A ) {

  var T = matrix( A[0].length, A.length );

  for ( var i = 0 ; i < A.length ; i++ )
    for ( var j = 0 ; j < A[0].length ; j++ )
      T[j][i] = A[i][j];

  return T;

}

function matrixAdd( A, B ) {

  if ( !Array.isArray(A) && !Array.isArray(B) ) throw( 'No matrices' );
  if ( !Array.isArray(A) ) A = matrix( B.length, B[0].length, A );
  if ( !Array.isArray(B) ) B = matrix( A.length, A[0].length, B );

  var C = matrix( A.length, A[0].length, 0 );

  for ( var i = 0 ; i < A.length ; i++ )
    for ( var j = 0 ; j < A[0].length ; j++ )
      C[i][j] = add( A[i][j], B[i][j] );

  return C;

}

function matrixSub( A, B ) {

  if ( !Array.isArray(A) && !Array.isArray(B) ) throw( 'No matrices' );
  if ( !Array.isArray(A) ) A = matrix( B.length, B[0].length, A );
  if ( !Array.isArray(B) ) B = matrix( A.length, A[0].length, B );

  var C = matrix( A.length, A[0].length, 0 );

  for ( var i = 0 ; i < A.length ; i++ )
    for ( var j = 0 ; j < A[0].length ; j++ )
      C[i][j] = sub( A[i][j], B[i][j] );

  return C;

}

function matrixMul( A, B ) {

  if ( !Array.isArray(A) && !Array.isArray(B) ) throw( 'No matrices' );
  if ( !Array.isArray(A) ) A = identity( B.length, A );
  if ( !Array.isArray(B) ) B = identity( A[0].length, B );
  if ( A[0].length !== B.length ) throw( 'Incompatible matrices' );

  var C = matrix( A.length, B[0].length, 0 );

  for ( var i = 0 ; i < A.length ; i++ )
    for ( var j = 0 ; j < B[0].length ; j++ )
      for ( var k = 0 ; k < A[0].length ; k++ )
        C[i][j] = add( C[i][j], mul( A[i][k], B[k][j] ) );

  return C;

}

