
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

