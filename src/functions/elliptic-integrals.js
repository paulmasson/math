
// Carlson symmetric integrals

function carlsonRC( x, y ) {

  if ( x < 0 || y < 0 || isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x);
    if ( !isComplex(y) ) y = complex(y);

    if ( x.re === y.re && x.im === y.im ) return inv( sqrt(x) );

    return div( arccos( sqrt( div(x,y) ) ), sqrt( sub(y,x) ) );

  }

  if ( x === y ) return 1 / Math.sqrt(x);

  if ( x < y )
    return Math.acos( Math.sqrt(x/y) ) / Math.sqrt(y-x);
  else
    return Math.acosh( Math.sqrt(x/y) ) / Math.sqrt(x-y);

}

function carlsonRD( x, y, z ) {

  return carlsonRJ( x, y, z, z );

}

function carlsonRF( x, y, z, tolerance=1e-10 ) {

  if ( isComplex(x) || isComplex(y) || isComplex(z) ) {

    var xm = x;
    var ym = y;
    var zm = z;

    var Am = A0 = div( add( x, y, z ), 3 );
    var Q = Math.pow( 3*tolerance, -1/6 )
            * Math.max( abs( sub(A0,x) ), abs( sub(A0,y) ), abs( sub(A0,z) ) );
    var g = .25;
    var pow4 = 1;
    var m = 0;

    while ( true ) {
      var xs = sqrt(xm);
      var ys = sqrt(ym);
      var zs = sqrt(zm);
      var lm = add( mul(xs,ys), mul(xs,zs), mul(ys,zs) );
      var Am1 = mul( add(Am,lm), g );
      xm = mul( add(xm,lm), g );
      ym = mul( add(ym,lm), g );
      zm = mul( add(zm,lm), g );
      if ( pow4 * Q < abs(Am) ) break;
      Am = Am1;
      m += 1;
      pow4 *= g;
    }

    var t = div( pow4, Am );
    var X = mul( sub(A0,x), t );
    var Y = mul( sub(A0,y), t );
    var Z = neg( add(X,Y) );
    var E2 = sub( mul(X,Y), mul(Z,Z) );
    var E3 = mul(X,Y,Z);

    return mul( pow( Am, -.5 ),
             add( 9240, mul(-924,E2), mul(385,E2,E2), mul(660,E3), mul(-630,E2,E3) ), 1/9240 );

  } else {

    if ( y === z ) return carlsonRC( x, y );
    if ( x === z ) return carlsonRC( y, x );
    if ( x === y ) return carlsonRC( z, x );

    // adapted from mpmath / elliptic.py

    var xm = x;
    var ym = y;
    var zm = z;

    var Am = A0 = (x + y + z) / 3;
    var Q = Math.pow( 3*tolerance, -1/6 )
            * Math.max( Math.abs(A0-x), Math.abs(A0-y), Math.abs(A0-z) );
    var g = .25;
    var pow4 = 1;
    var m = 0;

    while ( true ) {
      var xs = Math.sqrt(xm);
      var ys = Math.sqrt(ym);
      var zs = Math.sqrt(zm);
      var lm = xs*ys + xs*zs + ys*zs;
      var Am1 = (Am + lm) * g;
      xm = (xm + lm) * g;
      ym = (ym + lm) * g;
      zm = (zm + lm) * g;
      if ( pow4 * Q < Math.abs(Am) ) break;
      Am = Am1;
      m += 1;
      pow4 *= g;
    }

    var t = pow4 / Am;
    var X = (A0-x) * t;
    var Y = (A0-y) * t;
    var Z = -X-Y;
    var E2 = X*Y - Z**2;
    var E3 = X*Y*Z;

    return Math.pow( Am, -.5 )
           * ( 9240 - 924*E2 + 385*E2**2 + 660*E3 - 630*E2*E3 ) / 9240;

  }

}

function carlsonRG( x, y, z ) {

  return 1;

}

function carlsonRJ( x, y, z, p, tolerance=1e-10 ) {

  if ( isComplex(x) || isComplex(y) || isComplex(z) || isComplex(p) ) {

    var xm = x;
    var ym = y;
    var zm = z;
    var pm = p;

    var A0 = Am = div( add( x, y, z, mul(2,p) ), 5 );
    var delta = mul( sub(p,x), sub(p,y), sub(p,z) );
    var Q = Math.pow( .25*tolerance, -1/6 )
            * Math.max( abs( sub(A0,x) ), abs( sub(A0,y) ), abs( sub(A0,z) ), abs( sub(A0,p) ) );
    var m = 0;
    var g = .25;
    var pow4 = 1;
    var S = complex(0);

    while ( true ) {
      var sx = sqrt(xm);
      var sy = sqrt(ym);
      var sz = sqrt(zm);
      var sp = sqrt(pm);
      var lm = add( mul(sx,sy), mul(sx,sz), mul(sy,sz) );
      var Am1 = mul( add(Am,lm), g );
      xm = mul( add(xm,lm), g );
      ym = mul( add(ym,lm), g );
      zm = mul( add(zm,lm), g );
      pm = mul( add(pm,lm), g );
      var dm = mul( add(sp,sx), add(sp,sy), add(sp,sz) );
      var em = mul( delta, Math.pow( 4, -3*m ), inv(dm), inv(dm) );
      if ( pow4 * Q < abs(Am) ) break;
      var T = mul( carlsonRC( 1, add(1,em) ), pow4, inv(dm) );
      S = add( S, T );
      pow4 *= g;
      m += 1;
      Am = Am1;
    }

    var t = div( Math.pow( 2, -2*m ), Am );
    var X = mul( sub(A0,x), t );
    var Y = mul( sub(A0,y), t );
    var Z = mul( sub(A0,z), t );
    var P = div( add(X,Y,Z), -2 );
    var E2 = add( mul(X,Y), mul(X,Z), mul(Y,Z), mul(-3,P,P) );
    var E3 = add( mul(X,Y,Z), mul(2,E2,P), mul(4,P,P,P) );
    var E4 = mul( add( mul(2,X,Y,Z), mul(E2,P), mul(3,P,P,P) ), P );
    var E5 = mul(X,Y,Z,P,P);
    P = add( 24024, mul(-5148,E2), mul(2457,E2,E2), mul(4004,E3), mul(-4158,E2,E3), mul(-3276,E4), mul(2772,E5) );
    var v1 = mul( g**m, pow( Am, -1.5 ), P, 1/24024 );
    var v2 = mul(6,S);

    return add( v1, v2 );

  } else {

    // adapted from mpmath / elliptic.py

    var xm = x;
    var ym = y;
    var zm = z;
    var pm = p;

    var A0 = Am = (x + y + z + 2*p) / 5;
    var delta = (p-x) * (p-y) * (p-z);
    var Q = Math.pow( .25*tolerance, -1/6 )
            * Math.max( Math.abs(A0-x), Math.abs(A0-y), Math.abs(A0-z), Math.abs(A0-p) );
    var m = 0;
    var g = .25;
    var pow4 = 1;
    var S = 0;

    while ( true ) {
      var sx = Math.sqrt(xm);
      var sy = Math.sqrt(ym);
      var sz = Math.sqrt(zm);
      var sp = Math.sqrt(pm);
      var lm = sx*sy + sx*sz + sy*sz;
      var Am1 = (Am + lm) * g;
      xm = (xm + lm) * g;
      ym = (ym + lm) * g;
      zm = (zm + lm) * g;
      pm = (pm + lm) * g;
      var dm = (sp+sx) * (sp+sy) * (sp+sz);
      var em = delta * Math.pow( 4, -3*m ) / dm**2;
      if ( pow4 * Q < Math.abs(Am) ) break;
      var T = carlsonRC( 1, 1 + em ) * pow4 / dm;
      S += T;
      pow4 *= g;
      m += 1;
      Am = Am1;
    }

    var t = Math.pow( 2, -2*m ) / Am;
    var X = (A0-x) * t;
    var Y = (A0-y) * t;
    var Z = (A0-z) * t;
    var P = (-X-Y-Z) / 2;
    var E2 = X*Y + X*Z + Y*Z - 3*P**2;
    var E3 = X*Y*Z + 2*E2*P + 4*P**3;
    var E4 = ( 2*X*Y*Z + E2*P + 3*P**3 ) * P;
    var E5 = X*Y*Z*P**2;
    P = 24024 - 5148*E2 + 2457*E2**2 + 4004*E3 - 4158*E2*E3 - 3276*E4 + 2772*E5;
    var v1 = g**m * Math.pow( Am, -1.5 ) * P / 24024;
    var v2 = 6*S;

    return v1 + v2;

  }

}


// elliptic integrals

function ellipticF( x, m ) {

  if ( arguments.length === 1 ) {
    m = x;
    x = pi / 2;
  }

  if ( isComplex(x) || isComplex(m) ) {

    if ( !isComplex(x) ) x = complex(x);

    var period = complex(0);
    if ( Math.abs(x.re) > pi / 2 ) {
      var p = Math.round( x.re / pi );
      x.re = x.re - p * pi;
      period = mul( 2 * p, ellipticK( m ) );
    }

    return add( mul( sin(x), carlsonRF( mul(cos(x),cos(x)), sub( 1, mul(m,sin(x),sin(x)) ), 1 ) ), period );

  } else {

    if ( m > 1 && x > Math.asin( 1 / Math.sqrt(m) ) ) return ellipticF( C(x), m );

    var period = 0;
    if ( Math.abs(x) > pi / 2 ) {
      var p = Math.round( x / pi );
      x = x - p * pi;
      period = 2 * p * ellipticK( m );
    }

    return sin(x) * carlsonRF( cos(x)**2, 1 - m * sin(x)**2, 1 ) + period;

  }

}

function ellipticK( m ) {

  return ellipticF( m );

}

function ellipticE( x, m ) {

  if ( arguments.length === 1 ) {
    m = x;
    x = pi / 2;
  }

  if ( isComplex(x) || isComplex(m) ) {

    if ( !isComplex(x) ) x = complex(x);

    var period = complex(0);
    if ( Math.abs(x.re) > pi / 2 ) {
      var p = Math.round( x.re / pi );
      x.re = x.re - p * pi;
      period = mul( 2 * p,  ellipticE( m ) );
    }

    return add( mul( sin(x), carlsonRF( mul(cos(x),cos(x)), sub( 1, mul(m,sin(x),sin(x)) ), 1 ) ),
                mul( -1/3, m, pow(sin(x),3), carlsonRD( mul(cos(x),cos(x)), sub( 1, mul(m,sin(x),sin(x)) ), 1 ) ),
                period );

  } else {

    if ( m > 1 && x > Math.asin( 1 / Math.sqrt(m) ) ) return ellipticE( C(x), m );

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

}

function ellipticPi( n, x, m ) {

  if ( arguments.length === 2 ) {
    m = x;
    x = pi / 2;
  }

  if ( isComplex(n) || isComplex(x) || isComplex(m) ) {

    if ( !isComplex(x) ) x = complex(x);

    var period = complex(0);
    if ( Math.abs(x.re) > pi / 2 ) {
      var p = Math.round( x.re / pi );
      x.re = x.re - p * pi;
      period = mul( 2 * p, ellipticPi( n, m ) );
    }

    return add( mul( sin(x), carlsonRF( mul(cos(x),cos(x)), sub( 1, mul(m,sin(x),sin(x)) ), 1 ) ),
                mul( 1/3, n, pow(sin(x),3),
                  carlsonRJ( mul(cos(x),cos(x)), sub( 1, mul(m,sin(x),sin(x)) ), 1, sub( 1, mul(n,sin(x),sin(x)) ) ) ),
                period );

  } else {

    if ( n > 1 && x > Math.asin( 1 / Math.sqrt(n) ) ) return ellipticPi( n, C(x), m );

    if ( m > 1 && x > Math.asin( 1 / Math.sqrt(m) ) ) return ellipticPi( n, C(x), m );

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

}


function jacobiZeta( x, m ) {

  // using definition matching elliptic integrals
  // alternate definition replaces x with am(x,m)

  return sub( ellipticE( x, m ), mul( ellipticF(x,m), ellipticE(m), inv( ellipticK(m) ) ) );

}

