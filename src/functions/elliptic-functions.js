
function jacobiTheta( n, x, q ) {

  if ( abs(q) >= 1 ) throw 'Unsupported elliptic nome';

  switch( n ) {

    case 1:

      if ( q < 0 || isComplex(x) || isComplex(q) ) {

        var s = complex(0);
        for ( var i = 0 ; i < 100 ; i++ )
          s = add( s, mul( (-1)**i, pow( q, i*i+i ), sin( mul(2*i+1,x) ) ) );
        return mul( 2, pow( q, 1/4 ), s );

      }

      var s = 0;
      for ( var i = 0 ; i < 100 ; i++ ) s += (-1)**i * q**(i*i+i) * sin( (2*i+1) * x );
      return 2 * q**(1/4) * s;

    case 2:

      if ( q < 0 || isComplex(x) || isComplex(q) ) {

        var s = complex(0);
        for ( var i = 0 ; i < 100 ; i++ )
          s = add( s, mul( pow( q, i*i+i ), cos( mul(2*i+1,x) ) ) );
        return mul( 2, pow( q, 1/4 ), s );

      }

      var s = 0;
      for ( var i = 0 ; i < 100 ; i++ ) s += q**(i*i+i) * cos( (2*i+1) * x );
      return 2 * q**(1/4) * s;

    case 3:

      if ( isComplex(x) || isComplex(q) ) {

        var s = complex(0);
        for ( var i = 1 ; i < 100 ; i++ )
          s = add( s, mul( pow( q, i*i ), cos( mul(2*i,x) ) ) );
        return add( 1, mul(2,s) );

      }

      var s = 0;
      for ( var i = 1 ; i < 100 ; i++ ) s += q**(i*i) * cos( 2*i * x );
      return 1 + 2 * s;

    case 4:

      if ( isComplex(x) || isComplex(q) ) {

        var s = complex(0);
        for ( var i = 1 ; i < 100 ; i++ )
          s = add( s, mul( pow( neg(q), i*i ), cos( mul(2*i,x) ) ) );
        return add( 1, mul(2,s) );

      }

      var s = 0;
      for ( var i = 1 ; i < 100 ; i++ ) s += (-q)**(i*i) * cos( 2*i * x );
      return 1 + 2 * s;

    default:

      throw 'Undefined Jacobi theta index';

  }

}


function ellipticNome( m ) {

  if ( isComplex(m) ) return exp( div( mul( -pi, ellipticK( sub(1,m) ) ), ellipticK(m) ) );

  if ( m > 1 ) return ellipticNome( complex(m) );

  if ( m < 0 ) return -exp( -pi * ellipticK( 1/(1-m) ) / ellipticK( m/(m-1) ) );

  return exp( -pi * ellipticK(1-m) / ellipticK(m) );

}


function sn( x, m ) {

  var q = ellipticNome(m);

  if ( m > 1 || isComplex(x) || isComplex(m) ) {

    var t = div( x, pow( jacobiTheta(3,0,q), 2 ) );

    return mul( div( jacobiTheta(3,0,q), jacobiTheta(2,0,q) ),
                div( jacobiTheta(1,t,q), jacobiTheta(4,t,q) ) );

  }

  var t = x / jacobiTheta(3,0,q)**2;

  if ( m < 0 )
    return jacobiTheta(3,0,q) / jacobiTheta(4,t,q)
           * div( jacobiTheta(1,t,q), jacobiTheta(2,0,q) ).re;

  return jacobiTheta(3,0,q) / jacobiTheta(2,0,q)
         * jacobiTheta(1,t,q) / jacobiTheta(4,t,q);

}

function cn( x, m ) {

  var q = ellipticNome(m);

  if ( m > 1 || isComplex(x) || isComplex(m) ) {

    var t = div( x, pow( jacobiTheta(3,0,q), 2 ) );

    return mul( div( jacobiTheta(4,0,q), jacobiTheta(2,0,q) ),
                div( jacobiTheta(2,t,q), jacobiTheta(4,t,q) ) );

  }

  var t = x / jacobiTheta(3,0,q)**2;

  if ( m < 0 )
    return jacobiTheta(4,0,q) / jacobiTheta(4,t,q)
           * div( jacobiTheta(2,t,q), jacobiTheta(2,0,q) ).re;

  return jacobiTheta(4,0,q) / jacobiTheta(2,0,q)
         * jacobiTheta(2,t,q) / jacobiTheta(4,t,q);

}

function dn( x, m ) {

  var q = ellipticNome(m);

  if ( m > 1 || isComplex(x) || isComplex(m) ) {

    var t = div( x, pow( jacobiTheta(3,0,q), 2 ) );

    return mul( div( jacobiTheta(4,0,q), jacobiTheta(3,0,q) ),
                div( jacobiTheta(3,t,q), jacobiTheta(4,t,q) ) );

  }

  var t = x / jacobiTheta(3,0,q)**2;

  return jacobiTheta(4,0,q) / jacobiTheta(3,0,q)
         * jacobiTheta(3,t,q) / jacobiTheta(4,t,q);

}

function am( x, m ) {

  var K = ellipticK(m);
  var n = Math.floor( x / 2 / K );
  x = x - 2 * n * K;

  return Math.atan2( sn(x,m), cn(x,m) ) + n * pi;

}

