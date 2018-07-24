
function jacobiTheta( n, x, q ) {

  switch( n ) {

    case 1:

      var s = 0;
      for ( var i = 0 ; i < 100 ; i++ ) s += (-1)**i * q**(i*i+i) * sin( (2*i+1) * x );
      if ( q < 0 ) return mul( 2 * s, pow( q, 1/4 ) );
      return 2 * q**(1/4) * s;

    case 2:

      var s = 0;
      for ( var i = 0 ; i < 100 ; i++ ) s += q**(i*i+i) * cos( (2*i+1) * x );
      if ( q < 0 ) return mul( 2 * s, pow( q, 1/4 ) );
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

      throw 'Undefined Jacobi theta index';

  }

}


function ellipticNome( m ) {

  if ( m >= 1 ) throw 'Unsupported elliptic parameter';

  if ( m < 0 ) return -exp( -pi * ellipticK( 1/(1-m) ) / ellipticK( m/(m-1) ) );

  return exp( -pi * ellipticK(1-m) / ellipticK(m) );

}


function sn( x, m ) {

  var q = ellipticNome(m);

  var t = x / jacobiTheta(3,0,q)**2;

  return jacobiTheta(3,0,q) / jacobiTheta(2,0,q)
         * jacobiTheta(1,t,q) / jacobiTheta(4,t,q);

}

function cn( x, m ) {

  var q = ellipticNome(m);

  var t = x / jacobiTheta(3,0,q)**2;

  return jacobiTheta(4,0,q) / jacobiTheta(2,0,q)
         * jacobiTheta(2,t,q) / jacobiTheta(4,t,q);

}

function dn( x, m ) {

  var q = ellipticNome(m);

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

