
function zeta( x ) {

  // Borwein algorithm

  var n = 14; // from error bound for tolerance of 1e-10 

  function d( k ) {

    var d = 0;

    for ( var i = 0 ; i <= k ; i++ )
      d += 4**i * factorial( n+i-1 ) / factorial( n-i ) / factorial( 2*i );

    return n * d;

  }

  if ( isComplex(x) ) {

    // functional equation
    if ( x.re < 0 )
      return mul( pow(2,x), pow(pi,sub(x,1)), sin( mul(pi/2,x) ), gamma( sub(1,x) ), zeta( sub(1,x) ) );

    var s = complex(0);

    for ( var k = 0 ; k < n ; k++ )
      s = add( s, div( (-1)**k * ( d(k) - d(n) ), pow( k+1, x ) ) );

    return div( div( s, -d(n) ), sub( 1, pow( 2, sub(1,x) ) ) );

  }

  // functional equation
  if ( x < 0 ) return 2**x * pi**(x-1) * sin(pi*x/2) * gamma(1-x) * zeta(1-x);

  var s = 0;

  for ( var k = 0 ; k < n ; k++ )
    s += (-1)**k * ( d(k) - d(n) ) / (k+1)**x;

  return -s / d(n) / ( 1 - 2**(1-x) );

}

function eta( x ) { return ( 1 - 2**(1-x) ) * zeta(x); }

