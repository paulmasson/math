
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

  if ( Number.isInteger(n) && n >= 0 ) return polynomial( x, coefficients(n) );

  if ( isComplex(n) || isComplex(x) ) {

    var a = div( n, -2 );
    var b = div( sub(1,n), 2 );

    var s = sub( div( hypergeometric1F1( a, 1/2, pow(x,2) ), gamma( b ) ),
                 mul( 2, x, div( hypergeometric1F1( b, 3/2, pow(x,2) ), gamma( a ) ) ) );

    return mul( pow(2,n), sqrt(pi), s );

  }

  var s = hypergeometric1F1( -n/2, 1/2, x**2 ) / gamma( (1-n)/2 )
          - 2 * x * hypergeometric1F1( (1-n)/2, 3/2, x**2 ) / gamma( -n/2 );

  return 2**n * sqrt(pi) * s;

}


function laguerre( n, a, x ) {

  if ( arguments.length < 3 ) {
    x = a;
    a = 0
  }

  return mul( binomial( add(n,a), n ), hypergeometric1F1( neg(n), add(a,1), x ) ); 

}

