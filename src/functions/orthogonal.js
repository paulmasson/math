
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

  if ( isComplex(n) || isComplex(x) ) {

    if ( !isComplex(n) ) n = complex(n);
    if ( isPositiveIntegerOrZero(n) ) return polynomial( x, coefficients(n.re) );

    var a = div( n, -2 );
    var b = div( sub(1,n), 2 );

    var s = sub( div( hypergeometric1F1( a, .5, pow(x,2) ), gamma( b ) ),
                 mul( 2, x, div( hypergeometric1F1( b, 1.5, pow(x,2) ), gamma( a ) ) ) );

    return mul( pow(2,n), sqrt(pi), s );

  }

  if ( isPositiveIntegerOrZero(n) ) return polynomial( x, coefficients(n) );

  var s = hypergeometric1F1( -n/2, .5, x**2 ) / gamma( (1-n)/2 )
          - 2 * x * hypergeometric1F1( (1-n)/2, 1.5, x**2 ) / gamma( -n/2 );

  return 2**n * sqrt(pi) * s;

}


function laguerre( n, a, x ) {

  // explict recursion unnecessary: hypergeometric series handles integers

  if ( arguments.length < 3 ) {
    x = a;
    a = 0
  }

  return mul( binomial( add(n,a), n ), hypergeometric1F1( neg(n), add(a,1), x ) ); 

}


function chebyshevT( n, x ) {

  return cos( mul( n, arccos(x) ) );

}

function chebyshevU( n, x ) {

  return div( sin( mul( add(n,1), arccos(x) ) ), sin( arccos(x) ) );

}


function legendreP( l, m, x, renormalized=false ) {

  if ( arguments.length < 3 ) {
    x = m;
    m = 0;
  }

  if ( Number.isInteger(l) && Number.isInteger(m) && Math.abs(x) <= 1 ) {

    var mm = Math.abs(m);
    if ( mm > l ) throw Error( 'Invalid spherical harmonic indices' );

    if ( !renormalized ) {
      var norm = 1;
      for ( var i = l-m+1 ; i <= l+m ; i++ ) norm *= i;
      norm = Math.sqrt( 4 * pi * norm / (2*l+1) );
    }

    var legendre1 = (-1)**mm * Math.sqrt( (2*mm+1) / 4 / pi / factorial(2*mm) )
                    * factorial2( 2*mm-1 ) * ( 1 - x*x )**(mm/2);

    if ( mm === l ) 
      if ( renormalized ) return legendre1;
      else return norm * legendre1;

    var ll = mm + 1;
    var factor1 = Math.sqrt( 2*mm+3 );
    var legendre2 = factor1 * x * legendre1;

    if ( ll === l )
      if ( renormalized ) return legendre2;
      else return norm * legendre2;

    while ( ll < l ) {
      ll++
      var factor2 = Math.sqrt( ( 4*ll*ll - 1 ) / ( ll*ll - mm*mm ) );
      var legendre3 = factor2 * ( x*legendre2 - legendre1/factor1 );
      legendre1 = legendre2;
      legendre2 = legendre3;
      factor1 = factor2;
    }

    if ( renormalized ) return legendre3;
    else return norm * legendre3;

  }

  // dlmf.nist.gov/14.3.5
  if ( isPositiveInteger(m) )
    return mul( pow(-1,m), inv( gamma( add(m,1) ) ),
                gamma( add(l,m,1) ), inv( gamma( add(l,neg(m),1) ) ),
                pow( add(1,x), div(m,-2) ), pow( sub(1,x), div(m,2) ),
                hypergeometric2F1( neg(l), add(l,1), add(m,1), div(sub(1,x),2) ) );

  return mul( inv( gamma( sub(1,m) ) ),
              pow( add(1,x), div(m,2) ), pow( sub(1,x), div(m,-2) ),
              hypergeometric2F1( neg(l), add(l,1), sub(1,m), div(sub(1,x),2) ) );

}

function sphericalHarmonic( l, m, theta, phi ) {

  var renormalizedLegendre = legendreP( l, m, cos(theta), true );

  return mul( Math.sign(m)**m, renormalizedLegendre, exp( complex(0,m*phi) ) );

}

function legendreQ( l, m, x ) {

  if ( arguments.length < 3 ) {
    x = m;
    m = 0;
  }

  function difference( t ) {

    var t1 = mul( cos( mul(pi,t) ), legendreP(l,t,x) );

    var t2 = mul( gamma( add(l,t,1) ), inv( gamma( add(l,neg(t),1) ) ), legendreP(l,neg(t),x) );

    return sub( t1, t2 );

  }

  // l'Hopital's rule decent for small m, more accurate might be
  // functions.wolfram.com/HypergeometricFunctions/LegendreQ2General/26/01/02/0005/

  if ( isInteger(m) ) {

    // legendreP is pure real outside unit circle for even integers
   if ( abs(x) > 1 && !isComplex(m) && !( m & 1 ) ) m = complex(m);

    return mul( .5, pow(-1,m), diff( t => difference(t), m ) );

  }

  return mul( pi/2, inv( sin(mul(pi,m)) ), difference(m) );

}

