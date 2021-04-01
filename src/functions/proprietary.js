
// This file contains proprietary functions defined on analyticphysics.com
// Before each is the title of a presentation describing the function


// A Generalized Lambert Function of Two Arguments

function doubleLambert( n, x, y, tolerance=1e-10 ) {

  if ( isZero(y) ) return lambertW(n,x);
  if ( isZero(x) ) return neg( lambertW( -n, neg(y) ) );

  function asymptotic( n, x, y ) {

    return add( log(sqrt(x)), neg(log(sqrt(neg(y)))), complex(0,n*pi) );

  }

  function test( n, x, y ) {

    return div( asymptotic(n,x,y), add( 1, mul(2*(-1)**n,sqrt(x),sqrt(neg(y))) ) );

  }

  function start( n, x, y, tolerance=1e-3 ) {

    var a = asymptotic(n,x,y), w1, w2;

    var testValue = .9;

    if ( abs(test(n,x,y)) < testValue ) return a;

    if ( n === 0 ) console.log( 'Using Lambert W on principal branch' );

    if ( n & 1 ) {

      if ( n > 0 ) {
        w1 = lambertW( (n+1)/2, x, tolerance );
        w2 = neg( lambertW( -(n+1)/2, neg(y), tolerance ) );
      } else {
        w1 = lambertW( (n-1)/2, x, tolerance );
        w2 = neg( lambertW( -(n-1)/2, neg(y), tolerance ) );
      }

    } else {

      w1 = lambertW( n/2, x, tolerance );
      w2 = neg( lambertW( -n/2, neg(y), tolerance ) );

    }

    if ( Math.abs( w1.im - a.im ) < Math.abs( w2.im - a.im ) )
      return w1;
    else
      return w2;

  }

  var maxIter = 100;
  var root = start(n,x,y);

  for ( var i = 0; i < maxIter ; i++ ) {
    var N = add( mul(x,exp(neg(root))), mul(y,exp(root)), neg(root) );
    var D = add( mul(-1,x,exp(neg(root))), mul(y,exp(root)), -1 );
    var delta = div( N, D );
    root = sub( root, delta );
    if ( abs(delta) < tolerance ) return root;
  }

  throw Error( 'No double Lambert root found for x = ' + JSON.stringify(x)
               + ' and y = ' + JSON.stringify(y) );

}

