
function diff( f, x, n=1, method='ridders' ) {

  if ( isComplex(x) || isComplex(f(x)) ) {

    if ( !isComplex(f(x)) ) throw 'Function must handle complex math';

    var real = diff( t => f( mul(x,t) ).re, 1, n, method );
    var imag = diff( t => f( mul(x,t) ).im, 1, n, method );

    return div( complex( real, imag ), x );

  }

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

      throw 'Unsupported method';

  }

}

var D = diff;


function gradient( f, point ) {

  if ( f.length !== point.length ) throw 'Gradient point length differs from function';

  var result = [];

  for ( var i = 0 ; i < point.length ; i++ ) {

    var a = [].concat( point );

    result.push( diff( x => { a[i] = x; return f.apply( null, a ); }, a[i] ) );

  }

  return result;

}

