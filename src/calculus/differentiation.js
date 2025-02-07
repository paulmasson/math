
function diff( f, x, n=1, options={} ) {

  if ( !Number.isInteger(n) || n < 0 )
    throw Error( 'Only positive integer orders for differentiation' );

  if ( isComplex(x) || isComplex(f(x)) ) {

    if ( !isComplex(f(x)) ) throw Error( 'Function must handle complex math' );

    var absX = abs(x);
    var normed = absX === 0 ? complex(1) : div( x, absX );

    var real = diff( t => f( mul(normed,t) ).re, absX, n, options );
    var imag = diff( t => f( mul(normed,t) ).im, absX, n, options );

    return div( complex( real, imag ), normed );

  }

  var method = 'method' in options ? options.method : 'ridders';
  var epsilon = 'epsilon' in options ? options.epsilon : 1e-5;

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
      var h = epsilon**(1/(n+2));
      return difference();

    case 'ridders':

      var h = epsilon**(1/(n+2));
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

      throw Error( 'Unsupported differentiation method' );

  }

}

var D = diff;


function taylorSeries( f, x0, terms=5 ) {

  var c = [ f(x0) ];
  for ( var i = 1 ; i < terms ; i++ ) c.push( diff( f, x0, i ) );

  return function( x ) {

    var s = 0;
    for ( var i = 0 ; i < c.length ; i++ )
      s = add( s, mul( c[i], pow( sub(x,x0), i ), 1/factorial(i) ) );
    return s;

  }

}


function gradient( f, point ) {

  if ( f.length !== point.length ) throw Error( 'Gradient point length differs from function' );

  var result = [];

  for ( var i = 0 ; i < point.length ; i++ ) {

    var a = point.slice();

    result.push( diff( x => { a[i] = x; return f.apply( null, a ); }, a[i] ) );

  }

  return result;

}

function findExtremum( f, point, options={} ) {

  if ( !Array.isArray(point) ) point = [ point ];

  var sign = options.findMaximum ? 1 : -1;
  var tolerance = 'tolerance' in options ? options.tolerance : 1e-10;

  var maxIter = 1e4;
  var gamma = .01 * sign;
  var grad, step, test;

  for ( var i = 0 ; i < maxIter ; i++ ) {

    grad = gradient( f, point );
    test = true;

    for ( var j = 0 ; j < point.length ; j++ ) {
      step = gamma * grad[j];
      point[j] += step;
      test = test && step < tolerance;
    }

    if ( test )
     if ( point.length === 1 ) return point[0];
     else return point;

  }

  throw Error( 'No extremum found for tolerance ' + tolerance );

}

