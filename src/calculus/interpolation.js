
function polynomial( x, coefficients, derivative=false ) {

  // Horner's method with highest power coefficient first

  var p = coefficients[0];
  var q = 0;

  for ( var i = 1 ; i < coefficients.length ; i++ ) {
    if ( derivative ) q = p + q * x;
    p = coefficients[i] + p * x;
  }

  if ( derivative ) return { polynomial:p, derivative:q };
  else return p;

}


function findRoot( f, a, b, tolerance=1e-10, method='bisect' ) {

  switch( method ) {

    case 'bisect':

      var fa = f(a);
      var fb = f(b);

      if ( fa * f(b) >= 0 ) throw 'Change of sign necessary for bisection';

      var root, h;
      if ( fa < 0 ) {
        root = a;
        h = b - a;
      } else {
        root = b;
        h = a - b;
      }

      var maxIter = 100;
      for ( var i = 0; i < maxIter ; i++ ) {
        h /= 2;
        var mid = root + h;
        fmid = f(mid);
        if ( fmid <= 0 ) root = mid;
        if ( fmid === 0 || Math.abs(h) < tolerance ) return root;
      }

      throw 'No root found for tolerance ' + tolerance;

    default:

      throw 'Unsupported method';

  }

}


function spline( points, value='function' ) {

  // adapted from gsl / cspline.c and reference therein

  var a = [], b = [], c = [], d = [];

  for ( var i = 0 ; i < points.length ; i++ ) a[i] = points[i][1];

  c[0] = 0;
  c[ points.length - 1 ] = 0;

  var A = matrix( points.length - 2 );
  var y = vector( points.length - 2 );

  function h( i ) { return points[i+1][0] - points[i][0]; }

  for ( var i = 0 ; i < A.length ; i++ ) {
    A[i][i] = 2 * ( h(i) + h(i+1) );
    y[i] = 3 * ( a[i+2] - a[i+1] ) / h(i+1) - 3 * ( a[i+1] - a[i] ) / h(i);
  }
  for ( var i = 1 ; i < A.length ; i++ ) {
    A[i][i-1] = h(i); 
    A[i-1][i] = h(i); 
  }

  var x = luSolve( A, y );

  for ( var i = 0 ; i < x.length ; i++ ) c[i+1] = x[i];

  for ( var i = 0 ; i < c.length - 1 ; i++ ) {
    b[i] = ( a[i+1] - a[i] ) / h(i) - ( c[i+1] + 2*c[i] ) * h(i) / 3;
    d[i] = ( c[i+1] - c[i] ) / 3 / h(i);
  }

  switch( value ) {

    case 'function':

      return function( x ) {

        for ( var i = 0 ; i < points.length ; i++ )
          if ( x === points[i][0] ) return a[i];

        for ( var i = 0 ; i < points.length - 1 ; i++ )
          if ( x > points[i][0] && x < points[i+1][0] ) {
            var xi = points[i][0];
            return a[i] + b[i] * ( x - xi )
                   + c[i] * ( x - xi )**2 + d[i] * ( x - xi )**3;
          }

      }

    case 'derivative':

      return function( x ) {

        for ( var i = 0 ; i < points.length ; i++ )
          if ( x === points[i][0] ) return b[i];

        for ( var i = 0 ; i < points.length - 1 ; i++ )
          if ( x > points[i][0] && x < points[i+1][0] ) {
            var xi = points[i][0];
            return b[i] + 2 * c[i] * ( x - xi ) + 3 * d[i] * ( x - xi )**2;
          }

      }

    case 'integral':

      return function( x ) {

        var sum = 0;

        function F( x, i ) {
          var xi = points[i][0];
          return a[i] * ( x - xi ) + b[i] * ( x - xi )**2 / 2
                   + c[i] * ( x - xi )**3 / 3 + d[i] * ( x - xi )**4 / 4;
        }

        for ( var i = 0 ; i < points.length - 1 ; i++ )
          if ( x < points[i+1][0] ) {
            sum += F( x, i ) - F( points[i][0], i );
            break;
          } else sum += F( points[i+1][0], i ) - F( points[i][0], i );

        return sum;

      }

    default:

      throw 'Unsupported value';

  }

}



