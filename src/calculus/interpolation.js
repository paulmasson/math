
function polynomial( x, coefficients, derivative=false ) {

  // Horner's method with highest power coefficient first

  var p = coefficients[0];
  var q = 0;

  for ( var i = 1 ; i < coefficients.length ; i++ ) {
    if ( derivative ) q = add( p, mul( q, x ) );
    p = add( coefficients[i], mul( p, x ) );
  }

  if ( derivative ) return { polynomial: p, derivative: q };

  return p;

}


function partialBell( n, k, argumentArray ) {

  if ( n === 0 && k === 0 ) return 1;

  if ( n === 0 || k === 0 ) return 0;

  // evaluate recursively
  var s = 0;
  var p = 1;

  for ( var i = 1 ; i <= n - k + 1 ; i++ ) {

    s += p * argumentArray[i-1] * partialBell( n-i, k-1, argumentArray );
    p *= ( n - i ) / i;

  }

  return s;

}


function findRoot( f, start, options={} ) {

  var tolerance = 'tolerance' in options ? options.tolerance : 1e-10;
  var maxIter = 100;

  if ( Array.isArray(f) ) {

    if ( f.length !== start.length )
      throw Error( 'Mismatch between equations and starting point for root' );

    var root = start;

    for ( var i = 0; i < maxIter ; i++ ) {

      var J = [], F = [];

      for ( var j = 0 ; j < root.length ; j++ ) {
        J.push( gradient( f[j], root ) );
        F.push( f[j].apply( null, root ) );
      }

      var delta = luSolve( J, F );

      for ( var j = 0 ; j < root.length ; j++ ) root[j] -= delta[j];

      if ( delta.every( d => Math.abs(d) < tolerance ) ) return root;

    }

    throw Error( 'No root found for tolerance ' + tolerance );

  }

  if ( !Array.isArray(start) && !options.method ) options.method = 'newton';

  var method = 'method' in options ? options.method : 'bisect';

  switch( method ) {

    case 'bisect':

      var a = start[0];
      var b = start[1];

      var fa = f(a);
      var fb = f(b);

      if ( fa * f(b) >= 0 ) throw Error( 'Change of sign necessary for bisection' );

      var root, h;
      if ( fa < 0 ) {
        root = a;
        h = b - a;
      } else {
        root = b;
        h = a - b;
      }

      for ( var i = 0; i < maxIter ; i++ ) {
        h /= 2;
        var mid = root + h;
        fmid = f(mid);
        if ( fmid <= 0 ) root = mid;
        if ( fmid === 0 || Math.abs(h) < tolerance ) return root;
      }

      throw Error( 'No root found for tolerance ' + tolerance );

    case 'newton':

      var root = start;

      if ( isComplex( f(root) ) ) {

        for ( var i = 0; i < maxIter ; i++ ) {
          var delta = div( f(root), diff( f, root ) );
          root = sub( root, delta );
          if ( abs(delta) < tolerance ) return root;
        }

      } else {

        for ( var i = 0; i < maxIter ; i++ ) {
          var delta = f(root) / diff( f, root );
          root -= delta;
          if ( Math.abs(delta) < tolerance ) return root;
        }

      }

      throw Error( 'No root found for tolerance ' + tolerance );

    default:

      throw Error( 'Unsupported root finding method' );

  }

}


function findRoots( f, point, tolerance=1e-10 ) {

  console.log( 'Renamed to findRoot' );

}


function spline( points, value='function', tolerance=1e-10 ) {

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

        if ( x < points[0][0] || x > points[points.length-1][0] )
          throw Error( 'Argument outside spline input domain' );

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

        if ( x < points[0][0] || x > points[points.length-1][0] )
          throw Error( 'Argument outside spline input domain' );

        // method does not define b[points.length-1] so fudge endpoint
        if ( x === points[points.length-1][0] ) x -= tolerance;

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

        if ( x < points[0][0] || x > points[points.length-1][0] )
          throw Error( 'Argument outside spline input domain' );

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

      throw Error( 'Unsupported spline value' );

  }

}

function padeApproximant( f, n, d, center=0 ) {

  if ( !isPositiveInteger(n) || !isPositiveInteger(d) )
    throw Error( 'Pade indices must be positive integers' );

  var c = [];

  if ( typeof f === 'function' ) {
    c.push( f(center) );
    for ( var i = 1 ; i <= n+d ; i++ ) c.push( diff( f, center, i ) / factorial(i) );
  } else {
    if ( f.length < n+d+1 ) throw Error( 'Need n+d+1 Pade coefficients' );
    c = f;
  }

  var M = matrix( d ), v = vector( d );

  for ( var i = 0 ; i < d ; i++ ) {
    v[i] = -c[ n + 1 + i ];
    // avoid negative index leading to undefined entry
    for ( var j = 0 ; j < d && j <= n+i ; j++ ) M[i][j] = c[ n + i - j ];
  }

  // need to understand why Mathematica returns results for singular matrices
  try { var b = luSolve( M, v ); }
  catch { throw Error( 'Singular Pade matrix encountered' ); };
  b.unshift( 1 );

  var a = [ c[0] ];

  for ( var i = 1 ; i <= n ; i++ ) {
    var s = c[i];
    // cannot exceed i when d > n
    for( var j = 1 ; j <= d && j <= i ; j++ ) s += b[j] * c[i-j];
    a.push( s );
  }

  return { N: chop(a), D: chop(b) };

}

