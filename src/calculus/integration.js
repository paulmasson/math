
function integrate( f, [a,b], options={} ) {

  var method = 'method' in options ? options.method : 'adaptive-simpson';
  var tolerance = 'tolerance' in options ? options.tolerance : 1e-10;

  if ( isComplex(a) || isComplex(b) || isComplex(f(a)) ) {

    if ( !isComplex(a) ) a = complex(a);
    if ( !isComplex(b) ) b = complex(b);

    if ( a.re === b.re && a.im === b.im ) return complex(0);

    if ( !isComplex(f(a)) ) throw Error( 'Function must handle complex math' );

    function lerp( t ) { return add( mul( sub(b,a), t ), a ); }

    var real = integrate( t => f( lerp(t) ).re, [0,1], options );
    var imag = integrate( t => f( lerp(t) ).im, [0,1], options );

    return mul( sub(b,a), complex( real, imag ) );

  }

  if ( a === b ) return 0;

  if ( options.avoidEndpoints )
    if ( a < b ) { a += tolerance; b -= tolerance; }
    else { a -= tolerance; b += tolerance; }

  function nextEulerIteration() {

      h /= 2;
      var x = a + h;
      while ( x < b ) {
        // only add new function evaluations
        s += f(x);
        x += 2*h;
      }

  }

  switch( method ) {

    case 'euler-maclaurin':

      // Euler-Maclaurin summation formula

      var maxIter = 50;

      var h = ( b - a ) / 2;
      var s = ( f(a) + f(b) ) / 2 + f( (a+b)/2 );
      var result = h * s;
      var previous = result;

      for ( var i = 0 ; i < maxIter ; i++ ) {

        nextEulerIteration();
        result = h * s;
        if ( Math.abs( result - previous ) < tolerance * Math.abs(previous) ) return result;
        previous = result;

      }

      throw Error( 'Maximum interations reached' );

    case 'romberg':

      var error = Number.MAX_VALUE;
      var maxIter = 30;

      var h = ( b - a ) / 2;
      var s = ( f(a) + f(b) ) / 2 + f( (a+b)/2 );
      var result = h * s;

      var d = [];
      for ( var i = 0 ; i < maxIter ; i++ ) d.push( [] );

      // Richardson extrapolation of Euler-Maclaurin trapezoids
      d[0][0] = result;

      for ( var i = 1 ; i < maxIter ; i++ ) {

        nextEulerIteration();
        d[0][i] = h * s;

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

    case 'adaptive-simpson':

      // algorithm by Charles Collins

      var maxIter = 50;

      function adaptiveSimpson( a, b, fa, fm, fb, s, tolerance, depth ) {

        var h = b - a;
        var f1 = f( a + h/4 );
        var f2 = f( b - h/4 )

        if ( isNaN(f1) || isNaN(f2) ) throw Error( 'NaN encountered in integration' );

        var s1 = ( fa + 4*f1 + fm ) * h / 12;
        var s2 = ( fm + 4*f2 + fb ) * h / 12;
        var ss = s1 + s2;
        var error = ( ss - s ) / 15;

        if ( Math.abs(error) < tolerance  || depth > maxIter ) return ss + error;
        else {
          var m = a + h/2;
          return adaptiveSimpson( a, m, fa, f1, fm, s1, tolerance/2, depth+1 )
                 + adaptiveSimpson( m, b, fm, f2, fb, s2, tolerance/2, depth+1 );
        }

      }

      var fa = f(a);
      var fm = f( (a+b)/2 );
      var fb = f(a);
      var s = ( fa + 4*fm + fb ) * (b-a) / 6;
      var depth = 0;

      return adaptiveSimpson( a, b, fa, fm, fb, s, tolerance, depth );

    case 'tanh-sinh':

      // based on Borwein & Bailey, Experimentation in Mathematics

      var m = 10;
      var h = 1 / 2**m;
      var x = [], w = [];

      for ( var k = 0 ; k <= 20 * 2**m ; k++ ) {
        var t = k * h;
        x[k] = Math.tanh( Math.PI/2 * Math.sinh(t) );
        w[k] = Math.PI/2 * Math.cosh(t) / Math.cosh( Math.PI/2 * Math.sinh(t) )**2;
        if ( Math.abs(1-x[k]) < tolerance ) break;
      }

      var nt = k;
      var sum = 0;

      // rescale [a,b] to [-1,1]
      var len = ( b - a ) / 2;
      var mid = ( b + a ) / 2;

      for ( var k = 1 ; k <= m ; k++ ) {
        for ( var i = 0 ; i < nt ; i += 2**(m-k) ) {
          if ( i % 2**(m-k+1) !== 0 || k === 1 ) {
            if ( i === 0 ) sum += w[0] * f( mid );
            else sum += w[i] * ( f( mid - len*x[i] ) + f( mid + len*x[i] ) );
          }
        }
      }

      return len * h * sum;

    case 'gaussian':

      // based on Borwein & Bailey, Experimentation in Mathematics

      var m = 10;
      var x = [], w = [];

      var n = 3 * 2**m;

      for ( var j = 1 ; j <= n/2 ; j++ ) {

        var r = Math.cos( Math.PI * (j-.25) / (n+.5) );

        while ( true ) {

          var t1 = 1, t2 = 0;

          for ( var j1 = 1 ; j1 <= n ; j1++ ) {
            t3 = t2;
            t2 = t1;
            t1 = ( (2*j1-1) * r * t2 - (j1-1) * t3 ) / j1;
          }

          var t4 = n * ( r*t1 - t2 ) / ( r**2 - 1 );
          var delta = t1 / t4;
          r -= delta;

          if ( Math.abs( delta ) < tolerance ) break;

        }

        x[j] = r;
        w[j] = 2 / ( 1 - r**2 ) / t4**2

      }

      // rescale [a,b] to [-1,1]
      var len = ( b - a ) / 2;
      var mid = ( b + a ) / 2;

      var sum = 0;

      for ( var j = 1 ; j <= n/2 ; j++ ) 
        sum += w[j] * ( f( mid - len*x[j] ) + f( mid + len*x[j] ) );

      return len * sum;

    default:

      throw Error( 'Unsupported integration method' );

  }

}


function discreteIntegral( values, step ) {

  // composite trapezoidal rule

  var s = ( values[0] + values[ values.length - 1 ] ) / 2;

  for ( var i = 1 ; i < values.length - 1 ; i++ ) s += values[i];

  return s * step;

}


function summation( f, [a,b] ) {

  if ( isComplex( f(a) ) ) {

    var s = complex(0);

    for ( var i = a ; i <= b ; i++ ) s = add( s, f(i) );

    return s;

  } else {

    var s = 0;

    for ( var i = a ; i <= b ; i++ ) s += f(i);

    return s;

  }

}

