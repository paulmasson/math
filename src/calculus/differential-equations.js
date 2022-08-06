
function ode( f, y, [x0,x1], step=.001, method='runge-kutta' ) {

  if ( x1 < x0 ) {
    var compare = x => x >= x1;
    step *= -1;
  } else
    var compare = x => x <= x1;

  // vectorizing first-order real equation works because +[1] = 1
  // for complex case +[C(1)] = NaN, so explicit array references
  //    are necessary in the input function

  if ( !Array.isArray(y) ) {
    var g = f;
    f = (x,y) => [ g(x,y) ];
    y = [ y ];
  }

  // preparation for complex system
  if ( isComplex(x0) || isComplex(x1) || y.some( e => isComplex(e) )
         || f(x0,y).some( e => isComplex(e) ) ) {

    if ( !isComplex(x0) ) x0 = complex(x0);

    y.forEach( (e,i,a) => { if ( !isComplex(e) ) a[i] = complex(e); } );

    if ( f(x0,y).every( e => !isComplex(e) ) )
      throw Error( 'All functions must handle complex math' );

    var d = sub(x1,x0), absD = abs(d);
    step = mul( step, div( d, absD ) );
    var steps = Math.trunc( absD / abs(step) ), currentStep = 0;

  }

  var points = [ [x0].concat(y) ];
  var size = y.length;

  switch( method ) {

    case 'euler':

      if ( isComplex(x0) ) {

        for ( var x = add(x0,step) ; currentStep < steps ; x = add(x,step) ) {

          var k = f(x,y);

          for ( var i = 0 ; i < size ; i++ ) y[i] = add( y[i], mul( k[i], step ) );

          points.push( [x].concat(y) );

          currentStep++;

        }

        return points;

      } else {

        for ( var x = x0+step ; compare(x) ; x += step ) {

          var k = f(x,y);

          for ( var i = 0 ; i < size ; i++ ) y[i] += k[i] * step;

          points.push( [x].concat(y) );

        }

        return points;

      }

    case 'runge-kutta':

      if ( isComplex(x0) ) {

        var halfStep = div( step, 2 );

        for ( var x = add(x0,step) ; currentStep < steps ; x = add(x,step) ) {

          var y1 = [], y2 = [], y3 = [];

          var k1 = f(x,y);
          for ( var i = 0 ; i < size ; i++ ) y1.push( add( y[i], mul( k1[i], halfStep ) ) );
          var k2 = f( add( x, halfStep ), y1 );
          for ( var i = 0 ; i < size ; i++ ) y2.push( add( y[i], mul( k2[i], halfStep ) ) );
          var k3 = f( add( x, halfStep ), y2 );
          for ( var i = 0 ; i < size ; i++ ) y3.push( add( y[i], mul( k3[i], step ) ) );
          var k4 = f( add( x, step ), y3 );

          for ( var i = 0 ; i < size ; i++ )
            y[i] = add( y[i], mul( add( k1[i], mul(2,k2[i]), mul(2,k3[i]), k4[i] ), step, 1/6 ) );

          points.push( [x].concat(y) );

          currentStep++;

        }

        return points;

      } else {

        for ( var x = x0+step ; compare(x) ; x += step ) {

          var y1 = [], y2 = [], y3 = [];

          var k1 = f(x,y);
          for ( var i = 0 ; i < size ; i++ ) y1.push( y[i] + k1[i]*step/2 );
          var k2 = f( x+step/2, y1 );
          for ( var i = 0 ; i < size ; i++ ) y2.push( y[i] + k2[i]*step/2 );
          var k3 = f( x+step/2, y2 );
          for ( var i = 0 ; i < size ; i++ ) y3.push( y[i] + k3[i]*step );
          var k4 = f( x+step, y3 );

          for ( var i = 0 ; i < size ; i++ )
            y[i] += ( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] ) * step / 6;

          points.push( [x].concat(y) );

        }

        return points;

      }

    default:

      throw Error( 'Unsupported differential equation solver method' );

  }

}

