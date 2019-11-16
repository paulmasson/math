
function ode( f, y, [x0,x1], step=.001, method='runge-kutta' ) {

  if ( x1 < x0 ) {
    function compare( x ) { return x >= x1; };
    step *= -1;
  } else
    function compare( x ) { return x <= x1; };

  if ( f(x0,y)[0] === undefined ) {
    g = f;
    f = function(x,y) { return [ g(x,y) ]; };
    y = [ y ];
  }

  var points = [ [x0].concat(y) ];
  var size = y.length;

  switch( method ) {

    case 'euler':

      for ( var x = x0+step ; compare(x) ; x += step ) {

        var k = f(x,y);

        for ( var i = 0 ; i < size ; i++ ) y[i] += k[i] * step;

        points.push( [x].concat(y) );

      }

      return points;

    case 'runge-kutta':

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

    default:

      throw Error( 'Unsupported differential equation solver method' );

  }

}




