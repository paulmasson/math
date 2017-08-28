
function ode( f, y, [x0, x1], step=.001, method='runge-kutta' ) {

  if ( f(x0,y)[0] === undefined ) {
    g = f;
    f = function(x) { return [ g(x) ]; };
    y = [ y ];
  }

  var points = [ [x0].concat(y) ];

  switch( method ) {

    case 'euler':

      for ( var x = x0+step ; x <= x1 ; x += step ) {

        var k = [];
        for ( var i = 0 ; i < y.length ; i++ ) k.push( f(x,y)[i] * step );

        for ( var i = 0 ; i < y.length ; i++ ) y[i] += k[i];
        points.push( [x].concat(y) );

      }

      return points;

    case 'runge-kutta':

      for ( var x = x0+step ; x <= x1 ; x += step ) {

        var k1 = [], k2 = [], k3 = [], k4 = [];
        var y1 = [], y2 = [], y3 = [];

        for ( var i = 0 ; i < y.length ; i++ ) k1.push( f(x,y)[i] );
        for ( var i = 0 ; i < y.length ; i++ ) y1.push( y[i] + k1[i]*step/2 );
        for ( var i = 0 ; i < y.length ; i++ ) k2.push( f( x+step/2, y1 )[i] );
        for ( var i = 0 ; i < y.length ; i++ ) y2.push( y[i] + k2[i]*step/2 );
        for ( var i = 0 ; i < y.length ; i++ ) k3.push( f( x+step/2, y2 )[i] );
        for ( var i = 0 ; i < y.length ; i++ ) y3.push( y[i] + k3[i]*step );
        for ( var i = 0 ; i < y.length ; i++ ) k4.push( f( x+step, y3 )[i] );

        for ( var i = 0 ; i < y.length ; i++ )
          y[i] += ( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] ) * step / 6;

        points.push( [x].concat(y) );

      }

      return points;

    default:

      throw 'Unsupported method';

  }

}




