
function fourierSineCoefficient( f, n, period ) {

  if ( !Number.isInteger(n) ) throw 'Nonintegral Fourier index';

  if ( typeof f === 'function' ) {

    var T = period || 2*pi;

    return 2/T * integrate( t => f(t) * sin( 2*n*pi/T * t ), [0,T], 'tanh-sinh' );

  }

  if ( Array.isArray(f) ) {

    var s = 0, N = f.length;

    for ( var i = 0 ; i < N ; i++ ) s += f[i][1] * sin( 2*n*pi*i/N );

    return 2 * s / N;

  }

  throw 'Unsupported Fourier input';

}

function fourierCosineCoefficient( f, n, period ) {

  if ( !Number.isInteger(n) ) throw 'Nonintegral Fourier index';

  if ( typeof f === 'function' ) {

    var T = period || 2*pi;

    if ( n === 0 ) return 1/T * integrate( t => f(t), [0,T], 'tanh-sinh' );

    return 2/T * integrate( t => f(t) * cos( 2*n*pi/T * t ), [0,T], 'tanh-sinh' );

  }

  if ( Array.isArray(f) ) {

    var s = 0, N = f.length;

    if ( n === 0 ) {

      for ( var i = 0 ; i < N ; i++ ) s += f[i][1];

      return s / N;

    }

    for ( var i = 0 ; i < N ; i++ ) s += f[i][1] * cos( 2*n*pi*i/N );

    return 2 * s / N;

  }

  throw 'Unsupported Fourier input';

}

