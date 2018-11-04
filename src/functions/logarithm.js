
function exp( x ) {

  if ( isComplex(x) )

    return { re: Math.exp(x.re) * Math.cos(x.im),
             im: Math.exp(x.re) * Math.sin(x.im) };

  return Math.exp(x);

}


function log( x, base ) {

  if ( isComplex(x) ) {

    if ( isComplex(base) ) return div( log(x), log(base) );

    return { re: log( abs(x), base ), im: log( Math.E, base ) * arg(x) };

  }

  if ( x < 0 ) return log( complex(x), base );

  if ( base === undefined ) return Math.log(x);

  return Math.log(x) / Math.log(base);

}

var ln = log;


function lambertW( k, x ) {

  if ( arguments.length === 1 ) {
    x = k;
    k = 0;
  }

  if ( Math.abs( x + Math.exp(-1) ) < 1e-16 ) return -1;

  // inversion by root finding

  switch ( k ) {

    case 0:

      if ( x < -Math.exp(-1) ) throw 'Unsupported lambertW argument';

      return findRoot( w => w * Math.exp(w) - x, [-1,1000], { tolerance: 1e-16 } );

    case -1:

      if ( x < -Math.exp(-1) || x > 0 ) throw 'Unsupported lambertW argument';

      return findRoot( w => w * Math.exp(w) - x, [-1000,-1], { tolerance: 1e-16 } );

    default:

      throw 'Unsupported lambertW index';

  }

}

