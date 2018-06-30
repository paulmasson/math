
function exp( x ) {

  if ( isComplex(x) )

    return { re: Math.exp(x.re) * Math.cos(x.im),
             im: Math.exp(x.re) * Math.sin(x.im) };

  return Math.exp(x);

}


function log( x, base ) {

  if ( isComplex(x) ) {

    var r = sqrt( x.re * x.re + x.im * x.im );
    var phi = Math.atan2( x.im, x.re );

    return { re: log(r,base), im: log(Math.E,base) * phi };

  }

  if ( x < 0 ) return log( complex(x), base );

  if ( base === undefined ) return Math.log(x);

  return Math.log(x) / Math.log(base);

}

var ln = log;


function lambertW( x ) {

  // simple inversion by root finding

  return findRoot( w => w * Math.exp(w) - x, [-1,100], { tolerance: 1e-16 } );

}

