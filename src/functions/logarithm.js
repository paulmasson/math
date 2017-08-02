
function exp( x ) {

  if ( isComplex(x) )

    return { re: Math.exp(x.re) * Math.cos(x.im),
             im: Math.exp(x.re) * Math.sin(x.im) };

  else return Math.exp(x);

}


function log( x, base ) {

  if ( isComplex(x) ) {

    var r = sqrt( x.re * x.re + x.im * x.im );
    var phi = Math.atan2( x.im, x.re );

    return { re: log(r,base), im: log(e,base) * phi };

  } else if ( x < 0 ) return log( complex(x), base );

  else if ( base === undefined ) return Math.log(x);

  else return Math.log(x) / Math.log(base);

}

var ln = log;

