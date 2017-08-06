
function complex( x, y ) {

  var y = y || 0;
  return { re: x, im: y };

}

var C = complex;

// arrays of length greater than two are true for this test
// need a fast way to distinguish between array and dictionary

var isComplex = isNaN;


// JavaScript does not yet support operator overloading

function add( x, y ) {

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x,0);
    if ( !isComplex(y) ) y = complex(y,0);

    return { re: x.re + y.re, im: x.im + y.im };

  } else return x + y;

}

function sub( x, y ) {

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x,0);
    if ( !isComplex(y) ) y = complex(y,0);

    return { re: x.re - y.re, im: x.im - y.im };

  } else return x - y;

}

function mul( x, y ) {

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x,0);
    if ( !isComplex(y) ) y = complex(y,0);

    return { re: x.re * y.re - x.im * y.im,
             im: x.im * y.re + x.re * y.im };

  } else return x * y;

}

function div( x, y ) {

  // need to handle 0/0...

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x,0);
    if ( !isComplex(y) ) y = complex(y,0);

    var ySq = y.re * y.re + y.im * y.im;
    return { re: ( x.re * y.re + x.im * y.im ) / ySq,
             im: ( x.im * y.re - x.re * y.im ) / ySq };

  } else return x / y;

}

function pow( x, y ) {

  if ( isComplex(x) || isComplex(y) ) {

    if ( !isComplex(x) ) x = complex(x,0);
    if ( !isComplex(y) ) y = complex(y,0);

    var r = Math.sqrt( x.re * x.re + x.im * x.im );
    var phi = Math.atan2( x.im, x.re );

    var R = r**y.re * Math.exp( -phi * y.im );
    var Phi = phi * y.re + y.im * Math.log(r); // NaN at origin...

    return { re: R * Math.cos(Phi), im: R * Math.sin(Phi) };

  } else return x**y;


}

function root( x, y ) { return pow( x, div( 1, y ) ); }

function sqrt( x ) {

  if ( isComplex(x) ) {

    var R = ( x.re * x.re + x.im * x.im )**(1/4);
    var Phi = Math.atan2( x.im, x.re ) / 2;

    return { re: R * Math.cos(Phi), im: R * Math.sin(Phi) };

  } else return Math.sqrt(x);

}


