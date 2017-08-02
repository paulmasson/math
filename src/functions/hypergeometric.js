
function hypergeometric0F1( a, x ) {

  if ( Number.isInteger(a) && a < 0 ) throw( 'Hypergeometric function pole' );

  var s = 1;
  var p = 1;

  for ( var i = 1 ; i < 100 ; i++ ) {
    p *= x / a / i;
    s += p;
    a++;
  }

  return s;

}


function hypergeometric1F1( a, b, x ) {

  if ( Number.isInteger(b) && b < 0 ) throw( 'Hypergeometric function pole' );

  var s = 1;
  var p = 1;

  for ( var i = 1 ; i < 100 ; i++ ) {
    p *= x * a / b / i;
    s += p;
    a++;
    b++;
  }

  return s;

}


function hypergeometric2F1( a, b, c, x ) {

  if ( Number.isInteger(c) && c < 0 ) throw( 'Hypergeometric function pole' );

  if ( x === 1 ) return gamma(c) * gamma(c-a-b) / gamma(c-a) / gamma(c-b);

  var s = 1;
  var p = 1;

  for ( var i = 1 ; i < 1000 ; i++ ) {
    p *= x * a * b / c / i;
    s += p;
    a++;
    b++;
    c++;
  }

  return s;

}

