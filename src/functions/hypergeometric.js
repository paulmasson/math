
function hypergeometric0F1( a, x ) {

  if ( isComplex(a) || isComplex(x) ) {

    if ( Number.isInteger(a.re) && a.re <= 0 && a.im === 0 )
      throw 'Hypergeometric function pole';

    if ( !isComplex(a) ) a = complex(a,0);
    if ( !isComplex(x) ) x = complex(x,0);

    var s = complex(1);
    var p = complex(1);

    for ( var i = 1 ; i < 100 ; i++ ) {
      p = mul( p, div( div( x, a ), i ) );
      s = add( s, p );
      a = add( a, 1 );
    }

    return s;

  } else {

    if ( Number.isInteger(a) && a <= 0 ) throw 'Hypergeometric function pole';

    var s = 1;
    var p = 1;

    for ( var i = 1 ; i < 100 ; i++ ) {
      p *= x / a / i;
      s += p;
      a++;
    }

    return s;

  }

}


function hypergeometric1F1( a, b, x ) {

  if ( isComplex(a) || isComplex(b) || isComplex(x) ) {

    if ( Number.isInteger(b.re) && b.re <= 0 && b.im === 0 )
      throw 'Hypergeometric function pole';

    if ( !isComplex(a) ) a = complex(a,0);
    if ( !isComplex(b) ) b = complex(b,0);
    if ( !isComplex(x) ) x = complex(x,0);

    var s = complex(1);
    var p = complex(1);

    for ( var i = 1 ; i < 100 ; i++ ) {
      p = mul( p, div( div( mul( x, a ), b ), i ) );
      s = add( s, p );
      a = add( a, 1 );
      b = add( b, 1 );
    }

    return s;

  } else {

    if ( Number.isInteger(b) && b <= 0 ) throw 'Hypergeometric function pole';

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

}


function hypergeometric2F1( a, b, c, x ) {

  if ( isComplex(a) || isComplex(b) || isComplex(c) || isComplex(x) ) {

    if ( Number.isInteger(c.re) && c.re <= 0 && c.im === 0 )
      throw 'Hypergeometric function pole';

    if ( !isComplex(a) ) a = complex(a,0);
    if ( !isComplex(b) ) b = complex(b,0);
    if ( !isComplex(c) ) c = complex(c,0);
    if ( !isComplex(x) ) x = complex(x,0);

    var s = complex(1);
    var p = complex(1);

    for ( var i = 1 ; i < 1000 ; i++ ) {
      p = mul( p, div( div( mul( mul( x, a ), b ), c ), i ) );
      s = add( s, p );
      a = add( a, 1 );
      b = add( b, 1 );
      c = add( c, 1 );
    }

    return s;

  } else {

    if ( Number.isInteger(c) && c <= 0 ) throw 'Hypergeometric function pole';

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

}

