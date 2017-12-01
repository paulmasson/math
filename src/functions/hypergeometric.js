
function hypergeometric0F1( a, x, tolerance=1e-10 ) {

  var useAsymptotic = 100;

  if ( isComplex(a) || isComplex(x) ) {

    if ( !isComplex(a) ) a = complex(a,0);
    if ( !isComplex(x) ) x = complex(x,0);

    if ( Number.isInteger(a.re) && a.re <= 0 && a.im === 0 )
      throw 'Hypergeometric function pole';

    // asymptotic form as per Johansson
    if ( abs(x) > useAsymptotic ) {

      var b = sub( mul(2,a), 1 ); // do first
      var a = sub( a, 1/2 );
      var x = mul( 4, sqrt(x) );

      // copied from hypergeometric1F1
      var t1 = div( mul( gamma(b), pow( mul(-1,x), mul(-1,a) ) ), gamma( sub(b,a) ) );
      t1 = mul( t1, hypergeometric2F0( a, add( sub(a,b), 1 ), div(-1,x) ) );

      var t2 = div( mul( gamma(b), mul( pow( x, sub(a,b) ), exp( x ) ) ), gamma(a) );
      t2 = mul( t2, hypergeometric2F0( sub(b,a), sub(1,a), div(1,x) ) );

      return mul( exp( div(x,-2) ), add( t1, t2 ) );

    }

    var s = complex(1);
    var p = complex(1);
    var i = 1;

    while ( Math.abs(p.re) > tolerance * Math.abs(s.re)
            || Math.abs(p.im) > tolerance * Math.abs(s.im) ) {
      p = mul( p, div( div( x, a ), i ) );
      s = add( s, p );
      a = add( a, 1 );
      i++;
    }

    return s;

  } else {

    if ( Number.isInteger(a) && a <= 0 ) throw 'Hypergeometric function pole';

    if ( Math.abs(x) > useAsymptotic ) return hypergeometric0F1( a, complex(x) ).re;

    var s = 1;
    var p = 1;
    var i = 1;

    while ( Math.abs(p) > tolerance * Math.abs(s) ) {
      p *= x / a / i;
      s += p;
      a++;
      i++;
    }

    return s;

  }

}


function hypergeometric1F1( a, b, x, tolerance=1e-10 ) {

  var useAsymptotic = 30;

  if ( isComplex(a) || isComplex(b) || isComplex(x) ) {

    if ( !isComplex(a) ) a = complex(a,0);
    if ( !isComplex(b) ) b = complex(b,0);
    if ( !isComplex(x) ) x = complex(x,0);

    if ( Number.isInteger(b.re) && b.re <= 0 && b.im === 0 )
      throw 'Hypergeometric function pole';

    // Kummer transformation
    if ( x.re < 0 ) return mul( exp(x), hypergeometric1F1( sub(b,a), b, mul(x,-1) ) );

    // asymptotic form as per Johansson
    if ( abs(x) > useAsymptotic ) {

      var t1 = div( mul( gamma(b), pow( mul(-1,x), mul(-1,a) ) ), gamma( sub(b,a) ) );
      t1 = mul( t1, hypergeometric2F0( a, add( sub(a,b), 1 ), div(-1,x) ) );

      var t2 = div( mul( gamma(b), mul( pow( x, sub(a,b) ), exp( x ) ) ), gamma(a) );
      t2 = mul( t2, hypergeometric2F0( sub(b,a), sub(1,a), div(1,x) ) );

      return add( t1, t2 );

    }

    var s = complex(1);
    var p = complex(1);
    var i = 1;

    while ( Math.abs(p.re) > tolerance * Math.abs(s.re)
            || Math.abs(p.im) > tolerance * Math.abs(s.im) ) {
      p = mul( p, div( div( mul( x, a ), b ), i ) );
      s = add( s, p );
      a = add( a, 1 );
      b = add( b, 1 );
      i++;
    }

    return s;

  } else {

    if ( Number.isInteger(b) && b <= 0 ) throw 'Hypergeometric function pole';

    // Kummer transformation
    if ( x < 0 ) return exp(x) * hypergeometric1F1( b-a, b, -x );

    if ( Math.abs(x) > useAsymptotic ) return hypergeometric1F1( a, b, complex(x) ).re;

    var s = 1;
    var p = 1;
    var i = 1;

    while ( Math.abs(p) > tolerance * Math.abs(s) ) {
      p *= x * a / b / i;
      s += p;
      a++;
      b++;
      i++;
    }

    return s;

  }

}


function hypergeometric2F0( a, b, x, tolerance=1e-10 ) {

  if ( isComplex(a) || isComplex(b) || isComplex(x) ) {

    if ( !isComplex(a) ) a = complex(a,0);
    if ( !isComplex(b) ) b = complex(b,0);
    if ( !isComplex(x) ) x = complex(x,0);

    var s = complex(1);
    var p = complex(1), pLast = p;
    var i = 1;

    while ( Math.abs(p.re) > tolerance * Math.abs(s.re)
            || Math.abs(p.im) > tolerance * Math.abs(s.im) ) {
      p = mul( p, div( mul( mul( x, a ), b ), i ) );
      if ( abs(p) > abs(pLast) ) break; // prevent runaway sum
      s = add( s, p );
      a = add( a, 1 );
      b = add( b, 1 );
      i++;
      pLast = p;
    }

    return s;

  } else {

    var s = 1;
    var p = 1, pLast = p;
    var i = 1;

    while ( Math.abs(p) > tolerance * Math.abs(s) ) {
      p *= x * a * b / i;
      if ( Math.abs(p) > Math.abs(pLast) ) break; // prevent runaway sum
      s += p;
      a++;
      b++;
      i++;
      pLast = p;
    }

    return s;

  }

}


function hypergeometric2F1( a, b, c, x, tolerance=1e-10 ) {

  if ( abs(x) > 1 ) throw 'Not yet implemented'

  if ( isComplex(a) || isComplex(b) || isComplex(c) || isComplex(x) ) {

    if ( !isComplex(a) ) a = complex(a,0);
    if ( !isComplex(b) ) b = complex(b,0);
    if ( !isComplex(c) ) c = complex(c,0);
    if ( !isComplex(x) ) x = complex(x,0);

    if ( Number.isInteger(c.re) && c.re <= 0 && c.im === 0 )
      throw 'Hypergeometric function pole';

    var s = complex(1);
    var p = complex(1);
    var i = 1;

    while ( Math.abs(p.re) > tolerance * Math.abs(s.re)
            || Math.abs(p.im) > tolerance * Math.abs(s.im) ) {
      p = mul( p, div( div( mul( mul( x, a ), b ), c ), i ) );
      s = add( s, p );
      a = add( a, 1 );
      b = add( b, 1 );
      c = add( c, 1 );
      i++;
    }

    return s;

  } else {

    if ( Number.isInteger(c) && c <= 0 ) throw 'Hypergeometric function pole';

    if ( x === 1 ) return gamma(c) * gamma(c-a-b) / gamma(c-a) / gamma(c-b);

    var s = 1;
    var p = 1;
    var i = 1;

    while ( Math.abs(p) > tolerance * Math.abs(s) ) {
      p *= x * a * b / c / i;
      s += p;
      a++;
      b++;
      c++;
      i++;
    }

    return s;

  }

}

