
function hypergeometric0F1( a, x, tolerance=1e-10 ) {

  var useAsymptotic = 100;

  if ( isComplex(a) || isComplex(x) ) {

    if ( !isComplex(a) ) a = complex(a);
    if ( !isComplex(x) ) x = complex(x);

    if ( Number.isInteger(a.re) && a.re <= 0 && a.im === 0 )
      throw Error( 'Hypergeometric function pole' );

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

    while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {
      p = mul( p, div( div( x, a ), i ) );
      s = add( s, p );
      a = add( a, 1 );
      i++;
    }

    return s;

  } else {

    if ( Number.isInteger(a) && a <= 0 ) throw Error( 'Hypergeometric function pole' );

    // asymptotic form is complex
    if ( Math.abs(x) > useAsymptotic ) return hypergeometric0F1( a, complex(x) ).re;

    var s = 1;
    var p = 1;
    var i = 1;

    while ( Math.abs(p) > tolerance ) {
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

    if ( !isComplex(a) ) a = complex(a);
    if ( !isComplex(b) ) b = complex(b);
    if ( !isComplex(x) ) x = complex(x);

    if ( Number.isInteger(b.re) && b.re <= 0 && b.im === 0 )
      throw Error( 'Hypergeometric function pole' );

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

    while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {
      p = mul( p, div( div( mul( x, a ), b ), i ) );
      s = add( s, p );
      a = add( a, 1 );
      b = add( b, 1 );
      i++;
    }

    return s;

  } else {

    if ( Number.isInteger(b) && b <= 0 ) throw Error( 'Hypergeometric function pole' );

    // Kummer transformation
    if ( x < 0 ) return exp(x) * hypergeometric1F1( b-a, b, -x );

    // asymptotic form is complex
    if ( Math.abs(x) > useAsymptotic ) return hypergeometric1F1( a, b, complex(x) ).re;

    var s = 1;
    var p = 1;
    var i = 1;

    while ( Math.abs(p) > tolerance ) {
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

  var terms = 50;

  if ( isComplex(a) || isComplex(b) || isComplex(x) ) {

    if ( !isComplex(a) ) a = complex(a);
    if ( !isComplex(b) ) b = complex(b);
    if ( !isComplex(x) ) x = complex(x);

    var s = complex(1);
    var p = complex(1), pLast = p;
    var converging = false;
    var i = 1;

    while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {

      p = mul( p, div( mul( mul( x, a ), b ), i ) );

      if ( abs(p) > abs(pLast) && converging ) break; // prevent runaway sum
      if ( abs(p) < abs(pLast) ) converging = true;
      if ( i > terms ) throw Error( 'Not converging after ' + terms + ' terms' );

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
    var converging = false;
    var i = 1;

    while ( Math.abs(p) > tolerance ) {

      p *= x * a * b / i;

      if ( Math.abs(p) > Math.abs(pLast) && converging ) break; // prevent runaway sum
      if ( Math.abs(p) < Math.abs(pLast) ) converging = true;
      if ( i > terms ) throw Error( 'Not converging after ' + terms + ' terms' );

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

  if ( isComplex(a) || isComplex(b) || isComplex(c) || isComplex(x) ) {

    // choose smallest absolute value of transformed argument
    // transformations from dlmf.nist.gov/15.8

    var absArray = [ abs(x), abs(div(x,sub(x,1))), abs(inv(x)),
                     abs(inv(sub(1,x))), abs(sub(1,x)), abs(sub(1,inv(x))) ];

    var index = absArray.indexOf( Math.min.apply( null, absArray ) );

    switch( index ) {

      case 0:

        break;

      case 1:

        return mul( pow( sub(1,x), neg(a) ), hypergeometric2F1( a, sub(c,b), c, div(x,sub(x,1)) ) );

      case 2:

        var factor = div( sin( mul( pi, sub(b,a) ) ), mul( pi, gamma(c) ) );

        var t1 = mul( div( pow( neg(x), neg(a) ), mul( gamma(b), gamma(sub(c,a)), gamma(add(a,neg(b),1)) ) ),
                      hypergeometric2F1( a, add(a,neg(c),1), add(a,neg(b),1), inv(x) ) );

        var t2 = mul( div( pow( neg(x), neg(b) ), mul( gamma(a), gamma(sub(c,b)), gamma(add(b,neg(a),1)) ) ),
                      hypergeometric2F1( b, add(b,neg(c),1), add(b,neg(a),1), inv(x) ) );

        return div( sub( t1, t2 ), factor );

      case 3:

        var factor = div( sin( mul( pi, sub(b,a) ) ), mul( pi, gamma(c) ) );

        var t1 = mul( div( pow( sub(1,x), neg(a) ), mul( gamma(b), gamma(sub(c,a)), gamma(add(a,neg(b),1)) ) ),
                      hypergeometric2F1( a, sub(c,b), add(a,neg(b),1), inv(sub(1,x)) ) );

        var t2 = mul( div( pow( sub(1,x), neg(b) ), mul( gamma(a), gamma(sub(c,b)), gamma(add(b,neg(a),1)) ) ),
                      hypergeometric2F1( b, sub(c,a), add(b,neg(a),1), inv(sub(1,x)) ) );

        return div( sub( t1, t2 ), factor );

      case 4:

        var factor = div( sin( mul( pi, sub(c,add(a,b)) ) ), mul( pi, gamma(c) ) );

        var t1 = mul( inv( mul( gamma(sub(c,a)), gamma(sub(c,b)), gamma(add(a,b,neg(c),1)) ) ),
                      hypergeometric2F1( a, b, add(a,b,neg(c),1), sub(1,x) ) );

        var t2 = mul( div( pow( sub(1,x), sub(c,add(a,b)) ),
                           mul( gamma(a), gamma(b), gamma(add(c,neg(a),neg(b),1)) ) ),
                      hypergeometric2F1( sub(c,a), sub(c,b), add(c,neg(a),neg(b),1), sub(1,x) ) );

        return div( sub( t1, t2 ), factor );

      case 5:

        var factor = div( sin( mul( pi, sub(c,add(a,b)) ) ), mul( pi, gamma(c) ) );

        var t1 = mul( div( pow( x, neg(a) ), mul( gamma(sub(c,a)), gamma(sub(c,b)), gamma(add(a,b,neg(c),1)) ) ),
                      hypergeometric2F1( a, add(a,neg(c),1), add(a,b,neg(c),1), sub(1,inv(x)) ) );

        var t2 = mul( div( mul( pow( sub(1,x), sub(c,add(a,b)) ), pow( x, sub(a,c) ) ),
                           mul( gamma(a), gamma(b), gamma(add(c,neg(a),neg(b),1)) ) ),
                      hypergeometric2F1( sub(c,a), sub(1,a), add(c,neg(a),neg(b),1), sub(1,inv(x)) ) );

        return div( sub( t1, t2 ), factor );

    }

    if ( !isComplex(a) ) a = complex(a);
    if ( !isComplex(b) ) b = complex(b);
    if ( !isComplex(c) ) c = complex(c);
    if ( !isComplex(x) ) x = complex(x);

    if ( Number.isInteger(c.re) && c.re <= 0 && c.im === 0 )
      throw Error( 'Hypergeometric function pole' );

    var s = complex(1);
    var p = complex(1);
    var i = 1;

    while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {
      p = mul( p, div( div( mul( mul( x, a ), b ), c ), i ) );
      s = add( s, p );
      a = add( a, 1 );
      b = add( b, 1 );
      c = add( c, 1 );
      i++;
    }

    return s;

  } else {

    if ( x > 1 || x < 0 ) throw Error( 'Unsupported real hypergeometric argument' );

    if ( Number.isInteger(c) && c <= 0 ) throw Error( 'Hypergeometric function pole' );

    if ( x === 1 ) return gamma(c) * gamma(c-a-b) / gamma(c-a) / gamma(c-b);

    var s = 1;
    var p = 1;
    var i = 1;

    while ( Math.abs(p) > tolerance ) {
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

