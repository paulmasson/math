
function hypergeometric0F1( a, x, tolerance=1e-10 ) {

  var useAsymptotic = 100;

  if ( isComplex(a) || isComplex(x) ) {

    if ( !isComplex(a) ) a = complex(a);
    if ( !isComplex(x) ) x = complex(x);

    if ( Number.isInteger(a.re) && a.re <= 0 && a.im === 0 )
      throw Error( 'Hypergeometric function pole' );

    // asymptotic form as per Johansson arxiv.org/abs/1606.06977
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

    // asymptotic form as per Johansson arxiv.org/abs/1606.06977
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

    var s = complex(1);
    var p = complex(1), pLast = p;
    var converging = false;
    var i = 1;

    while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {

      p = mul( p, x, a, b, 1/i );

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
    // transformations from Abramowitz & Stegun p.559
    // fewer operations compared to dlmf.nist.gov/15.8

    var absArray = [ abs(x), abs(div(x,sub(x,1))), abs(sub(1,x)),
                     abs(inv(x)), abs(inv(sub(1,x))), abs(sub(1,inv(x))) ];

    var index = absArray.indexOf( Math.min.apply( null, absArray ) );

    switch( index ) {

      case 0:

        break;

      case 1:

        return mul( pow( sub(1,x), neg(a) ), hypergeometric2F1( a, sub(c,b), c, div(x,sub(x,1)) ) );

      case 2:

        var t1 = mul( gamma(c), gamma( sub(c,add(a,b)) ), 
                      inv( gamma(sub(c,a)) ), inv( gamma(sub(c,b)) ),
                      hypergeometric2F1( a, b, add(a,b,neg(c),1), sub(1,x) ) );

        var t2 = mul( pow( sub(1,x), sub(c,add(a,b)) ),
                      gamma(c), gamma( sub(add(a,b),c) ), inv( gamma(a) ), inv( gamma(b) ),
                      hypergeometric2F1( sub(c,a), sub(c,b), add(c,neg(a),neg(b),1), sub(1,x) ) );

        return add( t1, t2 );

      case 3:

        var t1 = mul( gamma(c), gamma(sub(b,a)), inv( gamma(b) ),
                      inv( gamma(sub(c,a)) ), pow( neg(x), neg(a) ),
                      hypergeometric2F1( a, add(1,neg(c),a), add(1,neg(b),a), inv(x) ) );

        var t2 = mul( gamma(c), gamma(sub(a,b)), inv( gamma(a) ),
                      inv( gamma(sub(c,b)) ), pow( neg(x), neg(b) ),
                      hypergeometric2F1( b, add(1,neg(c),b), add(1,neg(a),b), inv(x) ) );

        return add( t1, t2 );

      case 4:

        var t1 = mul( pow( sub(1,x), neg(a) ), gamma(c), gamma(sub(b,a)),
                      inv( gamma(b) ), inv( gamma(sub(c,a)) ),
                      hypergeometric2F1( a, sub(c,b), add(a,neg(b),1), inv(sub(1,x)) ) );

        var t2 = mul( pow( sub(1,x), neg(b) ), gamma(c), gamma(sub(a,b)),
                      inv( gamma(a) ), inv( gamma(sub(c,b)) ),
                      hypergeometric2F1( b, sub(c,a), add(b,neg(a),1), inv(sub(1,x)) ) );

        return add( t1, t2 );

      case 5:

        var t1 = mul( gamma(c), gamma( sub(c,add(a,b)) ), inv( gamma(sub(c,a)) ),
                      inv( gamma(sub(c,b)) ), pow( x, neg(a) ),
                      hypergeometric2F1( a, add(a,neg(c),1), add(a,b,neg(c),1), sub(1,inv(x)) ) );

        var t2 = mul( gamma(c), gamma( sub(add(a,b),c) ), inv( gamma(a) ), inv( gamma(b) ),
                      pow( sub(1,x), sub(c,add(a,b)) ), pow( x, sub(a,c) ),
                      hypergeometric2F1( sub(c,a), sub(1,a), add(c,neg(a),neg(b),1), sub(1,inv(x)) ) );

        return add( t1, t2 );

    }

    if ( !isComplex(c) ) c = complex(c);

    if ( Number.isInteger(c.re) && c.re <= 0 && c.im === 0 )
      throw Error( 'Hypergeometric function pole' );

    var s = complex(1);
    var p = complex(1);
    var i = 1;

    while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {
      p = mul( p, x, a, b, inv(c), 1/i );
      s = add( s, p );
      a = add( a, 1 );
      b = add( b, 1 );
      c = add( c, 1 );
      i++;
    }

    return s;

  } else {

    if ( Number.isInteger(c) && c <= 0 ) throw Error( 'Hypergeometric function pole' );

    // transformation from Abramowitz & Stegun p.559
    if ( x < -1 ) {

      var t1 = gamma(c) * gamma(b-a) / gamma(b) / gamma(c-a)
                 * (-x)**(-a) * hypergeometric2F1( a, 1-c+a, 1-b+a, 1/x );
      var t2 = gamma(c) * gamma(a-b) / gamma(a) / gamma(c-b)
                 * (-x)**(-b) * hypergeometric2F1( b, 1-c+b, 1-a+b, 1/x );

      return t1 + t2;

    }

    if ( x === -1 ) throw Error( 'Unsupported real hypergeometric argument' );

    if ( x === 1 ) return gamma(c) * gamma(c-a-b) / gamma(c-a) / gamma(c-b);

    if ( x > 1 ) return hypergeometric2F1( a, b, c, complex(x) );

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


function hypergeometricPFQ( A, B, x, tolerance=1e-10 ) {

  var complexArg = false;
  A.forEach( a => complexArg = complexArg || isComplex(a) );
  B.forEach( b => complexArg = complexArg || isComplex(b) );

  if ( complexArg || isComplex(x) ) {

    var s = complex(1);
    var p = complex(1);
    var i = 1;

    while ( Math.abs(p.re) > tolerance || Math.abs(p.im) > tolerance ) {
      p = mul( p, x, A.reduce( (x,y) => mul(x,y) ), inv( B.reduce( (x,y) => mul(x,y) ) ), 1/i );
      s = add( s, p );
      A.forEach( (e,i,a) => a[i] = add( a[i], 1 ) );
      B.forEach( (e,i,a) => a[i] = add( a[i], 1 ) );
      i++;
    }

    return s;

  } else {

    var s = 1;
    var p = 1;
    var i = 1;

    while ( Math.abs(p) > tolerance ) {
      p *= x * A.reduce( (x,y) => x*y ) / B.reduce( (x,y) => x*y ) / i;
      s += p;
      A.forEach( (e,i,a) => a[i]++ );
      B.forEach( (e,i,a) => a[i]++ );
      i++;
    }

    return s;

  }

}

