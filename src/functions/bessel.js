
function besselJ( n, x ) {

  return (x/2)**n * hypergeometric0F1( n+1, -.25*x**2 ) / gamma(n+1);

}

function besselY( n, x ) {

  return ( besselJ(n,x) * cos(n*pi) - besselJ(-n,x) ) / sin(n*pi);


}

function besselI( n, x ) {

  return (x/2)**n * hypergeometric0F1( n+1, .25*x**2 ) / gamma(n+1);

}

function besselK( n, x ) {

  return pi/2 * ( besselI(-n,x) - besselI(n,x) ) / sin(n*pi);


}


function airyAi( x ) {

  if ( x > 0 ) return 1/pi * sqrt(x/3) * besselK(1/3, 2/3*x**(3/2));

  if ( x === 0 ) return 1 / 3**(2/3) / gamma(2/3);

  if ( x < 0 ) return sqrt(-x) / 2 * (  besselJ(1/3, 2/3*(-x)**(3/2))
                                       - besselY(1/3, 2/3*(-x)**(3/2)) / sqrt(3) );

}

function airyBi( x ) {

  if ( x > 0 ) return sqrt(x/3) * ( besselI(1/3, 2/3*x**(3/2))
                                    + besselI(-1/3, 2/3*x**(3/2)) );

  if ( x === 0 ) return 1 / 3**(1/6) / gamma(2/3);

  if ( x < 0 ) return -sqrt(-x) / 2 * ( besselJ(1/3, 2/3*(-x)**(3/2)) / sqrt(3)
                                        + besselY(1/3, 2/3*(-x)**(3/2)) );

}
