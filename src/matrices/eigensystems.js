
function eigensystem( A, symmetric=true ) {

  if ( symmetric ) return tridiagonalQL( tridiagonalForm(A) );
  else throw( 'Unsupported system' );

}

// sourced from http://math.nist.gov/javanumerics/jama/
// no need to reinvent this wheel...

function tridiagonalForm( A ) {

  var n = A.length;
  var V = [];
  for ( var i = 0 ; i < n ; i++ ) V[i] = A[i].slice(); // deeper copy

  var d = vector( n );
  var e = vector( n );

  for ( var j = 0 ; j < n ; j++ ) d[j] = V[n-1][j];

  // Householder reduction to tridiagonal form
   
  for ( var i = n - 1 ; i > 0 ; i-- ) {
   
    // scale to avoid under/overflow
   
    var scale = 0;
    var h = 0;
    for ( var k = 0 ; k < i ; k++ ) scale += Math.abs(d[k]);

    if ( scale === 0 ) {
      e[i] = d[i-1];
      for ( var j = 0 ; j < i ; j++ ) {
        d[j] = V[i-1][j];
        V[i][j] = 0;
        V[j][i] = 0;
      }
    } else {

      // generate Householder vector

      for ( var k = 0 ; k < i ; k++ ) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      var f = d[i-1];
      var g = Math.sqrt(h);
      if ( f > 0 ) g = -g;
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for ( var j = 0; j < i; j++ ) e[j] = 0;

      // apply similarity transformation to remaining columns

      for ( var j = 0 ; j < i ; j++ ) {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for ( var k = j + 1 ; k <= i - 1 ; k++ ) {
          g += V[k][j] * d[k];
          e[k] += V[k][j] * f;
        }
        e[j] = g;
      }
      f = 0;
      for ( var j = 0 ; j < i ; j++ ) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      var hh = f / ( h + h );
      for ( var j = 0 ; j < i ; j++ ) {
        e[j] -= hh * d[j];
      }
      for ( var j = 0 ; j < i ; j++ ) {
        f = d[j];
        g = e[j];
        for ( var k = j ; k <= i - 1 ; k++ ) {
          V[k][j] -= f * e[k] + g * d[k];
        }
        d[j] = V[i-1][j];
        V[i][j] = 0;
      }

    }

    d[i] = h;

  }
   
  // accumulate transformations
   
  for ( var i = 0 ; i < n - 1 ; i++ ) {
    V[n-1][i] = V[i][i];
    V[i][i] = 1;
    var h = d[i+1];
    if ( h !== 0 ) {
      for ( var k = 0 ; k <= i ; k++ ) d[k] = V[k][i+1] / h;
      for ( var j = 0 ; j <= i ; j++ ) {
        var g = 0;
        for ( var k = 0 ; k <= i ; k++ ) g += V[k][i+1] * V[k][j];
        for ( var k = 0 ; k <= i ; k++ ) V[k][j] -= g * d[k];
      }
    }
    for ( var k = 0; k <= i; k++) V[k][i+1] = 0;
  }
  for ( var j = 0 ; j < n ; j++ ) {
    d[j] = V[n-1][j];
    V[n-1][j] = 0;
  }
  V[n-1][n-1] = 1;
  e[0] = 0;

  return { diagonal:d, offDiagonal:e, eigenvectors:V };

}


function tridiagonalQL( tridiagonalForm ) {

  var d = tridiagonalForm.diagonal;
  var n = d.length;
  var e = tridiagonalForm.offDiagonal;
  var V = tridiagonalForm.eigenvectors;

  function hypot( a, b) {
    var r;
    if ( Math.abs(a) > Math.abs(b) ) {
      r = b/a;
      r = Math.abs(a) * Math.sqrt( 1 + r*r );
    } else if (b != 0) {
      r = a/b;
      r = Math.abs(b) * Math.sqrt( 1 + r*r );
    } else r = 0;

    return r;

  }

  for ( var i = 1 ; i < n ; i++ ) e[i-1] = e[i];
  e[n-1] = 0;

  var f = 0;
  var tst1 = 0;
  var eps = Math.pow( 2, -52 );

  for ( var l = 0 ; l < n ; l++ ) {

    // find small subdiagonal element

    tst1 = Math.max( tst1, Math.abs(d[l]) + Math.abs(e[l]) );
    var m = l;
    while ( m < n ) {
      if ( Math.abs(e[m]) <= eps*tst1 ) break;
      m++;
    }

    // if m === l, d[l] is an eigenvalue, otherwise iterate

    if ( m > l ) {

      var iter = 0;
      do {

        iter = iter + 1;
        if ( iter > 1000 ) throw( 'Eigenvalues not converging...' );

        // compute implicit shift

        var g = d[l];
        var p = ( d[l+1] - g ) / ( 2 * e[l] );
        var r = hypot(p,1);
        if ( p < 0 ) r = -r;
        d[l] = e[l] / ( p + r );
        d[l+1] = e[l] * ( p + r );
        var dl1 = d[l+1];
        var h = g - d[l];
        for ( var i = l + 2 ; i < n ; i++ ) d[i] -= h;
        f = f + h;

        // implicit QL transformation

        p = d[m];
        var c = 1;
        var c2 = c;
        var c3 = c;
        var el1 = e[l+1];
        var s = 0;
        var s2 = 0;
        for ( var i = m - 1 ; i >= l ; i-- ) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = hypot(p,e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * ( c * g + s * d[i] );

          // accumulate transformation

          for ( var k = 0 ; k < n ; k++ ) {
            h = V[k][i+1];
            V[k][i+1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
   
        // check for convergence
   
      } while ( Math.abs(e[l]) > eps*tst1 );
    }
    d[l] = d[l] + f;
    e[l] = 0;
  }
   
  // sort eigenvalues and corresponding vectors

  for ( var i = 0 ; i < n - 1 ; i++ ) {
    var k = i;
    var p = d[i];
    for ( var j = i + 1 ; j < n ; j++ ) {
      if ( d[j] < p ) {
        k = j;
        p = d[j];
      }
    }
    if ( k != i ) {
      d[k] = d[i];
      d[i] = p;
      for ( var j = 0 ; j < n ; j++ ) {
        p = V[j][i];
        V[j][i] = V[j][k];
        V[j][k] = p;
      }
    }
  }

  return { eigenvalues:d, eigenvectors:V };

}


function hessenbergForm( A ) {

}

