
function vector( size, value=0 ) {

  var v = [];
  for ( var i = 0 ; i < size ; i++ ) v.push( value );

  return v;

}

function matrix( rows, columns, value=0 ) {

  var columns = columns || rows;

  var m = [];
  for ( var i = 0 ; i < rows ; i++ ) {
    m.push( [] );
    for ( var j = 0 ; j < columns ; j++ ) m[i].push( value );
  }

  return m;

}

function identity( rows, value=1 ) {

  var m = matrix( rows );
  for ( var i = 0 ; i < rows ; i++ ) m[i][i] = value;

  return m;

}

function transpose( A ) {

  var T = matrix( A[0].length, A.length );

  for ( var i = 0 ; i < A.length ; i++ )
    for ( var j = 0 ; j < A[0].length ; j++ )
      T[j][i] = A[i][j];

  return T;

}

function matrixAdd( A, B ) {

  if ( !Array.isArray(A) && !Array.isArray(B) ) throw 'No matrices';
  if ( !Array.isArray(A) ) A = matrix( B.length, B[0].length, A );
  if ( !Array.isArray(B) ) B = matrix( A.length, A[0].length, B );

  var C = matrix( A.length, A[0].length, 0 );

  for ( var i = 0 ; i < A.length ; i++ )
    for ( var j = 0 ; j < A[0].length ; j++ )
      C[i][j] = add( A[i][j], B[i][j] );

  return C;

}

function matrixSub( A, B ) {

  if ( !Array.isArray(A) && !Array.isArray(B) ) throw 'No matrices';
  if ( !Array.isArray(A) ) A = matrix( B.length, B[0].length, A );
  if ( !Array.isArray(B) ) B = matrix( A.length, A[0].length, B );

  var C = matrix( A.length, A[0].length, 0 );

  for ( var i = 0 ; i < A.length ; i++ )
    for ( var j = 0 ; j < A[0].length ; j++ )
      C[i][j] = sub( A[i][j], B[i][j] );

  return C;

}

function matrixMul( A, B ) {

  if ( !Array.isArray(A) && !Array.isArray(B) ) throw 'No matrices';
  if ( !Array.isArray(A) ) A = identity( B.length, A );
  if ( !Array.isArray(B) ) B = identity( A[0].length, B );
  if ( A[0].length !== B.length ) throw( 'Incompatible matrices' );

  var C = matrix( A.length, B[0].length, 0 );

  for ( var i = 0 ; i < A.length ; i++ )
    for ( var j = 0 ; j < B[0].length ; j++ )
      for ( var k = 0 ; k < A[0].length ; k++ )
        C[i][j] = add( C[i][j], mul( A[i][k], B[k][j] ) );

  return C;

}

