
var pi = Math.PI;

var eulerGamma = .5772156649015329;

var constants = {

  decimals: 50, precisionScale: 10n**50n,

  e: 271828182845904523536028747135266249775724709369995n,

  eulerGamma: 57721566490153286060651209008240243104215933593992n,

  pi: 314159265358979323846264338327950288419716939937510n

};

function getConstant( name ) {

  return constants[name] * precisionScale / constants.precisionScale;

}


var factorialCache = [ 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800 ];

