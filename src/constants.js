
var constants = {

  decimals: 50, precisionScale: 10n**50n,

  e: 271828182845904523536028747135266249775724709369995n,

  eulerGamma: 57721566490153286060651209008240243104215933593992n,

  pi: 314159265358979323846264338327950288419716939937510n

};

var pi = Number( constants.pi ) / 10**constants.decimals;

var eulerGamma = Number( constants.eulerGamma ) / 10**constants.decimals;


var factorialCache = [ 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800 ];

