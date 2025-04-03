// Compute eccentric anomaly from mean anomaly and eccentricity
class Kepler {
  constructor(e) {
    this.maxIter = 10;
    this.tol = 1e-4;
    this.e = e;
  }
  solve(M) {
    var E = (M<Math.PI ? M + this.e/2 : M - this.e/2);
    var dE = 1000;
    for (let i=0; i<this.maxIter && Math.abs(dE)>this.tol; ++i) {
      var f = E - this.e*Math.sin(E) - M;
      var df = 1 - this.e*Math.cos(E);
      dE = f/df;
      E -= dE;
    }
    if (Math.abs(dE) > tol) console.log(E, dE);
    return E;
  }
};

