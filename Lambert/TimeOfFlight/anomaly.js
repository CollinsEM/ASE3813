// Compute all anomaly types
class Anomaly {
  constructor(obj) {
    const maxIter=10;
    if (obj.M !== undefined && obj.e !== undefined) {
      const M = this.M = obj.M;
      const e = this.e = obj.e;
      var E = (M<Math.PI ? M + e/2 : M - e/2);
      for (let i=0; i<maxIter; ++i) {
        var f = E - e*Math.sin(E) - M;
        var df = 1 - e*Math.cos(E);
        var de = f/df;
        if (Math.abs(de) < 0.0001) {
          this.nitr = i;
          this.E = E;
          const theta = 2*Math.atan(Math.sqrt((1+e)/(1-e))*Math.tan(E/2));
          this.theta = Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
          console.log(theta, this.theta);
          break;
        }
        E -= de;
      }
    }
    else if (obj.E !== undefined && obj.e !== undefined) {
      const E = this.E = obj.E;
      const e = this.e = obj.e;
      this.M = E - e*Math.sin(E);
      const theta = 2*Math.atan(Math.sqrt((1+e)/(1-e))*Math.tan(E/2));
      this.theta = Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
      console.log(theta, this.theta);
    }
    else if (obj.theta !== undefined && obj.e !== undefined) {
      const theta = this.theta = obj.theta;
      const e = this.e = obj.e;
      const E = this.E = 2*Math.atan(Math.sqrt((1-e)/(1+e))*Math.tan(theta/2));
      this.M = E - e*Math.sin(E);
    }
  }
};

  // Compute true anomaly from time, period, and eccentricity
  var trueAnomaly = (t,T,ecc) => {
    const MA = 2*Math.PI*t/T;
    const EA = solveKeplersEq(MA,ecc);
    return Math.atan2(Math.sin(EA)*Math.sqrt(1-ecc*ecc),Math.cos(EA)-ecc);
  }
  
  // Compute true anomaly from mean anomaly and eccentricity
  var trueAnomalyFromMean = (MA,ecc) => {
    const EA = solveKeplersEq(MA,ecc);
    return Math.atan(Math.sqrt((1+ecc)/(1-ecc))*Math.tan(EA/2));
  }
  
  // Compute true anomaly from eccentric anomaly and eccentricity
  var trueAnomalyFromEccentric = (EA,ecc) => {
    return Math.atan(Math.sqrt((1+ecc)/(1-ecc))*Math.tan(EA/2));
  }
  
  // Compute eccentric anomaly from true anomaly and eccentricity
  var eccentricAnomalyFromTrue = (theta,ecc) => {
    return Math.atan(Math.sqrt((1-ecc)/(1+ecc))*Math.tan(theta/2));
  }
  
