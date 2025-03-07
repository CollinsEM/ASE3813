// quick 2D PGA test.
Algebra(2,0,1,()=>{
  // const params = {
  //   phi: 0
  // };
  // const GUI = lil.GUI;
  // const gui = new GUI();
  // gui.add( document, 'title' );
  // gui.add( params, 'phi' );

  // Using Geometric/Clifford algebra with signature (2,0,1), based on the work of Charles Gunn
  const point   = (x,y)=>!(1e0 + x*1e1 + y*1e2);
  const vector  = (x,y)=>!(x*1e1 + y*1e2);
  const dist    = (P,Q)=>((P.Normalized)&(Q.Normalized)).Length;
  const refl    = (l,m)=>(l.Normalized<<m.Normalized)*2*(m.Normalized)-l.Normalized;
  const project = (a,b)=>(a | b) / b;
  const reject  = (a,b)=>(a | b);
  const lerp    = (a,b,t)=>(1-t)*a + t*b;
  
  var toString = (v,n)=>v.toPrecision(n||1).toString();

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
  
  // Compute eccentric anomaly from mean anomaly and eccentricity
  var eccentricAnomaly = (M,e) => {
    const maxIter = 10;
    var E = (M<Math.PI ? M + e/2 : M - e/2);
    for (let i=0; i<maxIter; ++i) {
      var f = E - e*Math.sin(E) - M;
      var df = 1 - e*Math.cos(E);
      var de = f/df;
      if (Math.abs(de) < 0.0001) break;
      E -= de;
    }
    return E;
  }
  
  // Compute true anomaly from mean anomaly and eccentricity
  var trueAnomaly = (t,T,e) => {
    const M = 2*Math.PI*t/T;
    const E = eccentricAnomaly(M,e);
    return Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
  }
  
  // Compute true anomaly from mean anomaly and eccentricity
  var trueAnomalyFromMean = (M,e) => {
    const E = eccentricAnomaly(M,e);
    return Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
  }
  
  // Compute true anomaly from mean anomaly and eccentricity
  var trueAnomalyFromEccentric = (E,e) => {
    return Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
  }
  
  // Compute the location of a point on an ellipse
  // F0: primary focus location
  // a: semi-major axis length
  // e: eccentricity
  // t: true anomaly
  // p: semi-major axis
  // q: semi-parameter axis
  const Ellipse = (F0, a, e, t, p, q)=>{
    var r = a*(1-e*e)/(1+e*Math.cos(t));
    return F0 + r*Math.cos(t)*p + r*Math.sin(t)*q;
  }
  
  // Computes the location of a point on a circle of radius, r, centered at R0
  var Circle=(R0,r,theta)=>{
    return ()=>R0+vector(r*Math.cos(theta), r*Math.sin(theta));
  };

  // Closed curve
  var curve=(pts)=>pts.map((pt,i,arr)=>[pt, pts[i+1%arr.length]]);
  //var OpenCurve=(pts)=>pts.map((pt,i,arr)=>[(i?pts[i-1]:pt), pt]);

  var N=120, Nc=60;
  const mu = 1;
  
  // Random orbit parameters (for initialization)
  const rp  = 0.5 + 0.5*Math.random();
  const ra  = rp + 1.5*Math.random();
  const e0  = (ra - rp)/(ra + rp);
  const a0  = 0.5*(ra + rp);
  const h0  = Math.sqrt(mu*a0*(1-e0*e0));
  // Generate two random points on the above orbit (for initialization)
  const t1_ = 2*Math.PI*Math.random();
  const t2_ = 2*Math.PI*Math.random();
  const r1_ = a0*(1-e0*e0)/(1 + e0*Math.cos(t1_));
  const r2_ = a0*(1-e0*e0)/(1 + e0*Math.cos(t2_));
  
  // NOTE: The following points can be modified by dragging them
  // around the screen.
  
  // Initial location of the primary focus (draggable)
  // const F0  = O + vector(a0*e0, 0)
  const F1  = point(a0*e0, 0)
  // Initial locations of the two observations (draggable)
  const R1  = F1 + vector(r1_*Math.cos(t1_), r1_*Math.sin(t1_));
  const R2  = F1 + vector(r2_*Math.cos(t2_), r2_*Math.sin(t2_));

  // NOTE: All values inside of the following lambda will be
  // recomputed whenever any of the above points are modified by the
  // user.
  document.body.appendChild(this.graph(() => {
    // Invariants (for a given R1 and R2)
    const r1  = dist(F1, R1);
    const r2  = dist(F1, R2);
    const c   = dist(R1, R2);
    // Fundamental elliptic orbit
    const af  = (r1 + r2)/2;          // semi-major axis length
    const ef  = Math.abs(r1 - r2)/c;  // eccentricity
    const Pf  = (r1 < r2 ? (R1 - R2)/c : (R2 - R1)/c).Normalized; // perifocal unit vector (p)
    const O   = F1 - af*ef*Pf;        // ellipse origin location
    const F2  = O  - af*ef*Pf;        // secondary focus location
    //const O   = (F1 + F2).Normalized; // center of the orbital ellipse
    const Ef  = ef*Pf;                // eccentricity vector
    const Qf  = Pf*1e12;
    // Generate fundamental ellipse...
    // ... with uniform true anomaly
    // const FundEllipse = [...Array(N)].map((x,i)=>()=> {
    //   const th = trueAnomalyFromMean(2*Math.PI*i/N, ef);
    //   return Ellipse(F1, af, ef, th, Pf, Qf)
    // });
    // ... with uniform eccentric anomaly
    const FundEllipse = [...Array(N)].map((x,i)=>()=> {
      const th = trueAnomalyFromEccentric(2*Math.PI*i/N, ef);
      return Ellipse(F1, af, ef, th, Pf, Qf)
    });
    // ... with uniform mean anomaly
    // const FundEllipse = [...Array(N)].map((x,i)=>()=> {
    //   const th = trueAnomalyFromMean(2*Math.PI*i/N, ef);
    //   return Ellipse(F1, af, ef, th, Pf, Qf)
    // });
                                         
    const AuxCircle   = [...Array(Nc)].map((x,i)=>()=>Circle(O, af, 2*Math.PI*i/Nc));
   
    return [
      'r1: ' + r1.toFixed(3),
      'r2: ' + r2.toFixed(3),
      'c:  ' + c.toFixed(3),
      'a:  ' + af.toFixed(3),
      'e:  ' + ef.toFixed(3),
      0xAAAAAA,
      F1&F2,
      F1&F2|O,
      R1&R2,
      R1&R2|(F1+ef*Pf),
      0x888888,
      O,  'O',
      F2, 'F2',
      0x444444,
      F1, 'F1',
      0x44AA44,
      R1, 'R1',
      R2, 'R2',
      0x448844,
      [F1,R1], 'r1',
      [F1,R2], 'r2',
      [R1,R2], 'c',
      0xAAAAAA,
      [F2,R1], '2a-r1',
      [F2,R2], '2a-r2',
      0x888844,
      [F1,F1+0.25*Pf], 'p',
      [F1,F1+0.25*Qf], 'q',
      0xAA0000,
      ...curve(AuxCircle),
      [F1,(F1+ef*Pf)], 'e',
      0x00AA00,
      ...curve(FundEllipse),
    ]}, { grid: true,
          labels: false,
          lineWidth: 1,
          pointRadius: 0.5,
          fontSize: 0.5,
          scale: 0.5,
          width: window.innerWidth,
          height: window.innerHeight
        }));
});
