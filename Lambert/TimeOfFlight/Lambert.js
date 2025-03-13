function dot(a,b) {
  return a.e02*b.e02 + a.e01*b.e01;
}

const GUI = lil.GUI;

// Using Geometric/Clifford algebra with signature (2,0,1), based on the work of Charles Gunn
Algebra(2,0,1,()=>{

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
  const maxIter = 10;
  const tol = 1e-4;
  class Kepler {
    constructor(e) {
      this.e = e;
    }
    solve(M) {
      var E = (M<Math.PI ? M + this.e/2 : M - this.e/2);
      var dE = 1000;
      for (let i=0; i<maxIter && Math.abs(dE)>tol; ++i) {
        var f = E - this.e*Math.sin(E) - M;
        var df = 1 - this.e*Math.cos(E);
        dE = f/df;
        E -= dE;
      }
      if (Math.abs(dE) > tol) console.log(E, dE);
      return E;
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
  
  // Compute the location of a point on an ellipse
  // F0: primary focus location
  // a: semi-major axis length
  // e: eccentricity
  // t: true anomaly
  // p: semi-major axis
  // q: semi-parameter axis
  const Ellipse = (F0, a, e, t, P, Q)=>{
    var r = a*(1-e*e)/(1+e*Math.cos(t));
    return F0 + r*Math.cos(t)*P + r*Math.sin(t)*Q;
  }
  
  // Computes the location of a point on a circle of radius, r, centered at R0
  var Circle=(R0,r,theta)=>{
    return ()=>R0+vector(r*Math.cos(theta), r*Math.sin(theta));
  };

  // Closed curve
  var curve=(pts)=>pts.map((pt,i,arr)=>[pt, pts[i+1%arr.length]]);
  //var OpenCurve=(pts)=>pts.map((pt,i,arr)=>[(i?pts[i-1]:pt), pt]);

  var N=360, Nc=360;
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
  // const F0  = O + vector(a0*e0, 0);
  // const F1  = point(a0*e0, 0);
  const F1  = point(1, 0);
  // Initial locations of the two observations (draggable)
  // const R1  = F1 + vector(r1_*Math.cos(t1_), r1_*Math.sin(t1_));
  // const R2  = F1 + vector(r2_*Math.cos(t2_), r2_*Math.sin(t2_));
  const R1  = point( 1, 1);
  const R2  = point(-1, 1);
  // console.log(R1);
  // console.log(R2);
  // Rendered items
  var graph = [];
  
  class Controls extends GUI {
    constructor() {
      super();
      // Clear the graph
      graph = [];
      
      this.distUpdate();
      
      // Minimum eccentricity orbit
      this.eMin = this.addFolder('Minimum Eccentricity Orbit');
      this.eMin.show = true;
      this.eMin.add(this.eMin, 'show').listen()
        .onChange(()=>{ controls.eMinUpdate(); });
      this.eMin.a = 0;
      this.eMin.add(this.eMin, 'a').listen();
      this.eMin.e = 0;
      this.eMin.add(this.eMin, 'e').listen();
      this.eMinUpdate();
      
      // Minimum energy orbit
      this.aMin = this.addFolder('Minimum Energy Orbit');
      this.aMin.show = true;
      this.aMin.add(this.aMin, 'show').listen()
        .onChange(()=>{ controls.aMinUpdate(); });
      this.aMin.a = 0;
      this.aMin.add(this.aMin, 'a').listen();
      this.aMin.e = 0;
      this.aMin.add(this.aMin, 'e').listen();
      this.aMinUpdate();
      
      // Currently displayed orbit (initialized to eMin)
      const dtMin     = Math.min(this.eMin.dt, this.aMin.dt);
      const dtMax     = Math.max(this.eMin.dt, this.aMin.dt);
      // Bound phi so that eccentricity remains less than one (elliptic)
      const e         = this.eMin.e;
      const phiMin    = Math.ceil(-Math.atan(Math.sqrt(1-e*e)/e)*180/Math.PI);
      const phiMax    = Math.floor(Math.atan(Math.sqrt(1-e*e)/e)*180/Math.PI);
      this.curr       = this.addFolder('Current Orbit');
      this.curr.show  = true;
      this.curr.add(this.curr, 'show').listen()
        .onChange(()=>{ controls.currUpdate(); });
      this.curr.dt    = this.eMin.dt;
      this.dtSlider   = this.curr.add(this.curr, 'dt', dtMin, dtMax, 1).listen();
      this.curr.phi   = 0.0;
      this.phiSlider  = this.curr.add(this.curr, 'phi', phiMin, phiMax, 1).listen();
      this.curr.a     = this.eMin.a;
      this.curr.add(this.curr, 'a').listen();
      this.curr.e     = this.eMin.e;
      this.curr.add(this.curr, 'e').listen();
      this.currUpdate();
    }
    distUpdate() {
      // Invariants (for a given R1 and R2)
      const r1  = this.r1  = dist(F1, R1);
      const r2  = this.r2  = dist(F1, R2);
      const r12 = this.r12 = dist(R1, R2);
    }
    //---------------------------
    // Fundamental elliptic orbit
    //---------------------------
    eMinUpdate() {
      const r1  = this.r1;
      const r2  = this.r2;
      const r12 = this.r12;
      // Eccentricity
      const e   = this.eMin.e = Math.abs(r1 - r2)/r12;
      // Update the minimum value for the eccentricity slider on the current orbit
      if (this.phiCurr) {
        const phiMin = Math.ceil(-Math.atan(Math.sqrt(1-e*e)/e)*180/Math.PI);
        const phiMax = Math.floor(Math.atan(Math.sqrt(1-e*e)/e)*180/Math.PI);
        this.phiCurr.min(phiMin).max(phiMax);
        this.curr.phi = Math.min(phiMax,Math.max(phiMin, this.curr.phi));
      }
      // Semi-major axis length
      const a   = this.eMin.a = (r1 + r2)/2;
      // The semi-major axis for the fundamental orbit is parallel to
      // the chord line connecting the two observation points.  The
      // orientation of the perifocal coordinate axes is such that the
      // point closer to the focus is on the periapsis side.
      const P   = this.eMin.pVec  = (r1 < r2 ? (R1 - R2)/r12 : (R2 - R1)/r12).Normalized;
      // The secondary coordinate axis is rotated 90 degrees from the
      // primary coordinate axis.
      const Q   = this.eMin.qVec  = P*1e12;
      // Secondary focus location
      const F2  = this.eMin.F2    = F1 - 2*a*e*P;
      // Ellipse origin location
      const O   = this.eMin.O     = (F1 + F2).Normalized;
      // Semi-parameter
      const q   = this.eMin.q     = a*(1-e*e);
      // Eccentricity vector
      const eVec  = this.eMin.eVec  = e*P;
      const th1   = Math.acos(dot((R1-F1),eVec)/(r1*e));
      const th2   = Math.acos(dot((R2-F1),eVec)/(r2*e));
      const E1    = 2*Math.atan(Math.sqrt((1-e)/(1+e))*Math.tan(th1/2));
      const E2    = 2*Math.atan(Math.sqrt((1-e)/(1+e))*Math.tan(th2/2));
      const dt    = this.eMin.dt = Math.sqrt(a*a*a/mu)*(E2 - E1 - e*(Math.sin(E2) - Math.sin(E1)));
      // Generate orbit ellipse...
      if (this.eMin.show) {
        // const K = new Kepler(e);
        this.eMin.orbit = [...Array(N)].map((x,i)=>()=> {
          // ... with uniform mean anomaly
          // const M = 2*Math.PI*i/N;
          // const E = K.solve(M);
          // const th = Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
          // ... with uniform eccentric anomaly
          // const E = 2*Math.PI*i/N;
          // const th = Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
          // ... with uniform true anomaly
          const th = 2*Math.PI*i/N;
          return Ellipse(F1, a, e, th, P, Q);
        });
      }
    }
    aMinUpdate() {
      const r1  = this.r1;
      const r2  = this.r2;
      const r12 = this.r12;
      //---------------------------
      // Minimum energy orbit
      //---------------------------
      // Semi-major axis length
      const a     = this.aMin.a     = (r1 + r2 + r12)/4;
      // Update the minimum value for the semi-major axis slider on the current orbit
      if (this.aCurr) this.aCurr.min(a);
      // Secondary focus location
      const F2    = this.aMin.F2    = lerp(R1, R2, (2*a-r1)/(4*a-r1-r2));
      // eccentricity of minimum energy orbit
      const e     = this.aMin.e     = dist(F1,F2)/(2*a);
      // periapsis direction (p)
      const P     = this.aMin.pVec  = (F1 - F2)/(2*e*a);
      // semi-parameter (semi-minor axis) direction (q)
      const Q     = this.aMin.qVec  = P*1e12;
      // Ellipse origin location
      const O     = this.aMin.O     = (F1 + F2).Normalized;
      // Semi-parameter
      const q     = this.aMin.q     = a*(1-e*e);
      // Eccentricity vector
      const eVec  = this.aMin.eVec  = e*P;
      const th1   = Math.acos(dot((R1-F1),eVec)/(r1*e));
      const th2   = Math.acos(dot((R2-F1),eVec)/(r2*e));
      const E1    = 2*Math.atan(Math.sqrt((1-e)/(1+e))*Math.tan(th1/2));
      const E2    = 2*Math.atan(Math.sqrt((1-e)/(1+e))*Math.tan(th2/2));
      const dt    = this.aMin.dt = Math.sqrt(a*a*a/mu)*(E2 - E1 - e*(Math.sin(E2) - Math.sin(E1)));
      // Generate orbit ellipse...
      if (this.aMin.show) {
        // const K = new Kepler(e);
        this.aMin.orbit= [...Array(N)].map((x,i)=>()=> {
          // ... with uniform mean anomaly
          // const M = 2*Math.PI*i/N;
          // const E = K.solve(M);
          // const th = Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
          // ... with uniform eccentric anomaly
          // const E = 2*Math.PI*i/N;
          // const th = Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
          // ... with uniform true anomaly
          const th = 2*Math.PI*i/N;
          return Ellipse(F1, a, e, th, P, Q)
        });
      }
    }
    currUpdate() {
      const r1  = this.r1;
      const r2  = this.r2;
      const r12 = this.r12;
      // Apse line rotation from fundamental ellipse
      const phi = this.curr.phi;
      // Eccentricity
      const ep  = this.eMin.e;
      const eq  = ep*Math.tan(phi*Math.PI/180);
      const e   = this.curr.e = Math.sqrt(ep*ep + eq*eq);
      // Eccentricity vector
      const eVec= this.curr.eVec  = ep*this.eMin.pVec + eq*this.eMin.qVec;
      // Semi-major axis length
      const a   = this.curr.a = (r1 + dot((R1-F1),eVec))/(1-e*e);
      // const a   = this.curr.a = Math.max(this.aMin.a, this.curr.a);
      // Primary axis direction
      const P   = this.curr.pVec  = (eVec/e).Normalized;
      // Secondary axis direction for the current orbit
      const Q   = this.curr.qVec  = P*1e12;
      // Secondary focus location
      const F2  = this.curr.F2    = F1 - 2*a*e*P;
      // Ellipse origin location
      const O   = this.curr.O     = (F1 + F2).Normalized;
      // Semi-parameter
      const q   = this.curr.q     = a*(1-e*e);
      const th1   = Math.acos(dot((R1-F1),eVec)/(r1*e));
      const th2   = Math.acos(dot((R2-F1),eVec)/(r2*e));
      const E1    = 2*Math.atan(Math.sqrt((1-e)/(1+e))*Math.tan(th1/2));
      const E2    = 2*Math.atan(Math.sqrt((1-e)/(1+e))*Math.tan(th2/2));
      // console.log(phi, th1, th2, E1, E2);
      const dt    = this.curr.dt = Math.sqrt(a*a*a/mu)*(E2 - E1 - e*(Math.sin(E2) - Math.sin(E1)));
      // Generate orbit ellipse...
      if (this.curr.show) {
        // const K = new Kepler(e);
        this.curr.orbit= [...Array(N)].map((x,i)=>()=> {
          // ... with uniform mean anomaly
          // const M = 2*Math.PI*i/N;
          // const E = K.solve(M);
          // const th = Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
          // ... with uniform eccentric anomaly
          // const E = 2*Math.PI*i/N;
          // const th = Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
          // ... with uniform true anomaly
          const th = 2*Math.PI*i/N;
          return Ellipse(F1, a, e, th, P, Q)
        });
      }
    }
    graphUpdate() {
      // ... with uniform eccentric anomaly
      const eMinEllipse = [...Array(N)].map((x,i)=>()=> {
        const th = trueAnomalyFromEccentric(2*Math.PI*i/N, ef);
        return Ellipse(F1, af, ef, th, Pf, Qf)
      });
      graph.push('r1: '  + this.r1.toFixed(3));
      graph.push('r2: '  + this.r2.toFixed(3));
      graph.push('r12: ' + this.r12.toFixed(3));
      graph.push(0xAA0000, 'a: '  + this.eMin.a.toFixed(3) + ', e: '  + this.eMin.e.toFixed(3) + ', dt: ' + this.eMin.dt.toFixed(3));
      graph.push(0x00AA00, 'a: '  + this.aMin.a.toFixed(3) + ', e: '  + this.aMin.e.toFixed(3) + ', dt: ' + this.aMin.dt.toFixed(3));
      graph.push(0x0000AA, 'a: '  + this.curr.a.toFixed(3) + ', e: '  + this.curr.e.toFixed(3) + ', dt: ' + this.curr.dt.toFixed(3));
      graph.push(0x000000, F1, 'F1');
      graph.push(0x444444, R1, 'R1');
      graph.push(0x444444, R2, 'R2');
      graph.push(0xAAAAAA, R1&R2);
      graph.push(0xAAAAAA, R1&R2|(F1+this.eMin.eVec));
      if (this.eMin.show) {
        graph.push(0xFFAAAA, this.eMin.O,  'O');
        graph.push(0xFFAAAA, this.eMin.F2, 'F2');
        graph.push(0xFFAAAA, F1 & this.eMin.F2);
        graph.push(0xAA0000, [F1,F1+this.eMin.eVec]);
        graph.push(0xAA0000, ...curve(this.eMin.orbit));
      }
      if (this.aMin.show) {
        graph.push(0x88AA88, this.aMin.O,  'O');
        graph.push(0x88AA88, this.aMin.F2, 'F2');
        graph.push(0x88AA88, F1 & this.aMin.F2);
        graph.push(0x00AA00, [F1,F1+this.aMin.eVec]);
        graph.push(0x00AA00, ...curve(this.aMin.orbit));
      }
      if (this.curr.show) {
        graph.push(0xAAAAFF, this.curr.O,  'O');
        graph.push(0xAAAAFF, this.curr.F2, 'F2');
        graph.push(0x8888AA, F1 & this.curr.F2);
        graph.push(0x0000AA, [F1,F1+this.curr.eVec]);
        graph.push(0x0000AA, ...curve(this.curr.orbit));
      }
    }
      
    //   Pf: undefined,
    //   Qf: undefined,
    //   pf: 0, // component of eccentricity vector in Pf direction
    //   qf: 0, // component of eccentricity vector in Qf direction
    //   eVec: pf*Pf + qf*Qf,
    //   phi: Math.atan2(eq,ep);
    //   a: 2.0,
    //   e: 0.5
    // };
  };
  
  const controls = new Controls();
  
  // NOTE: All values inside of the following lambda will be
  // recomputed whenever any of the above points are modified by the
  // user.
  document.body.appendChild(this.graph(() => {
    graph = [];
    controls.distUpdate();
    controls.eMinUpdate();
    controls.aMinUpdate();
    controls.currUpdate();
    controls.graphUpdate();
    
    // const orbit = [...Array(N)].map((x,i)=>()=> {
    //   // const th = trueAnomalyFromEccentric(2*Math.PI*i/N, e);
    //   const th = 2*Math.PI*i/N;
    //   return Ellipse(F1, a, e, th, P, Q)
    // });
    // // ... with uniform mean anomaly
    // // const FundEllipse = [...Array(N)].map((x,i)=>()=> {
    // //   const th = trueAnomalyFromMean(2*Math.PI*i/N, ef);
    // //   return Ellipse(F1, af, ef, th, Pf, Qf)
    // // });
                                         
    // const auxCircle   = [...Array(Nc)].map((x,i)=>()=>Circle(O, a, 2*Math.PI*i/Nc));
   
    // let items = [
      // 'phi: ' + params.phi.toFixed(3),
      // 'e: '   + e.toFixed(3),
      // 'a: '   + a.toFixed(3),
      // 'tf: '  + tFlight.toFixed(3),
      // // 'p1,q1: ' + p1.toFixed(3) + ',' + q1.toFixed(3),
      // // 'p2,q2: ' + p2.toFixed(3) + ',' + q2.toFixed(3),
      // 0x888888, O,  'O',
      // 0x888888, F2, 'F2',
      // 0xAAAAAA, F1&F2, F1&F2|O,
      // 0xAAAAAA, R1&R2, R1&R2|(F1+ef*Pf),
      // 0x884488, [F1,F1+0.25*P], 'p',  [F1,F1+0.25*Q], 'q',
      // 0xAA0000, ...curve(auxCircle),
      // 0xAA0000, [F1,(F1+eVec)], 'e',
      // 0x0000AA, ...curve(orbit),
      // 0xAAAAAA, [F2,R1], '2a-r1',
      // 0xAAAAAA, [F2,R2], '2a-r2',
      // 0x448844, [F1,R1], 'r1',
      // 0x448844, [F1,R2], 'r2',
      // 0x448844, [R1,R2], 'c',
    //   0x444444, F1, 'F1',
    //   0x44AA44, R1, 'R1',
    //   0x44AA44, R2, 'R2' ];
    return graph;
  }, {
    animate: true,
    grid: true,
    labels: false,
    lineWidth: 1,
    pointRadius: 0.75,
    fontSize: 0.75,
    scale: 0.5,
    width: 0.99*Math.min(window.innerWidth, window.innerHeight),
    height: 0.99*Math.min(window.innerWidth, window.innerHeight)
  }));
});
