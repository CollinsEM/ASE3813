const GUI = lil.GUI;

// Using Geometric/Clifford algebra with signature (2,0,1), based on the work of Charles Gunn
Algebra(2,0,1,()=>{

  const point   = (x,y)=>!(1e0 + x*1e1 + y*1e2);
  const vector  = (x,y)=>!(x*1e1 + y*1e2);
  const dist    = (P,Q)=>((P.Normalized)&(Q.Normalized)).Length;
  const refl    = (l,m)=>(l.Normalized<<m.Normalized)*2*(m.Normalized)-l.Normalized;
  const dot     = (a,b)=>(a.e02*b.e02 + a.e01*b.e01);
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
  const F1  = point(0, 0);
  // Initial locations of the two observations (draggable)
  // const R1  = F1 + vector(r1_*Math.cos(t1_), r1_*Math.sin(t1_));
  // const R2  = F1 + vector(r2_*Math.cos(t2_), r2_*Math.sin(t2_));
  // const R1  = point( 1, 1);
  // const R2  = point(-1, 1);
  // console.log(R1);
  // console.log(R2);
  // Rendered items
  var graph = [];
  
  class Controls extends GUI {
    constructor() {
      super();
      // Clear the graph
      graph = [];
      
      // this.distUpdate();
      
      // Orbit 1
      this.orbit1 = this.addFolder('Orbit 1');
      this.orbit1.show = true;
      this.orbit1.add(this.orbit1, 'show').listen()
        .onChange((value)=>{ if (value) controls.updateOrbit1(); });
      this.orbit1.a = 1.0;
      this.orbit1.add(this.orbit1, 'a', 1.0, 5.0, 0.05).listen();
      this.orbit1.e = 0.0;
      this.orbit1.add(this.orbit1, 'e', 0.0, 0.8, 0.01).listen();
      this.orbit1.th = 0;
      this.orbit1.add(this.orbit1, 'th', 0, 360, 5).listen()
        .onChange((value)=>{ this.orbit2.th = (value + 180)%360 });
      this.updateOrbit1();
      this.orbit1.close();
      
      // Orbit 2
      this.orbit2 = this.addFolder('Orbit 2');
      this.orbit2.show = true;
      this.orbit2.add(this.orbit2, 'show').listen()
        .onChange((value)=>{ if (value) controls.updateOrbit2(); });
      this.orbit2.a = 2.5;
      this.orbit2.add(this.orbit2, 'a', 2.5, 8.0, 0.05).listen();
      this.orbit2.e = 0.0;
      this.orbit2.add(this.orbit2, 'e', 0.0, 0.8, 0.01).listen();
      this.orbit2.th = 180;
      this.orbit2.add(this.orbit2, 'th', 0, 360, 5).listen()
        .onChange((value)=>{ this.orbit1.th = (value + 180)%360 });
      this.updateOrbit2();
      this.orbit2.close();
      
      // Transfer Orbit
      this.transfer = this.addFolder('Orbit 3');
      this.transfer.show = true;
      this.transfer.add(this.transfer, 'show').listen()
        .onChange((value)=>{ if (value) controls.updateTransfer(); });
      const rp = this.orbit1.a*(1-this.orbit1.e);
      const ra = this.orbit2.a*(1+this.orbit2.e);
      this.transfer.a = 0.5*(rp+ra);
      this.transfer.add(this.transfer, 'a').listen();
      this.transfer.e = (ra-rp)/(ra+rp);
      this.transfer.add(this.transfer, 'e').listen();
      this.transfer.th = 0;
      this.transfer.add(this.transfer, 'th', 0, 360, 5).listen();
      this.updateTransfer();
      this.transfer.close();

      this.close();
    }
    //---------------------------
    // Plot orbit points
    //---------------------------
    plotOrbit(obj) {
      // const K = new Kepler(e);
      obj.orbit= [...Array(N)].map((x,i)=>()=> {
        // ... with uniform mean anomaly
        // const M = 2*Math.PI*i/N;
        // const E = K.solve(M);
        // const th = Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
        // ... with uniform eccentric anomaly
        // const E = 2*Math.PI*i/N;
        // const th = Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
        // ... with uniform true anomaly
        const th = 2*Math.PI*i/N;
        return Ellipse(F1, obj.a, obj.e, th, obj.pVec, obj.qVec);
      });
    }
    //---------------------------
    // Orbit 1
    //---------------------------
    updateOrbit1() {
      // const r1  = this.r1;
      // const r2  = this.r2;
      // const r12 = this.r12;
      // Eccentricity
      const e   = this.orbit1.e;
      // Semi-major axis length
      const a   = this.orbit1.a;
      // True anomaly of departure location
      const th  = this.orbit1.th*Math.PI/180;
      // Distance of body from planet at departure location
      const r   = a*(1-e*e)/(1+e*Math.cos(th));
      // Primary coordinate axis
      const P   = this.orbit1.pVec  = vector(1, 0);
      // Secondary coordinate axis, rotated 90 degrees from pVec.
      const Q   = this.orbit1.qVec  = P*1e12;
      // Secondary focus location
      const F2  = this.orbit1.F2    = F1 - 2*a*e*P;
      // Ellipse origin location
      const O   = this.orbit1.O     = (F1 + F2).Normalized;
      // Semi-parameter
      const q   = this.orbit1.q     = a*(1-e*e);
      // Eccentricity vector
      const eVec= this.orbit1.eVec  = e*P;
      // Location of the departure point
      const A   = this.A = Ellipse(F1, a, e, th, P, Q);
      // Specific mechanical energy: E = v*v/2 - mu/r = -mu/(2*a)
      // const E   = -mu/(2*a);
      // Orbital velocity at departure point on Orbit 1
      // const vSq = (2*mu/r - mu/a);
      // Specific angular momentum h*h/mu = a*(1-e*e)
      const h   = Math.sqrt(a*mu*(1-e*e));
      // Perpendicular velocity at departure point
      const vp  = h/r;
      // Radial velocity at departure point
      const vr  = (mu/h)*e*Math.sin(th);
      // Radial direction
      const rHat = (A-F1)/r;
      // Perpendicular direction
      const thHat = (rHat*1e12).Normalized;
      // Velocity vector
      this.vA = vr*rHat + vp*thHat;
      // Generate orbit ellipse...
      if (this.orbit1.show) this.plotOrbit(this.orbit1);
    }
    //---------------------------
    // Orbit 2
    //---------------------------
    updateOrbit2() {
      // const r1  = this.r1;
      // const r2  = this.r2;
      // const r12 = this.r12;
      // Eccentricity
      const e   = this.orbit2.e;
      // Semi-major axis length
      const a   = this.orbit2.a;
      // True anomaly of departure location
      const th  = this.orbit2.th*Math.PI/180;
      // Distance of body from planet at departure location
      const r   = a*(1-e*e)/(1+e*Math.cos(th));
      // Primary coordinate axis
      const P   = this.orbit2.pVec  = vector(1, 0);
      // Secondary coordinate axis, rotated 90 degrees from pVec.
      const Q   = this.orbit2.qVec  = P*1e12;
      // Secondary focus location
      const F2  = this.orbit2.F2    = F1 - 2*a*e*P;
      // Ellipse origin location
      const O   = this.orbit2.O     = (F1 + F2).Normalized;
      // Semi-parameter
      const q   = this.orbit2.q     = a*(1-e*e);
      // Eccentricity vector
      const eVec= this.orbit2.eVec  = e*P;
      // Location of the departure point
      const B   = this.B = Ellipse(F1, a, e, th, P, Q);
      // Specific mechanical energy: E = v*v/2 - mu/r = -mu/(2*a)
      // const E   = -mu/(2*a);
      // Orbital velocity at departure point on Orbit 1
      // const vSq  = (2*mu/r - mu/a);
      // Specific angular momentum h*h/mu = a*(1-e*e)
      const h   = Math.sqrt(a*mu*(1-e*e));
      // Perpendicular velocity at departure point
      const vp = h/r;
      // Radial velocity at departure point
      const vr = (mu/h)*e*Math.sin(th);
      // Radial direction
      const rHat = (B-F1)/r;
      // Perpendicular direction
      const thHat = (rHat*1e12).Normalized;
      // Velocity vector
      this.vB = vr*rHat + vp*thHat;
      // console.log(vr1, vp1, rHat, thHat);
      // Generate orbit ellipse...
      if (this.orbit2.show) this.plotOrbit(this.orbit2);
    }
    //---------------------------
    // Transfer Orbit
    //---------------------------
    updateTransfer() {
      // Orbit 1 - departure point
      const a1   = this.orbit1.a;
      const e1   = this.orbit1.e;
      const h1   = Math.sqrt(mu*a1*(1-e1*e1));
      const thA  = this.orbit1.th;
      const rA   = a1*(1-e1*e1)/(1+e1*Math.cos(thA));
      const v1A  = this.v1A  = Math.sqrt(2*mu/rA - mu/a1);
      const vr1A = this.vr1A = (mu/h1)*e1*Math.sin(thA);
      const vp1A = this.vp1A = Math.sqrt(v1A*v1A - vr1A*vr1A);
      // Orbit 2 - arrival point
      const a2   = this.orbit2.a;
      const e2   = this.orbit2.e;
      const h2   = Math.sqrt(mu*a2*(1-e2*e2));
      const thB  = this.orbit2.th;
      const rB   = a2*(1-e2*e2)/(1+e2*Math.cos(thB));
      const v2B  = this.v2B  = Math.sqrt(2*mu/rB - mu/a2);
      const vr2B = this.vr2B = (mu/h2)*e2*Math.sin(thB);
      const vp2B = this.vp2B = Math.sqrt(v2B*v2B - vr2B*vr2B);
      // Transfer orbit points
      const rp = Math.min(rA, rB);
      const ra = Math.max(rA, rB);
      // Semi-major axis length
      const a   = this.transfer.a = 0.5*(ra+rp);
      // Eccentricity
      const e   = this.transfer.e = (ra-rp)/(ra+rp);
      // True anomaly of C
      const th  = this.transfer.th*Math.PI/180;
      // Distance of body from planet at C
      const r   = a*(1-e*e)/(1+e*Math.cos(th));
      // Primary coordinate axis
      const P   = this.transfer.pVec  = (rA < rB ? (this.A - this.B) : (this.B - this.A))/(rA + rB);
      // Secondary coordinate axis, rotated 90 degrees from pVec.
      const Q   = this.transfer.qVec  = P*1e12;
      // Secondary focus location
      const F2  = this.transfer.F2    = F1 - 2*a*e*P;
      // Ellipse origin location
      const O   = this.transfer.O     = (F1 + F2).Normalized;
      // Semi-parameter
      const q   = this.transfer.q     = a*(1-e*e);
      // Eccentricity vector
      const eVec= this.transfer.eVec  = e*P;
      // Location of the departure point
      const C   = this.C = Ellipse(F1, a, e, th, P, Q);
      // Specific mechanical energy: E = v*v/2 - mu/r = -mu/(2*a)
      // const E   = -mu/(2*a);
      // Orbital velocity at departure point on Orbit 1
      // const vSq  = (2*mu/r - mu/a);
      // Specific angular momentum h*h/mu = a*(1-e*e)
      const h   = Math.sqrt(a*mu*(1-e*e));
      // Perpendicular velocity at C
      const vp = h/r;
      // Radial velocity at C
      const vr = (mu/h)*e*Math.sin(th);
      // Radial direction
      const rHat = (C-F1)/r;
      // Perpendicular direction
      const thHat = (rHat*1e12).Normalized;
      // Velocity vector
      this.vC = vr*rHat + vp*thHat;
      this.transfer.vp = h/rp;
      this.transfer.va = h/ra;
      this.vtA = h/rA;
      this.vtB = h/rB;
      this.dvA = this.vtA - this.v1A;
      this.dvB = this.v2B - this.vtB;
      // Generate orbit ellipse...
      if (this.transfer.show) this.plotOrbit(this.transfer);
    }
    graphUpdate() {
      graph.push('v1A: '  + this.v1A.toFixed(3));
      graph.push('vtA: '  + this.vtA.toFixed(3));
      graph.push('v2B: '  + this.v2B.toFixed(3));
      graph.push('vtB: '  + this.vtB.toFixed(3));
      graph.push('dvA: '  + this.dvA.toFixed(3));
      graph.push('dvB: '  + this.dvB.toFixed(3));
      // graph.push('r1: '  + this.r1.toFixed(3));
      // graph.push('r2: '  + this.r2.toFixed(3));
      // graph.push('r12: ' + this.r12.toFixed(3));
      // graph.push(0xAA0000, 'a: '  + this.orbit1.a.toFixed(3) + ', e: '  + this.orbit1.e.toFixed(3));
      // graph.push(0x00AA00, 'a: '  + this.aMin.a.toFixed(3) + ', e: '  + this.aMin.e.toFixed(3) + ', dt: ' + this.aMin.dt.toFixed(3));
      // graph.push(0x0000AA, 'a: '  + this.curr.a.toFixed(3) + ', e: '  + this.curr.e.toFixed(3) + ', dt: ' + this.curr.dt.toFixed(3));
      graph.push(0x000000, F1, 'F1');
      if (this.orbit1.show) {
        graph.push(0xAA0000, ...curve(this.orbit1.orbit));
        graph.push(0x444444, this.A, 'A');
        graph.push(0x884444, [this.A, this.A+this.vA]);
        // graph.push(0xFFAAAA, this.orbit1.O,  'O');
        // graph.push(0xFFAAAA, this.orbit1.F2, 'F2');
        // graph.push(0xFFAAAA, F1 & this.orbit1.F2);
        // graph.push(0xAA0000, [F1,F1+this.orbit1.eVec]);
      }
      if (this.orbit2.show) {
        graph.push(0x00AA00, ...curve(this.orbit2.orbit));
        graph.push(0x444444, this.B, 'B');
        graph.push(0x448844, [this.B, this.B+this.vB]);
        // graph.push(0xAAFFAA, this.orbit2.O,  'O');
        // graph.push(0xAAFFAA, this.orbit2.F2, 'F2');
        // graph.push(0xAAFFAA, F1 & this.orbit2.F2);
        // graph.push(0x00AA00, [F1,F1+this.orbit2.eVec]);
      }
      if (this.transfer.show) {
        graph.push(0x0000AA, ...curve(this.transfer.orbit));
        graph.push(0x444444, this.C, 'C');
        graph.push(0x448844, [this.C, this.C+this.vC]);
        // graph.push(0xAAFFAA, this.transfer.O,  'O');
        // graph.push(0xAAFFAA, this.transfer.F2, 'F2');
        // graph.push(0xAAFFAA, F1 & this.transfer.F2);
        // graph.push(0x00AA00, [F1,F1+this.transfer.eVec]);
      }
      // if (this.aMin.show) {
      //   graph.push(0x88AA88, this.aMin.O,  'O');
      //   graph.push(0x88AA88, this.aMin.F2, 'F2');
      //   graph.push(0x88AA88, F1 & this.aMin.F2);
      //   graph.push(0x00AA00, [F1,F1+this.aMin.eVec]);
      //   graph.push(0x00AA00, ...curve(this.aMin.orbit));
      // }
      // if (this.curr.show) {
      //   graph.push(0xAAAAFF, this.curr.O,  'O');
      //   graph.push(0xAAAAFF, this.curr.F2, 'F2');
      //   graph.push(0x8888AA, F1 & this.curr.F2);
      //   graph.push(0x0000AA, [F1,F1+this.curr.eVec]);
      //   graph.push(0x0000AA, ...curve(this.curr.orbit));
      // }
    }
  };
  
  const controls = new Controls();
  
  // NOTE: All values inside of the following lambda will be
  // recomputed whenever any of the above points are modified by the
  // user.
  document.body.appendChild(this.graph(() => {
    graph = [];
    controls.updateOrbit1();
    controls.updateOrbit2();
    controls.updateTransfer();
    controls.graphUpdate();
    return graph;
  }, {
    animate: true,
    grid: true,
    labels: false,
    lineWidth: 1,
    pointRadius: 0.75,
    fontSize: 0.75,
    scale: 0.5,
    width: 0.98*Math.min(window.innerWidth, window.innerHeight),
    height: 0.98*Math.min(window.innerWidth, window.innerHeight)
  }));
});
