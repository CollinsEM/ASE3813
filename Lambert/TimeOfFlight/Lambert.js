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

  const toString= (v,n)=>v.toPrecision(n||1).toString();
  
  // Compute the location of a point on an ellipse
  // F0: primary focus location
  // a: semi-major axis length
  // e: eccentricity
  // t: true anomaly
  // p: semi-major axis
  // q: semi-parameter axis
  // const Ellipse = (F0, a, e, t, P, Q)=>{
  //   var r = a*(1-e*e)/(1+e*Math.cos(t));
  //   return F0 + r*Math.cos(t)*P + r*Math.sin(t)*Q;
  // }
  
  // Computes the location of a point on a circle of radius, r, centered at R0
  var Circle=(R0,r,theta)=>{
    return ()=>R0+vector(r*Math.cos(theta), r*Math.sin(theta));
  };

  // Closed curve
  var curve=(pts)=>pts.map((pt,i,arr)=>[pt, pts[i+1%arr.length]]);
  //var OpenCurve=(pts)=>pts.map((pt,i,arr)=>[(i?pts[i-1]:pt), pt]);

  var N=360, Nc=360;
  // Planet constants
  const mu = 1;
  const planet_radius = 0.1;
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
      this.eMin.close();
      
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
      this.aMin.add(this.aMin, 'phi').listen();
      this.aMin.close();
      
      // Currently displayed orbit (initialized to eMin)
      this.curr       = this.addFolder('Current Orbit');
      this.curr.show  = true;
      this.curr.add(this.curr, 'show').listen()
        .onChange(()=>{ controls.currUpdate(); });
      // Bound phi so that eccentricity remains less than one (elliptic)
      const e         = this.eMin.e;
      const phiMin    = -89;//Math.ceil(-Math.atan(Math.sqrt(0.99999999-e*e)/e)*180/Math.PI);
      const phiMax    =  89;//Math.floor(Math.atan(Math.sqrt(1-e*e)/e)*180/Math.PI);
      this.curr.phi   = 0.0;
      this.phiSlider  = this.curr.add(this.curr, 'phi', phiMin, phiMax, 0.1).listen();
      // const dtMin     = Math.min(this.eMin.dt, this.aMin.dt);
      // const dtMax     = Math.max(this.eMin.dt, this.aMin.dt);
      this.curr.dt    = this.eMin.dt;
      this.dtSlider   = this.curr.add(this.curr, 'dt').listen();
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
      this.Planet = [...Array(Nc)].map((x,i)=>()=>Circle(F1, planet_radius, 2*Math.PI*i/Nc));
    }
    drawOrbit(obj) {
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
        var r = obj.q/(1+obj.e*Math.cos(th));
        return F1 + r*Math.cos(th)*obj.pVec + r*Math.sin(th)*obj.qVec;
        // return Ellipse(F1, obj.a, obj.e, th, obj.pVec, obj.qVec);
      });
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
      // if (this.phiCurr) {
      //   const phiMin = Math.ceil(-Math.atan(Math.sqrt(1-e*e)/e)*180/Math.PI);
      //   const phiMax = Math.floor(Math.atan(Math.sqrt(1-e*e)/e)*180/Math.PI);
      //   this.phiCurr.min(phiMin).max(phiMax);
      //   this.curr.phi = Math.min(phiMax,Math.max(phiMin, this.curr.phi));
      // }
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
      if (this.eMin.show) this.drawOrbit(this.eMin);
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
      // Rotation angle of apse line from fundamental ellipse
      const phi   = this.aMin.phi   = Math.atan2(dot(this.aMin.pVec,this.eMin.qVec),dot(this.aMin.pVec,this.eMin.pVec))*180/Math.PI;
      // Semi-parameter
      const q     = this.aMin.q     = a*(1-e*e);
      // Eccentricity vector
      const eVec  = this.aMin.eVec  = e*P;
      const th1   = Math.acos(dot((R1-F1),eVec)/(r1*e));
      const th2   = Math.acos(dot((R2-F1),eVec)/(r2*e));
      const E1    = 2*Math.atan(Math.sqrt((1-e)/(1+e))*Math.tan(th1/2));
      const E2    = 2*Math.atan(Math.sqrt((1-e)/(1+e))*Math.tan(th2/2));
      const dt   = this.aMin.dt = Math.sqrt(a*a*a/mu)*(E2 - E1 - e*(Math.sin(E2) - Math.sin(E1)));
      // Generate orbit ellipse...
      if (this.aMin.show) this.drawOrbit(this.aMin);
    }
    currUpdate() {
      const r1  = this.r1;
      const r2  = this.r2;
      const r12 = this.r12;
      // Apse line rotation from fundamental ellipse
      const phi = this.curr.phi;
      // Eccentricity
      const ep  = this.eMin.e;
      // const eq  = Math.min(ep*Math.tan(phi*Math.PI/180), Math.sqrt(0.99999999-ep*ep));
      const eq  = ep*Math.tan(phi*Math.PI/180);
      const e   = this.curr.e = Math.sqrt(ep*ep + eq*eq);
      // Bound phi so that eccentricity remains less than one (elliptic)
      // const phiMin = Math.ceil(-Math.atan(Math.sqrt(0.99999999-e*e)/e)*180/Math.PI);
      // const phiMax = Math.floor(Math.atan(Math.sqrt(0.99999999-e*e)/e)*180/Math.PI);
      // this.phiSlider.min(phiMin).max(phiMax);
      // this.curr.phi = Math.min(phiMax,Math.max(phiMin, this.curr.phi));
      // Range of dt values between the eMin and aMin orbits
      // const dtMin     = Math.min(this.eMin.dt, this.aMin.dt);
      // const dtMax     = Math.max(this.eMin.dt, this.aMin.dt);
      // this.dtSlider.min(dtMin).max(dtMax);
      // Eccentricity vector
      const eVec= this.curr.eVec  = ep*this.eMin.pVec + eq*this.eMin.qVec;
      // Semi-major axis length
      // const a   = this.curr.a = Math.max(this.aMin.a, (r1 + dot((R1-F1),eVec))/(1-e*e));
      const a   = this.curr.a = (r1 + dot((R1-F1),eVec))/(1-e*e);
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
      const h   = Math.sqrt(mu*q);
      const th1   = Math.acos(dot((R1-F1),eVec)/(r1*e));
      const th2   = Math.acos(dot((R2-F1),eVec)/(r2*e));
      if (a*(1-e) < planet_radius) {
        this.curr.dtString = "INVALID";
      }
      else {
        if (e < 1) {
          const E1  = 2*Math.atan(Math.sqrt((1-e)/(1+e))*Math.tan(th1/2));
          const E2  = 2*Math.atan(Math.sqrt((1-e)/(1+e))*Math.tan(th2/2));
          
          const M1  = E1 - e*Math.sin(E1);
          const M2  = E2 - e*Math.sin(E2);
          const dt  = this.curr.dt = Math.sqrt(a*a*a/mu)*(M2 - M1);
        }
        else if (e == 1) {
          const M1  = Math.tan(th1/2)/2 + Math.pow(Math.tan(th1/2),3)/6;
          const M2  = Math.tan(th2/2)/2 + Math.pow(Math.tan(th2/2),3)/6;
          const dt  = this.curr.dt = (q*q/h)*(M2 - M1);
        }
        else if (e > 1) {
          const F1  = 2*Math.atanh(Math.sqrt((e-1)/(e+1))*Math.tan(th1/2));
          const F2  = 2*Math.atanh(Math.sqrt((e-1)/(e+1))*Math.tan(th2/2));
          const M1  = e*Math.sinh(F1) - F1;
          const M2  = e*Math.sinh(F2) - F2;
          const dt  = this.curr.dt = Math.sqrt(-a*a*a/mu)*(M2 - M1);
        }
        this.curr.dtString = this.curr.dt.toFixed(3);
      }
      // Generate orbit ellipse...
      if (this.curr.show) this.drawOrbit(this.curr);
    }
    graphUpdate() {
      // ... with uniform eccentric anomaly
      // const eMinEllipse = [...Array(N)].map((x,i)=>()=> {
      //   const th = trueAnomalyFromEccentric(2*Math.PI*i/N, ef);
      //   return Ellipse(F1, af, ef, th, Pf, Qf)
      // });
      graph.push('r1: '  + this.r1.toFixed(3));
      graph.push('r2: '  + this.r2.toFixed(3));
      graph.push('r12: ' + this.r12.toFixed(3));
      graph.push(0xAA0000, 'a: '  + this.eMin.a.toFixed(3) + ', e: '  + this.eMin.e.toFixed(3) + ', dt: ' + this.eMin.dt.toFixed(3));
      graph.push(0x00AA00, 'a: '  + this.aMin.a.toFixed(3) + ', e: '  + this.aMin.e.toFixed(3) + ', dt: ' + this.aMin.dt.toFixed(3));
      graph.push(0x0000AA, 'a: '  + this.curr.a.toFixed(3) + ', e: '  + this.curr.e.toFixed(3) + ', dt: ' + this.curr.dtString);
      graph.push(0x000000, F1, 'F1');
      graph.push(0x444444, R1, 'R1');
      graph.push(0x444444, R2, 'R2');
      graph.push(0xAA0000, ...curve(this.Planet)),
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
