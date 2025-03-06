// quick 2D PGA test.
Algebra(2,0,1,()=>{
  // const params = {
  //   phi: 0
  // };
  // const GUI = lil.GUI;
  // const gui = new GUI();
  // gui.add( document, 'title' );
  // gui.add( params, 'phi', -90, 90 );

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
    const ef  = Math.abs(r1 - r2)/c;  // eccentricity of the fundamental ellipse
    const Pf  = (r1 < r2 ? (R1 - R2)/c : (R2 - R1)/c).Normalized; // perifocal unit vector (p)
    const Ef  = ef*Pf;                // eccentricity vector
    const Qf  = Pf*1e12;              // semi-parameter (semi-minor axis) direction (q)
    // Minimum energy orbit
    const am  = (r1 + r2 + c)/4;
    const F2  = lerp(R1, R2, (2*am-r1)/(4*am-r1-r2));  // secondary focus location
    const O   = (F1 + F2).Normalized;       // ellipse origin location
    const em  = dist(F1,F2)/(2*am);         // eccentricity of minimum energy orbit
    const Em  = (F1 - F2)/(2*am);           // eccentricity vector for minimum energy orbit
    const Pm  = Em/em;                      // periapsis direction (p)
    const Qm  = Pm*1e12;                    // semi-parameter (semi-minor axis) direction (q)
    // perifocal unit vector (p)
    // const phi = Math.acos(Pf.Dot(Pm));
    // Generate fundamental ellipse...
    // ... with uniform true anomaly
    // const FundEllipse = [...Array(N)].map((x,i)=>()=> {
    //   const th = trueAnomalyFromMean(2*Math.PI*i/N, ef);
    //   return Ellipse(F1, af, ef, th, Pf, Qf)
    // });
    // ... with uniform eccentric anomaly
    // const eMinEllipse = [...Array(N)].map((x,i)=>()=> {
    //   const th = trueAnomalyFromEccentric(2*Math.PI*i/N, ef);
    //   return Ellipse(F1, af, ef, th, Pf, Qf)
    // });
    const aMinEllipse = [...Array(N)].map((x,i)=>()=> {
      const th = trueAnomalyFromEccentric(2*Math.PI*i/N, em);
      return Ellipse(F1, am, em, th, Pm, Qm)
    });
    // ... with uniform mean anomaly
    // const FundEllipse = [...Array(N)].map((x,i)=>()=> {
    //   const th = trueAnomalyFromMean(2*Math.PI*i/N, ef);
    //   return Ellipse(F1, af, ef, th, Pf, Qf)
    // });
                                         
    // const eAuxCircle   = [...Array(Nc)].map((x,i)=>()=>Circle(O, af, 2*Math.PI*i/Nc));
    const aAuxCircle   = [...Array(Nc)].map((x,i)=>()=>Circle(O, am, 2*Math.PI*i/Nc));
   
    return [
      'r1:  ' + r1.toFixed(3),
      'r2:  ' + r2.toFixed(3),
      'c:   ' + c.toFixed(3),
      // 'af:  ' + af.toFixed(3),
      // 'ef:  ' + ef.toFixed(3),
      'am:  ' + am.toFixed(3),
      'em:  ' + em.toFixed(3),
      //'phi: ' + (phi*180/Math.PI).toFixed(3),
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
      [F2,R1],// 'r1',
      [F2,R2],// 'r2',
      0x884488,
      [F1,F1+0.25*Pm], 'pm',
      [F1,F1+0.25*Qm], 'qm',
      // 0x888844,
      // [F1,F1+0.25*Pf], 'pf',
      // [F1,F1+0.25*Qf], 'qf',
      0xAA0000,
      ...curve(aAuxCircle),
      [F1,(F1+em*Pm)], 'em',
      // 0xAA00AA,
      // [F1,(F1+ef*Pf)], 'ef',
      // 0x00AA00,
      // ...curve(eMinEllipse),
      0x0000AA,
      ...curve(aMinEllipse),
    ]}, { grid: true, scale: 1, width: window.innerWidth, height: window.innerHeight }));
});

// // https://enkimute.github.io/ganja.js/examples/coffeeshop.html#DSc4ck8rX
// Algebra(2,0,1,()=>{
//   // Some constants. (scientific notation with uppercase E not overloaded)
//   const G          = 6.6703E-11; // N*m2/kg2 
//   const mEarth     = 5.97237E24; // kg
//   const mMoon      = 7.342E22;   // kg
//   const vMoon      = 1085;       // m/s
//   const dEarthMoon = 362400E3;   // m
//   const CoG        = dEarthMoon*mMoon/(mEarth+mMoon); // CoG of earth/moon in earth coords.
  
//   // Using Geometric/Clifford algebra with signature (2,0,1), based on the work of Charles Gunn
//   var point   = (x,y)=>!(1e0 + x*1e1 + y*1e2);
//   var vector  = (x,y)=>!(x*1e1 + y*1e2);
//   var dist    = (P,Q)=>((P.Normalized)&(Q.Normalized)).Length;
//   var refl    = (l,m)=>(l.Normalized<<m.Normalized)*2*(m.Normalized)-l.Normalized;
//   var project = (a,b)=>(a | b) / b;
//   var reject  = (a,b)=>(a | b);
  
//   var toString = (v,n)=>v.toPrecision(n||1).toString();
  
//   // Compute eccentric anomaly from mean anomaly and eccentricity
//   const maxIter=10;
//   var nitr=0;
//   var eccentricAnomaly = (M,e) => {
//     var E = (M<Math.PI ? M + e/2 : M - e/2);
//     for (let i=0; i<maxIter; ++i) {
//       var f = E - e*Math.sin(E) - M;
//       var df = 1 - e*Math.cos(E);
//       var de = f/df;
//       if (Math.abs(de) < 0.0001) {nitr=i; break;}
//       E -= de;
//     }
//     return E;
//   }
  
//   // Compute true anomaly from mean anomaly and eccentricity
//   var trueAnomaly = (t,T,e) => {
//     const M = 2*Math.PI*t/T;
//     const E = eccentricAnomaly(M,e);
//     return Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
//   }
  
//   // Compute true anomaly from mean anomaly and eccentricity
//   var trueAnomalyFromMean = (M,e) => {
//     const E = eccentricAnomaly(M,e);
//     return Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
//   }
  
//   // Compute true anomaly from mean anomaly and eccentricity
//   var trueAnomalyFromEccentric = (E,e) => {
//     return Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
//   }
  
//   var mu = 0.00001;
//   var rp = 0.5;           // orbital radius at periapsis
//   var ra = 2.0;           // orbital radius at apoapsis
//   var a0 = 3.0; //0.5*(ra+rp);    // semi-major axis
//   var e0 = 0.5; //(ra-rp)/(ra+rp);// eccentricity
//   // when theta = 0, r = rp = h^2/(mu*(1+e*cos(theta)))
//   //var h = Math.sqrt(rp*mu*(1+e0)); // specific angular momentum
//   var F1=point( a0*e0,0);         // primary focus
//   var F2=point(-a0*e0,0);         // secondary focus
//   // var O = ()=>0.5*(F1,F2);        // origin
//   // var p = ()=>(F1-F2)/dist(F1,F2);// eccentricity vector
//   // var q = ()=>(F1&F2)|O;
//   //var Radius=(t)=>a*(1-e*e)/(1+e*Math.cos(t));   // equation of the ellipse
//   //var TH = [ 0, Math.PI, 2*Math.PI/3, 4*Math.PI/3, 5*Math.PI/3 ];
//   //var R = TH.map((t,i)=>Radius(t));
//   //var A=F+vector(R[0]*Math.cos(TH[0]),R[0]*Math.sin(TH[0]));
//   //var B=F+vector(R[1]*Math.cos(TH[1]),R[1]*Math.sin(TH[1]));
//   //var C=F+vector(R[2]*Math.cos(TH[2]),R[2]*Math.sin(TH[2]));
//   //var D=F+vector(R[3]*Math.cos(TH[3]),R[3]*Math.sin(TH[3]));
//   //var E=F+vector(R[4]*Math.cos(TH[4]),R[4]*Math.sin(TH[4]));

//   // Closed curve
//   var curve=(pts)=>pts.map((pt,i,arr)=>[pt, pts[i+1%arr.length]]);
//   //var OpenCurve=(pts)=>pts.map((pt,i,arr)=>[(i?pts[i-1]:pt), pt]);

//   // Computes the location of a point on a circle of radius, r, centered at R0
//   var Circle=(R0,r,theta)=>{
//     return ()=>R0+vector(r*Math.cos(theta), r*Math.sin(theta));
//   };
//   var Ellipse=(R0, a, e, theta, p, q)=>{
//     var r = a*(1-e*e)/(1+e*Math.cos(theta));
//     //return R0 + r*Math.cos(theta)*p + r*Math.sin(theta)*q;
//     return ()=>R0 + vector(r*Math.cos(theta), r*Math.sin(theta));
//   }
//   var N=100, Nc=60
  
//   //var circle_pts = [...Array(Nc)].map((x,i)=>()=>O+vec(r*Math.cos(i*2*Math.PI/Nc), r*Math.sin(i*2*Math.PI/Nc)));
//   //var circle = curve(circle_pts);
  
//   //var conic_pts = [...Array(N)].map((x,i)=>()=>Conic(i*Math.PI/N));
//   //var conic = curve(conic_pts);
//   //var conic_pts = [...Array(N)].map((x,i)=>()=>Radius(i*Math.PI/N));
//   //var conic = curve(conic_pts);
  
//   //var poncelet=(P,l)=>{   // Poncelet map
//   //  var pascal=()=>((P&D)^(B&C))&(l^(B&E)),
//   //    Q=()=>((pascal^(D&E))&C)^l;
//   //  return [Q,refl(l,Q&O)]
//   //};
//   const th1 = Math.random()*2*Math.PI;
//   const R1  = Ellipse(F1, a0, e0, th1, vector(1,0), vector(0,1));
//   const th2 = Math.random()*2*Math.PI;
//   const R2  = Ellipse(F1, a0, e0, th2, vector(1,0), vector(0,1));
  
//   document.body.appendChild(this.graph(()=>{
//     var e   = e0;
//     var O   = (F1+F2).Normalized;       // origin
//     var p   = (F1-F2)/dist(F1,F2);      // periapsis direction
//     var q   = (F1&F2)|O;                // semi-parameter direction
//     var f   = dist(F1,F2);            // focal distance
//     var a   = f/(2*e);                // semi-major axis length
//     var T   = 2*Math.PI*Math.sqrt(a*a*a/mu); //20000;
//     var n   = 2*Math.PI/T;
//     var ev  = (F1-F2)/f;              // eccentricity vector (unit)
//     var Lp  = F1&F2;                  // Major axis
//     var Lq  = reject(F1&F2,O);        // Minor axis
//     var p   = vector(1,0); //(F1-F2)/dist(F1,F2);    // major axis unit vector
//     var q   = vector(0,1);
//     var c0  = [...Array(Nc)].map((x,i)=>()=>Circle(O, a, 2*Math.PI*i/Nc));
//     var C0  = curve(c0);              // Circle of radius a centered on O
//     //var c00 = [...Array(Nc)].map((x,i)=>()=>Circle(O, 2*a, 2*Math.PI*i/Nc));
//     //var C00 = curve(c00);              // Circle of radius a centered on F1
//     //var c1  = [...Array(Nc)].map((x,i)=>()=>Circle(F1, a, 2*Math.PI*i/Nc));
//     //var C1  = curve(c1);              // Circle of radius 2a centered on F2
//     var c11 = [...Array(Nc)].map((x,i)=>()=>Circle(F1, 2*a, 2*Math.PI*i/Nc));
//     var C11 = curve(c11);              // Circle of radius 2a centered on F2
//     //var c2  = [...Array(Nc)].map((x,i)=>()=>Circle(F2, a, 2*Math.PI*i/Nc));
//     //var C2  = curve(c2);              // Ellipse computed using the orbit eqn.
//     var c22 = [...Array(Nc)].map((x,i)=>()=>Circle(F2, 2*a, 2*Math.PI*i/Nc));
//     var C22 = curve(c22);              // Ellipse computed using the orbit eqn.
//     var orbitPts = [...Array(N)].map((x,i)=>()=>Ellipse(F1, a, e, 2*Math.PI*i/N), p, q);
//     var orbitCrv = curve(orbitPts);             // Ellipse computed using the orbit eqn.

//     var t   = 0.25*performance.now()%T;              // current time
//     var M   = n*t;                              // Mean anomaly
//     var E   = eccentricAnomaly(M, e);           // Eccentric anomaly
//     var theta = trueAnomaly(t, T, e);           // True anomaly
//     //var theta = trueAnomalyFromMean(M, e);      // True anomaly
//     //var theta = trueAnomalyFromEccentric(E, e); // True anomaly
    
//     var P0  = Circle(O, 2*a, M);     // Point on circlen centered at F1
//     var P1  = Circle(F1, 2*a, M);     // Point on circlen centered at F1
//     var P2  = Circle(F2, 2*a, M);     // Point on circlen centered at F1
//     var PM  = Circle(O, a, M);
//     var PE  = Circle(O, a, E);
//     var PE2  = Circle(O, 2*a, E);
    
//     var r   = a*(1-e*e)/(1+e*Math.cos(theta));
//     var PR  = Ellipse(F1, a, e, theta, p, q);
//     var L11  = F1&P1;                 // Line through F1 and P1
//     var L12  = F1&P2;                 // Line through F1 and P2
//     var L21  = F2&P1;                 // Line through F2 and P1
//     var L22  = F2&P2;                 // Line through F2 and P2
//     var m11 = 0.5*(F1+P1);            // Midpoint on [F1,P1]
//     var m12 = 0.5*(F1+P2);            // Midpoint on [F1,P2]
//     var m21 = 0.5*(F2+P1);            // Midpoint on [F2,P1]
//     var m22 = 0.5*(F2+P2);            // Midpoint on [F2,P2]
//     var L1  = reject(L12, m21);       // Line perpendicular to L12 through m21
//     var L2  = reject(L21, m12);       // Line perpendicular to L21 through m12
//     var Q1  = PR + (2*a-r)*(PR-F1)/r; // Point that is 2a away from F1
//     var Q2  = PR + r*(PR-F2)/(2*a-r); // Point that is 2a away from F2
//     var Z1 = 0.5*(F2+Q1);
//     var Z2 = 0.5*(F1+Q2);
//     var H1 = (Q1&Q2)^(Z1&Z2);
//     //var Y1 = (H1&PE)^(F1&Q2);
//     //var Y2 = (H1&PE)^(F2&Q1);
    
//     //var P4 = L3 ^ L2;               // Point on the ellipse (L2 intersects L3)
//     //var P4 = project(P3, F1&F2);    // Project P3 onto the semi-major axis
//     //var b  = 1/Math.sqrt(1-e*e);    // Semi-minor axis length
//     //var P5 = (1-b)*P4 + b*P3;       // Project P3 onto a circle of radius, a
//     return [
//       //"Ellipse construction",
//       "e = " + toString(e,2),
//       "a = " + toString(a,2),
//       "T = " + toString(T,5),
//       "M = " + toString(180*M/Math.PI,3),
//       "E = " + toString(180*E/Math.PI,3),
//       "theta = " + toString(180*theta/Math.PI,3),
//       //"F2:Y2:Q1 = " + toString(dist(F2,Y2)/dist(F2,Q1),3),
//       //"F1:Y1:Q2 = " + toString(dist(F1,Y1)/dist(F1,Q2),3),
//       0xAA0088, ...C0, //...C00, //P0, 'P0',
//       //0x444444, ...C11, //...C1, //P1, 'P1', [F1,P1], 0xAAAAAA, [F1,P2], m12, 'm', //'m12', m11, 'm11', // Note: m12 = m21
//       //0x444444, ...C22, //...C2, //P2, 'P2', [F2,P2], 0xAAAAAA, [F2,P1], // m21, 'm21', m22, 'm22',
//       0x00AA88, ...orbitCrv,
//       0x444444, 
//       O, //'O', 
//       F1, //'F1', 
//       F2, //'F2',
//       //0x884488, 
//       //[O,O+p], 'p',
//       //[O,O+q], 'q',
//       // O&PM,
//       //[PE,PM],
//       //PE&PM,
//       //F1&PM,
//       //F2&PM,
//       0x000044, 
//       [O,PM], 'M', 
//       //PM, //O&PM,
//       0xAA0088, 
//       [O,PE], 'E', 
//       //PE, //PE2, //O&PE,
//       0x44AA44, 
//       [F1,PR], 'r', 
//       //[PR,Q1], '2a-r', Q1, 'Q1',
//       //0xAAAA00,
//       //[F2,PR], '2a-r',
//       //[PR,Q2], 'r', Q2, 'Q2',
//       //0x00AA88, 
//       PR,
//       //[O,PR],
//       0xAA44AA,
//       R1, 
//       R2,
//       [R1,R2], 'c',
//       //[F2,Q1], [F1,Q2], 
//       PR&PE,
//       //(PR&PE)^(O&PM),
//       //reject(PR&PE,PE),
//       //reject(PR&PE,PE)^(O&PM),
//       //[Z1,Z2],
//       //Z1&Z2,
//       //[Z1,Q2],
//       //[Z2,Q1],
//       //[Q1,Q2], //Q1&Q2,
//       //Z1, 'Z1', Z2, 'Z2',
//       //reject(O&PE,PE),
//       // reject(O&PE,PE2),
//       //H1,
//       //H1&PE,
//       //Y1, Y2,
//       //H1&PM,
//       //(H1&PE)^(O&PM),
//       //[O,Q1],[O,Q2],
//       //0xAA8800, [O,P5], P5, 'P5', 
//       //0xAA8800, //P2, //'P2',
//       //0x88AA00, O&P5, //[O,Q],
//       //...circle,
//       //P1, 'P1', P2, 'P2', //P3, 'P3', P4, 'P4', //P5, 'P5', (P3&P4),
//       //(Q2&Q3),
//       //[F1,P1],[F2,P1],//[O,P3],
//       //0x448800, L1, '(F1&P2)|m',
//       //L1^L22,
//       //0x884400, L2, '(F2&P1)|m', 
//       //L2^L11,
//       //0x444444, PR&PE,
//       //P0&PE,
//       //P0&PR,
//       //P1&PR, 
//       //P1&PE,
//       //P2&PR, 
//       //P2&PE,
//       //PM&PE,
//       //PM&PR,
//       //...polygon
//     ]}, {animate:true, grid:true, lineWidth:1, scale: 1, width: window.innerWidth, height: window.innerHeight}));
    
// })
