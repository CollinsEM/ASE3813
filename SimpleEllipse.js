// Create a Clifford Algebra with 3,1 metric for 2D CGA. 
Algebra(3,1,()=>{ 

  // The conformal model adds in more element types. (circles, point-pairs)
  // We no longer work in a dual space. (so ^ = join and & = meet)
  // Vectors are points, Bivectors are point pairs, Trivectors are lines/circles

  // We don't work directly in the e3/e4 basis, but instead rotate it so we have
  // two null vectors to work with (called origin and infinite)
  const ni = 1e4+1e3;           // n-infinity
  const no = .5e4-.5e3;         // n-origin
  
  const mu = 1.0;
  // const T = 2*Math.PI*Math.sqrt(a0*a0*a0/mu);
  const T = 10000;                  // Orbital period
  const n = 2*Math.PI/T;            // Mean motion
  
  // Define points, lines, circles using the null basis.  
  var point     = (x,y)=>no + x*1e1 + y*1e2 + 0.5*(x*x+y*y)*ni;
  var vec       = (x,y)=>!(x*1e1 + y*1e2);
  var line      = (a,b,c)=>!(a*1e1 + b*1e2 + c*ni);
  var circle    = (x,y,r)=>!(point(x,y) - r**2/2*ni);
  var ptRadius  = (c,r)=>!(c - r**2/2*ni);
  var randPoint = ()=>point(2*Math.random()-1,2*Math.random()-1);
  // Distances and Angles. 
  var dist=(x,y)=>(2*(x<<y).Length)**0.5;
  var angle=(x,y)=>Math.acos(!x.Normalized<<!y.Normalized);

  // Translation
  // v: vector displacement
  // returns a motor that will translate its target by v.
  var translate = (v)=>(1-0.5*v^ni);
  
  // Rotation
  // P: center of rotation
  // a: angle of rotation
  // returns a rotor that will rotate the target by a radians around P
  var rotate    = (P,a)=>Math.cos(a/2) - Math.sin(a/2)*(1e12-P<<1e12^ni);
  var a0 = 1.5;
  var e0 = 0.5;
  var l0 = a0*(1-e0*e0);
  var h = 1.0;

  // The location of the central body is at the primary focus
  var F1 = point( a0*e0,0);
  var F2 = point(-a0*e0,0);

  // The origin & eccentricity vector
  var O  = ()=>(F1+F2);

  // Define two points for which we want to locate an orbit passing through
  var th1 = Math.random()*2*Math.PI, r1 = l0/(1+e0*Math.cos(th1));
  var x1 = a0*e0 + r1*Math.cos(th1), y1 = r1*Math.sin(th1);
  // var x1 = a0*(Math.random()-0.5), y1 = (Math.random()-0.5);
  var R1 = point(x1, y1);               // First point relative to origin
  
  var th2 = Math.random()*2*Math.PI, r2 = l0/(1+e0*Math.cos(th2));
  var x2 = a0*e0 + r2*Math.cos(th2), y2 = r2*Math.sin(th2);
  var R2 = point(x2, y2);               // Second point relative to origin
  
  // Compute a rough estimate for dt
  var s1 = Math.sqrt(x1*x1 + y1*y1);
  var s2 = Math.sqrt(x2*x2 + y2*y2);
  var dE = Math.acos((x1*x2 + y1*y2)/(s1*s2));
  var dt0 = dE/n;
  
  // Find the orbit passing through R1 and R2 with time of flight, dt.
  // R1: initial position vector (relative to the primary focus)
  // R2: final position vector (relative to the primary focus)
  // dt: time of flight between R1 and R2
  var findOrbit = (R1,R2,dt)=> {
    // Radial distances from the primary focus
    var r1 = R1.Length;
    var r2 = R2.Length;
    // Length of the chord connecting R1 and R2
    var r12 = dist(R1,R2);
    // Semi-perimeter of the triangle formed by R1, R2, and the primary focus
    var s = (r1+r2+r12)/2;
    // Minimum eccentricity passing through R1 and R2
    var e0 = r12/(r1+r2);
    // Length of the semi-major axis of the minimum energy ellipse
    var a0 = s/2;
    
    // Compute difference in true anomaly
    var theta = ()=>Math.acos((R1<<R2)/(r1*r2));
    // Establish the direction perpendicular to the orbit
    var c12 = !(R1^R2);
    var A   = Math.sin(theta)*Math.sqrt((r1*r2)/(1-Math.cos(theta)));
    
    var alpha = 1/a0;
    var xi0 = Math.sqrt(mu)*Math.abs(alpha)*dt0;// \sqrt(mu) |alpha| dt
    var z0 = alpha*Math.pow(xi0,2); // alpha xi^2

    // Stumpff functions and their derivatives
    var C  = (z)=>(1/2)-z/24 + (z*z)/720 - (z*z*z)/40320;
    var S  = (z)=>(1/6)-z/120 + (z*z)/5040 - (z*z*z)/362880;
    var dC = (z)=>(1/(2*z))*(1-z*S(z)-2*C(z));
    var dS = (z)=>(1/(2*z))*(C(z)-3*S(z));

    // Lagrange coefficients
    var y  = (z)=>r1 + r2 + A*(z*S(z)-1)/Math.sqrt(C(z));
    var dy = (z)=>(A/(2*Math.sqrt(Math.pow(C(z),3))))*((1-z*S(z))*dC(z) + 2*(S(z) + z*dS(z))*C(z));
    var F = (z)=>Math.sqrt(Math.pow(y(z)/C(z),3))*S(z)+A*Math.sqrt(y(z))-Math.sqrt(mu)*t;
    var dF = (z)=>(1/(2*Math.sqrt(Math.pow(C(z),3))))*(3*y(z)*dC(z) + (1-z*S(z))*dy(z));

    var F = (z,dt,dth,r1,r2)=>Math.sqrt(Math.pow(y(z,dth,r1,r2)/C(z),3))*S(z)+A(dth,r1,r2)*Math.sqrt(y(z,dth,r1,r2))-Math.sqrt(mu)*dt;
    var dF = (z,dth,r1,r2)=>(1/(2*Math.sqrt(y(z,dth,r1,r2)*Math.pow(C(z),5))))*((2*C(z)*dS(z)-3*dC(z)*S(z))*Math.pow(y(z,dth,r1,r2),2) + (A(dth,r1,r2)*Math.sqrt(Math.pow(C(z),5))+3*C(z)*S(z)*y(z,dth,r1,r2))*dy(z,dth,r1,r2));
    
    var z = z0;
    for (let i=0; i<10; ++i) {
      var dz = -F(z,dt0,dTheta(),r1(),r2())/dF(z,dTheta(),r1(),r2());
      z += dz;
    }
  }
  
  //var F = 1 - (mu/h*h)*r2*(1-Math.cos(dTheta));
  //var G = (r1*r2/h)*(1 - Math.sin(dTheta));
  //var dF = (mu/h)*((1-Math.cos(dTheta))/Math.sin(dTheta))*((mu/(h*h))*(1-Math.cos(dTheta))-(1/r1)-(1/r2));
  //var dG = 1 - (mu*r1/(h*h))*(1-Math.cos(dTheta));

  //var r1 = ()=>R1-F1;                   // First point relative to focus
  //var r2 = ()=>R2-F2;                   // Second point relative to focus
  
  var S1 = ()=>(O-R1);                  // Reflect first point through origin
  var S2 = ()=>(O-R2);                  // Reflect second point through origin
  
  var T1 = ()=>F1-(R1-F1);
  var T2 = ()=>F1-(R2-F1);
  
  var b1 = ()=>dist(F2,R1);
  var b2 = ()=>dist(F2,R2);
  
  // Circles centered on primary points passing through primary focus
  var B1 = ()=>ptRadius(R1,dist(F1,R1));
  var B2 = ()=>ptRadius(R2,dist(F1,R2));
  
  // Circles centered on primary points passing through secondary focus
  var C1 = ()=>ptRadius(R1,dist(F2,R1));
  var C2 = ()=>ptRadius(R2,dist(F2,R2));
  
  // The distance between F1 and F2 is 2ae
  
  // Eccentricity is (r2-r1)/(r1*cos(theta1)-r2*cos(theta2))
  //var a = ()=>
  var e = ()=>e0;
  var l = ()=> a0*(1-e0*e0);
  var ev = ()=>(F1-O).Normalized;
  var ecc = ()=>[O,F1/a0];
  var nitr = 0, maxIter=10, de = 0;
  // Compute eccentric anomaly from mean anomaly and eccentricity
  var eccentricAnomaly = (M,e) => {
    var E = M;
    for (let i=0; i<maxIter; ++i) {
      var f = E - e*Math.sin(E) - M;
      var df = 1 - e*Math.cos(E);
      var de = f/df;
      if (Math.abs(de) < 0.0001) {nitr=i; break;}
      E -= de;
    }
    return E;
  }
  
  // Compute true anomaly from mean anomaly and eccentricity
  var trueAnomaly = (M,e) => {
    const E = eccentricAnomaly(M,e);
    return Math.atan2(Math.sin(E)*Math.sqrt(1-e*e),Math.cos(E)-e);
  }
  
  // Compute the location of a point on the origin
  var orbitPoint = (theta,e,a) => {
    var p0 = ()=>O;                   // origin
    var tr = ()=>translate(a*e*1e1);  // translator: shift origin to focus
    var r0 = ()=>rotate(p0,phi);      // rotor: rotate the apse line
    var p1 = tr>>>p0;                 // translate the origin to the focus along x-axis
    // r0 needs to be one of the following:
    var ro = tr*rotate(p0,theta); // rotate by theta around p0 then translate
    //var ro = rotate(p1,theta)*tr; // translate then rotate by theta around p1
    var r  = a*(1-e*e)/(1+e*Math.cos(theta)); // orbital radius
    var p2 = point(r,0);          // unrotated radial vector
    return ro>>>p2; // apply rotor to p2
  }
  
  var curve=(pts)=>pts.map((pt,i,arr)=>[pt, pts[i+1%arr.length]]);
  
  var N = 100;
  var conic_pts = [...Array(N)].map((x,i)=>()=>orbitPoint(2*i*Math.PI/N,e0,a0));
  var conic = curve(conic_pts);
  
  var toString = (v,n)=>v.toPrecision(n||1).toString();
  
  // Generate graph
  document.body.appendChild(this.graph(()=>{
    var t = performance.now()%T;    // Time t in [0,T-1]
    var M = n*t;                    // Mean anomaly
    var E = eccentricAnomaly(M, e0);// Eccentric anomaly
    var theta = trueAnomaly(M,e0);  // True anomaly
    //var r     = ()=>a0*(1-e0*e0)/(1+e0*Math.cos(theta));
    var r     = l0/(1+e0*Math.cos(theta));
    //var tr    = ()=>translate( Math.sin(theta)*1e1 );
    var tr    = translate(a0*e0*1e1);
    //var tr    = ()=>translate( ev );
    //var tr    = ()=>translate( F1-O );
    var p1    = tr>>>point( 0, 0 );
    var ro    = tr*rotate(p1,theta);
    //var P     = ()=>ro>>>(point( r, 0 ));
    var P     = orbitPoint(theta, e0, a0);
    //var phi   = angle(r1,r2);
    return [
      //  z0.toString(),
      //  z.toString(),
      //  "dot(R1,R2) = " + R1dotR2().toString(),
      //  "dTheta = " + dTheta().toString(),
      //  "A = " + A(dTheta(),r1(),r2()).toString(),
      "nitr = " + nitr,
      "de = " + toString(de,4),
      "dE = " + toString(dE,4),
      //"dt = " + dt.toString(),
      "t = " + t,
      "r = " + toString(r,4),
      "M = " + toString(M,4),
      "E = " + toString(E,4),
      "theta = " + toString(theta,4),
      //phi.toString(),
      //angle(,).toString(),
      //angle(1e1,r1).toString(),
      0x00AA88, ...conic,
      //0xFF8888, [R1,R2],
      0xFF8888, [F1,P],
      0x444444, P,'P', // A,'A', B,'B', C,'C', D,'D', E,'E', O,'O',
      // 0xFF0000, ecc,
      0xAA8888, B1, // "C1",                // circles
      0xAA8888, B2, // "C2",
      0xAA88AA, C1, // "C1",                // circles
      0xAA88AA, C2, // "C2",
      //0x44AA44, pp1, "pp1", //Y, "Y", p4,                 // lines
      0x4444FF, [F1,R1], [F1,R2], 
      0x44FFFF, [F2,R1], [F2,R2],
      //0x44FF44, [F1,F2],
      //0xAAAAAA, [R1,S1], [R2,S2],
      0x444444, O, "O",                   // origin
      0x44AA44, R1, "R1", R2, "R2",         // primary points
      //0x4444AA, S1, "S1", // "O-(A-O)",
      //0x4444AA, S2, "S2", // "O-(B-O)",
      //0x4444AA, T1, "F-(A-F)", T2, "F-(B-F)",       // reflected points
      //0x4444AA, [R1,T1], [R2,T2],
      0x666666, F1, "F", F2, "F'",        // focii
    ]},{conformal:true,grid:true,animate:true,width:640,height:480}));                 // conformal flag!  

});
