window.addEventListener('load', init);

function init() {
// Create a Clifford Algebra with 4,1 metric for 3D CGA.
var Alg = Algebra(4,1,()=>{

  // We start by defining a null basis, and upcasting for points
  var ni = 1e4+1e5, no = .5e5-.5e4, nino = ni^no,
      up = (x)=>no+x+.5*(x*x)*ni,
      sphere = (P,r)=>!(P-r**2*.5*ni),
      plane  = (v,h=0)=>!(v-h*ni);

  // Project and reject.      
  var project_point_on_round            = (point,sphere)=>-point^ni<<sphere<<sphere,
      project_point_on_flat             = (point,plane)=>up(-point<<plane<<plane^nino*nino),
      plane_through_point_tangent_to_x  = (point,x)=>point^ni<<x*point^ni;

  const inc = Math.PI*23.4/180.0;
  const sin = Math.sin(inc);
  const cos = Math.cos(inc);
  
  // Next we'll define some objects.
  var sol       = up(0e1+0e2+0e3);                          // Sol pos (origin)
  var Sol       = sphere(sol,0.1);                          // Sol (sun)
  var ecPlane   = plane(1e2,0);                             // Ecliptic Plane
  
  var ec1       = up(1e1+1e2);
  var ec2       = up(1e1-1e2);
  var ec3       = up(-1e1-1e2);
  var ec4       = up(-1e1+1e2);
  var ecliptic1 = [ ec1, ec2, ec3 ];
  var ecliptic2 = [ ec1, ec3, ec4 ];
  // Principal directions of the ecliptic reference frame
  const u1 = 1e1, u2 = 1e3, u3 = 1e2;
  // Principal directions of the equatorial reference frame
  const v1 = 1e1, v2 = -sin*1e2 + cos*1e3, v3 = cos*1e2 + sin*1e3;
  
  var orbit     = sphere(sol,1.0)&ecPlane;                  // 1 AU orbit in ecliptic

  var x         = (t)=>1e1*Math.cos(t);
  var y         = (t)=>1e3*Math.sin(t);
  var earth     = (t)=>up(x(t)+y(t)); // Earth pos
  var Earth     = (t)=>sphere(earth(t),0.05);               // Earth
  var eqPlane   = plane(v3);
  var vernal    = ecPlane&eqPlane;                          // Vernal equinox direction
  
  // Points around the edge of the ecliptic plane centered on the Sun
  var ec1       = up(+0.5*u1+0.5*u2);
  var ec2       = up(+0.5*u1-0.5*u2);
  var ec3       = up(-0.5*u1-0.5*u2);
  var ec4       = up(-0.5*u1+0.5*u2);
  
  // Points around the edge of the equatorial plane centered on the Sun
  var eq1       = up(+0.5*v1+0.5*v2);
  var eq2       = up(+0.5*v1-0.5*v2);
  var eq3       = up(-0.5*v1-0.5*v2);
  var eq4       = up(-0.5*v1+0.5*v2);
  // Ecliptic plane pole points (centered on Sol)
  var s4a       = up(0e1-u3);
  var s4b       = up(0e1+u3);
  // Equatorial plane pole points (centered on Sol)
  var s5a       = up(0e1-v3);
  var s5b       = up(0e1+v3);
  
  // Geocentric reference frame prinicpal axes (GCRF)
  var p1a       = (t)=>up(x(t)+y(t)-0.4*v1);
  var p1b       = (t)=>up(x(t)+y(t)+0.4*v1);
  var p2a       = (t)=>up(x(t)+y(t)-0.4*v2);
  var p2b       = (t)=>up(x(t)+y(t)+0.4*v2);
  var p3a       = (t)=>up(x(t)+y(t)-0.4*v3);
  var p3b       = (t)=>up(x(t)+y(t)+0.4*v3);
  // Points around the edge of the ecliptic plane centered on the Earth
  var ec5       = (t)=>up(x(t)+y(t)+0.2*u1+0.2*u2);
  var ec6       = (t)=>up(x(t)+y(t)+0.2*u1-0.2*u2);
  var ec7       = (t)=>up(x(t)+y(t)-0.2*u1-0.2*u2);
  var ec8       = (t)=>up(x(t)+y(t)-0.2*u1+0.2*u2);
  // Points around the edge of the equatorial plane centered on the Earth
  var eq5       = (t)=>up(x(t)+y(t)+0.2*v1+0.2*v2);
  var eq6       = (t)=>up(x(t)+y(t)+0.2*v1-0.2*v2);
  var eq7       = (t)=>up(x(t)+y(t)-0.2*v1-0.2*v2);
  var eq8       = (t)=>up(x(t)+y(t)-0.2*v1+0.2*v2);
  
  // Graph the 3D items. (hex numbers are html5 colors, two extra first bytes = alpha)
  document.body.appendChild(this.graph(()=>{
    var t = performance.now()/4000;
    return [
      //0xE0008800, ecPlane,                                    // Ecliptic plane
      //0xE0008800, ecliptic1, ecliptic2,
      //0xF0008800, eqPlane,                                    // Equatorial plane
      0x00000000, orbit,                                        // 1 AU orbit
      0x00FFFF00, Sol,                                          // Sol (sun)
      0x000000FF, Earth(t),                                     // Earth
      0x00000000, vernal,
      0xAAFFFF00, [ec1, ec2, ec3],[ec1, ec3, ec4],              // Eclipitic plane
      0xAA0000FF, [eq1, eq2, eq3],[eq1, eq3, eq4],              // Equatorial plane
      0xAA000000, [s4a, s4b],                                   // Ecliptic pole
      0xAA000000, [s5a, s5b],                                   // Equatorial pole
      0x88FFFF00, [ec5(t), ec6(t), ec7(t)], [ec5(t), ec7(t), ec8(t)],
      0x880000FF, [eq5(t), eq6(t), eq7(t)], [eq5(t), eq7(t), eq8(t)],
      0xAA000000, [p1a(t), p1b(t)], [p2a(t), p2b(t)], [p3a(t), p3b(t)], // GCRF axes
    ];
  },{animate:true,conformal:true,gl:true,grid:false,labels:true,lineWidth:2})); 
});
  var canvas = document.getElementsByTagName('canvas');
  canvas[0].width = window.innerWidth;
  canvas[0].height = window.innerHeight;
  console.log(canvas);
}
