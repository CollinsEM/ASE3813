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

  // Next we'll define some objects.
  var ptO       = up(0e1+0e2+0e3);                          // Sol pos (origin)
  var Sol       = sphere(ptO,0.1);                          // Sol (sun)
  var ecPlane   = plane(1e3,0);                             // Ecliptic Plane
  var orbit     = sphere(ptO,1.0)&ecPlane;                  // 1 AU orbit
  var ptE       = (t)=>up(1e1*Math.cos(t)+1e2*Math.sin(t)); // Earth pos
  var Earth     = (t)=>sphere(ptE(t),0.05);                 // Earth
  var eqPlane   = plane(0.397148e1+0.917755e3);
  var v1        = 1e2;
  var v2        = -0.917755e1+0.397148e3;
  var v3        = -0.397148e1-0.917755e3;
  var p1a       = (t)=>up(1e1*Math.cos(t)+1e2*Math.sin(t)-0.2*v1);
  var p1b       = (t)=>up(1e1*Math.cos(t)+1e2*Math.sin(t)+0.2*v1);
  var p2a       = (t)=>up(1e1*Math.cos(t)+1e2*Math.sin(t)-0.2*v2);
  var p2b       = (t)=>up(1e1*Math.cos(t)+1e2*Math.sin(t)+0.2*v2);
  var p3a       = (t)=>up(1e1*Math.cos(t)+1e2*Math.sin(t)-0.2*v3);
  var p3b       = (t)=>up(1e1*Math.cos(t)+1e2*Math.sin(t)+0.2*v3);
  var p4a       = (t)=>up(0e1-v3);
  var p4b       = (t)=>up(0e1+v3);
  var eq1       = (t)=>up(1e1*Math.cos(t)+1e2*Math.sin(t)-0.1e1+0.0397148e3+0.1e2);
  var eq2       = (t)=>up(1e1*Math.cos(t)+1e2*Math.sin(t)-0.1e1+0.0397148e3-0.1e2);
  var eq3       = (t)=>up(1e1*Math.cos(t)+1e2*Math.sin(t)+0.1e1-0.0397148e3-0.1e2);
  var eq4       = (t)=>up(1e1*Math.cos(t)+1e2*Math.sin(t)+0.1e1-0.0397148e3+0.1e2);
  var eq5       = (t)=>up(1e1*Math.cos(t)+1e2*Math.sin(t)+0.397148e1+0.917755e3);
  var eq6       = (t)=>up(1e1*Math.cos(t)+1e2*Math.sin(t)-0.397148e1-0.917755e3);
  var vernal    = ecPlane&eqPlane;
  // Our camera position and orientation
  
  // Graph the 3D items. (hex numbers are html5 colors, two extra first bytes = alpha)
  document.body.appendChild(this.graph(()=>{
    var t = performance.now()/4000;
    return [
      0xE0008800, ecPlane,                                  // Ecliptic plane
      0xF0008800, eqPlane,                                  // Equatorial plane
      0x00000000, orbit,                                    // 1 AU orbit
      0x00FFFF00, Sol,                                      // Sol (sun)
      0x000000FF, Earth(t),                                 // Earth
      0x00000000, vernal,
      0xAA000000, [eq1(t), eq2(t), eq3(t)],[eq1(t), eq3(t), eq4(t)],
      //0xAA000000, [eq5(t), eq6(t)], 
      0xAA000000, [p1a(t), p1b(t)], [p2a(t), p2b(t)], [p3a(t), p3b(t)], [p4a(t), p4b(t)], 
      0xAA000000, p1b(t), p2b(t), p3b(t), p4b(t)
    ];
  },{animate:true,conformal:true,gl:true,grid:false,labels:true,lineWidth:1})); //,camera})); 
});
  var canvas = document.getElementsByTagName('canvas');
  canvas[0].width = 1024;
  canvas[0].height = 768;
  console.log(canvas);
}
