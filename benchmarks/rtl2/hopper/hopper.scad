$fn=100;

module hopper(ht, funnel_h, funnel_r, orifice_r, hopper_center_yz) {
    funnel_extension = orifice_r*funnel_h/(funnel_r - orifice_r);
    overlap = 1e-6*funnel_h;

    translate([-funnel_extension, hopper_center_yz, hopper_center_yz]) {
       rotate([0,90,0]) {
          union(){
             cylinder(h = funnel_h + funnel_extension, r1 = 0, r2 = funnel_r, center = false);
             translate([0,0,funnel_h + funnel_extension - overlap]){
                cylinder(h = ht, r = funnel_r, center = false);
             }
          }
       }
    }
}

module exa_hopper(vis_mode=true) {
    domain = [0.008, 0.002, 0.002];
    ht = 0.009;
    funnel_h = 0.00199;
    funnel_r=0.00099;
    orifice_r=0.00049;
    hopper_center_yz=0.0010;
    
    scl = vis_mode ? 1e4 : 1.0; //scale for visualization
    
    scale(scl) {
       hopper(ht = ht, 
       funnel_h = funnel_h, 
       funnel_r = funnel_r, 
       orifice_r = orifice_r, 
       hopper_center_yz = hopper_center_yz);
    }
    
    if(vis_mode) {
        color("gray", 0.5) cube(size = scl * domain);
    }
}

exa_hopper(vis_mode=true);
