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
    domain = [0.004, 0.001, 0.001];
    ht = 0.005;
    funnel_h = 0.000999;
    funnel_r=0.000499;
    orifice_r=0.000249;
    hopper_center_yz=0.0005;
    
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
