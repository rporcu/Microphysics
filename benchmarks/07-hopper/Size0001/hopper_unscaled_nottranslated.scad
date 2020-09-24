$fn=100;

module hopper(ht, funnel_h, funnel_r, orifice_r) {
    funnel_extension = orifice_r*funnel_h/(funnel_r - orifice_r);
    overlap = 1e-6*funnel_h;

    translate([-funnel_extension, 0, 0]) {
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
    ht = 50;
    funnel_h = 9.99;
    funnel_r= 4.99;
    orifice_r= 2.49;
    
    hopper(ht = ht, funnel_h = funnel_h, funnel_r = funnel_r, orifice_r = orifice_r);
}

exa_hopper();
