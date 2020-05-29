$fa=1;
$fs=0.4;

module hopper(ht=1000, funnel_h=0.000999, funnel_r=0.000499, orifice_r=0.000249, hopper_center_yz=0.0005) {

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
    if(vis_mode){
        hopper(ht=50, funnel_h=9.99, funnel_r=4.99, orifice_r=2.49, hopper_center_yz=5);
    }
    else {
        difference(){
            ht=1000;
            cube(size=2*(ht+1), center=true);
            hopper(ht=ht, funnel_h=0.000999, funnel_r=0.000499, orifice_r=0.000249, hopper_center_yz=0.0005);
        }
    }
}

exa_hopper(vis_mode=true);
