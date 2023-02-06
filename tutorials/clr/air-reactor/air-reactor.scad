//!OpenSCAD

// constant dims (in m)
domain_overlap = 0.1;
air_reactor_diameter = 0.1524;
air_reactor_height = 0.7720;
air_reactor_c2c_height = 0.0762;
riser_diameter = 0.0635;
riser_height = 2.8602;
outlet_height = 3.6322;
outlet_length = 0.1392;
offset = 0.1392;
pipe_overlap = 1e-4*air_reactor_height;

translate([0, offset, offset]){ 
    rotate([0, 90, 0]){
        union(){
            // air reactor
            translate([0, 0, -domain_overlap]){
               cylinder($fn=50, 
                        r1=0.5*air_reactor_diameter, 
                        r2=0.5*air_reactor_diameter, 
                        h=air_reactor_height+domain_overlap+pipe_overlap, 
                        center=false); 
             }
            // air-reactor to riser neck 
            translate([0, 0, air_reactor_height]){
               cylinder($fn=50, 
                        r1=0.5*air_reactor_diameter, 
                        r2=0.5*riser_diameter, 
                        h=air_reactor_c2c_height, 
                        center=false); 
            }        
            // riser
            translate([0, 0, air_reactor_height+air_reactor_c2c_height-pipe_overlap]){
               cylinder($fn=50, 
                        r1=0.5*riser_diameter, 
                        r2=0.5*riser_diameter, 
                        h=riser_height,
                        center=false); 
            }
            //outlet
            translate([0, 0, outlet_height]){
                rotate([0,-90,0]){
                   cylinder($fn=50, 
                            r1=0.5*riser_diameter, 
                            r2=0.5*riser_diameter, 
                            h=outlet_length+domain_overlap,
                            center=false);
                } 
            }        
        }//union
    }
}
