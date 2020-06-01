module exa_cyclone()
{
   union() {
      translate([17,20,5]) cube(size = [18,40,7], center = false);
      translate([0,20,20]){
         union() {
            difference() {
               translate([15,0,0]) rotate(a=[0,90,0]) cylinder(40, 15, 15, center = true);
               translate([35,0,0]) rotate(a=[0,90,0]) cylinder(40, 8, 8, center = true);
            }
            translate([35,0,0]) rotate(a=[0,90,0]) cylinder(40, 5.5, 5.5, center = true);
         }
      }
   }
}
exa_cyclone();
