twist = 45;
hole_scale = 0.5;
slices = 50;
module twist_box(degrees, x_scale = 1, y_scale = 1, z_scale = 1, offset=[0,0,0]) {
    difference() {
        linear_extrude(height = z_scale > 1 ? z_scale : 1, center = true, convexity = 10, twist = degrees, slices=slices) {
            translate(offset)
            scale([x_scale, y_scale]) {
                square(1, center=true);
            }
        }
        //translate([2,0,0]){
        if(z_scale < 1) {
            translate([-10,-10, z_scale/2]) cube([20,20,1]);
            translate([-10,-10, -z_scale/2 - 1]) cube([20,20,1]);
        }
 
    }
}


module twisted_menger(twist_degrees, depth=1, max_depth=2) {
    difference() {
        twist_box(twist_degrees);
        

        twist_box(twist_degrees, 1.5, hole_scale, hole_scale); // x
        twist_box(twist_degrees, hole_scale, 1.5, hole_scale); // y
        twist_box(twist_degrees, hole_scale, hole_scale, 1.001); // z
        
        //twisted_menger
    }
    
}


    twisted_menger(twist);
//scale(hole_scale/4) twisted_menger(twist);

//Stack
rotate([0, 0, -twist]) {
    translate([0, 0, 1]) {
//        twisted_menger(twist);      
        
        rotate([0, 0, -twist]) {
            translate([0, 0, 1]) {
//                twisted_menger(twist);      
                
                rotate([0, 0, -twist]) {
                    translate([0, 0, 1]) {
//                        twisted_menger(twist);      
                    }
                }
            }
        }
    }
}