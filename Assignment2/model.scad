twist = 30;
hole_scale = 0.5;
n_wide = 3;
n_long = 3;
n_tall = 3;

slices = 20;

module twist_box(degrees, x_scale = 1, y_scale = 1, z_scale = 1, offset=[0,0,0]) {
    difference() {
        linear_extrude(height = z_scale > 1 ? z_scale : 1, center = true, convexity = 10, twist = degrees, slices=slices) {
            translate(offset)
            scale([x_scale, y_scale]) {
                square(1, center=true);
            }
        }

        if(z_scale < 1) {
            translate([-2,-2, z_scale/2] + offset) cube([4,4,1]);
            translate([-2,-2, -z_scale/2 - 1] + offset) cube([4,4,1]);
        }
 
    }
}

module twist_cross(twist_degrees) {
    twist_box(twist_degrees, 1.5, hole_scale, hole_scale); // x
    twist_box(twist_degrees, hole_scale, 1.5, hole_scale); // y
    twist_box(twist_degrees, hole_scale, hole_scale, 1.001); // z
}

module twisted_menger_unit(twist_degrees, offset=[0,0,0], depth=1, max_depth=2,hole_scale=0.5) {
    difference() {
        twist_box(twist_degrees, offset=offset);
        
        twist_box(twist_degrees, 1.5, hole_scale, hole_scale, offset); // x
        twist_box(twist_degrees, hole_scale, 1.5, hole_scale, offset); // y
        twist_box(twist_degrees, hole_scale, hole_scale, 1.001, offset); // z
    }
}


module twisted_menger(twist_degrees, hole_scale, n_wide, n_long, n_tall) {
    union() {
        for(j=[0:n_tall - 1]) {
            translate([0,0,j]) {
                rotate([0,0,-twist * j]) {
                    for(i=[0:n_wide-1]) {
                        for(k=[0:n_long-1]) {
                            twisted_menger(twist, [-(n_wide - 1)/2 + i, -(n_long - 1)/2 + k,0]);
                        }
                    }
                }
            }
        }
    }
}

twisted_menger(twist, hole_scale, n_wide, n_long, n_tall);