fn circumsphere_3d(vertices: [[f64; 3]; 4]) -> ([f64; 3], f64) {
    let a = [
        vertices[1][0] - vertices[0][0],
        vertices[1][1] - vertices[0][1],
        vertices[1][2] - vertices[0][2],
    ];
    let b = [
        vertices[2][0] - vertices[0][0],
        vertices[2][1] - vertices[0][1],
        vertices[2][2] - vertices[0][2],
    ];
    let c = [
        vertices[3][0] - vertices[0][0],
        vertices[3][1] - vertices[0][1],
        vertices[3][2] - vertices[0][2],
    ];

    let length_a = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    let length_b = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
    let length_c = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];

    let a_cross_b = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let b_cross_c = [
        b[1] * c[2] - b[2] * c[1],
        b[2] * c[0] - b[0] * c[2],
        b[0] * c[1] - b[1] * c[0],
    ];
    let c_cross_a = [
        c[1] * a[2] - c[2] * a[1],
        c[2] * a[0] - c[0] * a[2],
        c[0] * a[1] - c[1] * a[0],
    ];

    let denominator = 0.5 / (a[0] * b_cross_c[0] + a[1] * b_cross_c[1] + a[2] * b_cross_c[2]);

    let center = [
        (length_a * b_cross_c[0] + length_b * c_cross_a[0] + length_c * a_cross_b[0]) * denominator,
        (length_a * b_cross_c[1] + length_b * c_cross_a[1] + length_c * a_cross_b[1]) * denominator,
        (length_a * b_cross_c[2] + length_b * c_cross_a[2] + length_c * a_cross_b[2]) * denominator,
    ];

    let radius = (center[0] - vertices[0][0]) * (center[0] - vertices[0][0])
        + (center[1] - vertices[0][1]) * (center[1] - vertices[0][1])
        + (center[2] - vertices[0][2]) * (center[2] - vertices[0][2]);

    (center, radius.sqrt())
}

pub fn circumsphere<const N:usize, const M:usize>(vertices: [[f64; N]; M]) -> ([f64; N], f64) {
    let mut center = [0.0; N];
    let mut radius = 0.0;
    if N == 3 {
        let vertices_3d = [
            [vertices[0][0], vertices[0][1], vertices[0][2]],
            [vertices[1][0], vertices[1][1], vertices[1][2]],
            [vertices[2][0], vertices[2][1], vertices[2][2]],
            [vertices[3][0], vertices[3][1], vertices[3][2]],
        ];
        let (center_3d, radius_3d) = circumsphere_3d(vertices_3d);
        for i in 0..N {
            center[i] = center_3d[i];
        }
        radius = radius_3d;
    }
    (center, radius)
}

pub fn in_sphere<const N:usize>(vertex: [f64; N], center: [f64; N], radius: f64) -> bool {
    let mut distance: f64 = 0.0;
    for i in 0..N {
        distance += (center[i] - vertex[i]) * (center[i] - vertex[i]);
    }
    distance < radius * radius
}