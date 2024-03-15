
fn circumsphere_2d(vertices: [[f64; 2]; 3]) -> ([f64; 2], f64) {
    let [x1, y1] = vertices[0];
    let [x2, y2] = vertices[1];
    let [x3, y3] = vertices[2];

    let D = [(x1 + x2) / 2.0, (y1 + y2) / 2.0];
    let E = [(x2 + x3) / 2.0, (y2 + y3) / 2.0];

    let mAB = (y2 - y1) / (x2 - x1);
    let mBC = (y3 - y2) / (x3 - x2);

    let mD = -1. / mAB;
    let mE = -1. / mBC;

    let X = (mD * D[0] - mE * E[0] + E[1] - D[1]) / (mD - mE);
    let Y = mD * (X - D[0]) + D[1];

    let R = ((X - x1) * (X - x1) + (Y - y1) * (Y - y1)).sqrt();

    return ([X, Y], R);
}

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

pub fn circumsphere<const N: usize, const M: usize>(vertices: [[f64; N]; M]) -> ([f64; N], f64) {
    let mut center = [0.0; N];
    let mut radius = 0.0;
    if N == 2 {
        let vertices_2d = [
            [vertices[0][0], vertices[0][1]],
            [vertices[1][0], vertices[1][1]],
            [vertices[2][0], vertices[2][1]],
        ];
        let (center_2d, radius_2d) = circumsphere_2d(vertices_2d);
        for i in 0..N {
            center[i] = center_2d[i];
        }
        radius = radius_2d;
    }
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

pub fn in_sphere<const N: usize>(vertex: [f64; N], center: [f64; N], radius: f64) -> bool {
    let mut distance: f64 = 0.0;
    for i in 0..N {
        distance += (center[i] - vertex[i]) * (center[i] - vertex[i]);
    }
    distance < radius * radius
}

pub fn bounding_sphere<const N: usize>(vertices: Vec<[f64; N]>) -> ([f64; N], f64) {
    // Compute a bounding sphere for a given set of points.
    // It is not a minimal bounding sphere, but we don't need it to be minimum anyway
    let mut inside: Vec<bool> = vec![false; vertices.len()];
    let mut lower_corner: [f64; N] = [0.0; N];
    let mut upper_corner: [f64; N] = [0.0; N];
    let mut center: [f64; N] = [0.0; N];
    let mut radius: f64 = 0.0;
    for i in 0..N {
        lower_corner[i] = vertices.iter().map(|x| x[i]).fold(f64::INFINITY, f64::min);
        upper_corner[i] = vertices.iter().map(|x| x[i]).fold(f64::NEG_INFINITY, f64::max);
        center[i] = (upper_corner[i] + lower_corner[i])/2.0;
    }
    let upper_distance = upper_corner.iter().zip(center.iter()).map(|(x, y)| (x - y)*(x - y)).sum::<f64>().sqrt();
    let lower_distance = lower_corner.iter().zip(center.iter()).map(|(x, y)| (x - y)*(x - y)).sum::<f64>().sqrt();
    radius = if upper_distance > lower_distance {upper_distance} else {lower_distance};
    // Check if all points are inside the sphere.
    for (id, point) in vertices.iter().enumerate() {
        if inside[id] == false{
            if in_sphere(*point, center, radius) {
                inside[id] = true;
            }
            else {
                radius = radius *1.5;
            }
        }
    }
    (center, radius)
}
