// around_point.rs
use crate::types::Float;

pub mod types;

#[derive(Debug, Clone)]
pub struct AroundPoint {
    /// For each ray, we write 3 outputs (RGB). Length must be `angles.len()` * 3.
    pub values: Box<[Float]>,
    /// Ray angles in DEGREES, measured counterclockwise from +Y (forward).
    pub angles: Box<[Float]>,
    /// Max ray length used for intersection tests.
    pub length: Float,
}

/// Shapes in the *boat-local* frame where +Y is forward and the origin is the boat position.
/// We include color so `AroundPoint` can write RGBs directly.
pub enum TargetShape {
    Circle {
        center: (Float, Float),
        radius: Float,
    },
}
pub struct RGB {
    pub red: Float,
    pub green: Float,
    pub blue: Float,
}

impl AroundPoint {
    #[inline]
    fn compute_t_hit(
        proj: Float,
        perp: Float,
        r2: Float,
        radius: Float,
        length: Float,
    ) -> Option<Float> {
        // early-outs along the ray axis
        if proj < -radius || proj > length + radius {
            return None;
        }
        let d2 = perp * perp;
        if d2 > r2 {
            return None;
        }
        let thc = (r2 - d2).sqrt();
        let t0 = proj - thc;
        let t1 = proj + thc;
        if t1 < 0.0 || t0 > length {
            return None;
        }
        let t_hit = if t0 <= 0.0 { 0.0 } else { t0.min(length) };
        Some(t_hit)
    }

    /// `origin` is unused by the math (targets are already boat-local), but kept for API symmetry.
    /// `targets` must be pre-transformed into the boat-local frame (+Y forward).
    pub fn update(&mut self, _origin: (Float, Float), targets: &[(TargetShape, RGB)]) {
        // Ensure values buffer size matches ray count * 3 (RGB)
        let num_rays = self.angles.len();
        let num_outputs = 3usize;
        if self.values.len() != num_rays * num_outputs {
            self.values = vec![0.0; num_rays * num_outputs].into_boxed_slice();
        }
        // Clear previous frame
        for v in &mut self.values {
            *v = 0.0;
        }

        if num_rays == 0 || targets.is_empty() {
            return;
        }

        // Precompute sin/cos for each angle (angles are in degrees)
        let mut sins = vec![0.0 as Float; num_rays];
        let mut coss = vec![0.0 as Float; num_rays];
        for (i, &deg) in self.angles.iter().enumerate() {
            let r = deg.to_radians();
            sins[i] = r.sin();
            coss[i] = r.cos();
        }

        // Broad-phase cone bound: max |theta|
        let theta_max_abs = self
            .angles
            .iter()
            .fold(0.0_f64 as Float, |m, &a| m.max(a.abs()));
        let tan_theta_max = theta_max_abs.to_radians().tan();
        let length = self.length;
        let max_d2_pad = (length) * (length); // we'll add per-target radius below

        // Track best t per ray
        let mut best_t = vec![Float::INFINITY; num_rays];

        for (tgt, rgb) in targets {
            let TargetShape::Circle { center, radius } = *tgt;

            let (x, y) = center;
            let r2 = radius * radius;

            // Broad-phase culls (boat-local frame):
            // 1) Behind the origin more than radius
            if y <= -radius {
                continue;
            }
            // 2) Outside max range + radius
            let d2 = x.mul_add(x, y * y);
            if d2 >= (2.0 * length).mul_add(radius, radius.mul_add(radius, max_d2_pad)) {
                continue;
            }
            // 3) Outside the cone + radius
            if x.abs() > (y + radius).mul_add(tan_theta_max, radius) {
                continue;
            }

            // Per-ray precise tests (ONLY thetas in self.angles; no ± pairing)
            for j in 0..num_rays {
                // Ray oriented by +theta (CCW from +Y)
                // Transform target center into the ray's local coords:
                // proj = -x*sin(theta) + y*cos(theta)
                // perp =  x*cos(theta) + y*sin(theta)
                let s = sins[j];
                let c = coss[j];
                let proj = (-x).mul_add(s, y * c);
                let perp = x.mul_add(c, y * s);

                if let Some(t) = Self::compute_t_hit(proj, perp, r2, radius, length) {
                    if t < best_t[j] {
                        best_t[j] = t;
                        let base = j * num_outputs;
                        // Write RGB for the closest hit so far on this ray
                        self.values[base] = rgb.red;
                        self.values[base + 1] = rgb.blue;
                        self.values[base + 2] = rgb.green;
                    }
                }
            }
        }
    }
}
