use std::ops::{Add, AddAssign, Sub, SubAssign, DivAssign, MulAssign};
use std::cmp::Ordering;
use std::fmt;

#[derive(Copy, Clone)]
pub struct Oklab { pub l: f32, pub a: f32, pub b: f32 }

impl Oklab {
    pub fn to_linear_srgb(&self) -> LinearSrgb {
        let l_ = self.l + 0.3963377774f32*self.a + 0.2158037573f32*self.b;
        let m_ = self.l - 0.1055613458f32*self.a - 0.0638541728f32*self.b;
        let s_ = self.l - 0.0894841775f32*self.a - 1.2914855480f32*self.b;

        let l = l_*l_*l_;
        let m = m_*m_*m_;
        let s = s_*s_*s_;

        return LinearSrgb {
            r: 4.0767416621f32*l - 3.3077115913f32*m + 0.2309699292f32*s,
            g: -1.2684380046f32*l + 2.6097574011f32*m - 0.3413193965f32*s,
            b: -0.0041960863f32*l - 0.7034186147f32*m + 1.7076147010f32*s,
        };
    }

    // put into doubles. panics if colors is empty
    pub fn centroid(colors: &[Self]) -> Self {
        assert!(colors.len() > 0);
        let c0 = colors[0];
        let mut l = c0.l as f64; let mut a = c0.a as f64;
        let mut b = c0.b as f64;
        for i in 1..colors.len() {
            let c = colors[i];
            l += c.l as f64; a += c.a as f64; b += c.b as f64;
        }
        let lenf = colors.len() as f64;
        l /= lenf; a /= lenf; b /= lenf;
        Self { l: l as f32, a: a as f32, b: b as f32 }
    }

    pub fn sq_dist(&self, other: &Self) -> f32 {
        let l_d = other.l - self.l;
        let a_d = other.a - self.a;
        let b_d = other.b - self.b;
        l_d*l_d + a_d*a_d + b_d*b_d
    }

    pub const WHITE: Oklab = Oklab { l: 1.0, a: 0.0, b: 0.0 };
    pub const BLACK: Oklab = Oklab { l: 0.0, a: 0.0, b: 0.0 };
}

impl Add for Oklab {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            l: self.l + other.l,
            a: self.a + other.a,
            b: self.b + other.b,
        }
    }
}

impl AddAssign for Oklab {
    fn add_assign(&mut self, other: Self) {
        self.l += other.l;
        self.a += other.a;
        self.b += other.b;
    }
}

impl Sub for Oklab {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            l: self.l - other.l,
            a: self.a - other.a,
            b: self.b - other.b,
        }
    }
}

impl SubAssign for Oklab {
    fn sub_assign(&mut self, other: Self) {
        self.l -= other.l;
        self.a -= other.a;
        self.b -= other.b;
    }
}

impl MulAssign<f32> for Oklab {
    fn mul_assign(&mut self, rhs: f32) {
        self.l *= rhs;
        self.a *= rhs;
        self.b *= rhs;
    }
}

// todo: mulassign and div

impl DivAssign<f32> for Oklab {
    fn div_assign(&mut self, rhs: f32) {
        self.l /= rhs;
        self.a /= rhs;
        self.b /= rhs;
    }
}

impl Ord for Oklab {
    fn cmp(&self, other: &Self) -> Ordering {
        self.l.total_cmp(&other.l)
            .then(self.a.total_cmp(&other.a))
            .then(self.b.total_cmp(&other.b))
    }
}

impl PartialOrd for Oklab {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Oklab {
    fn eq(&self, other: &Self) -> bool {
        matches!(self.l.total_cmp(&other.l), Ordering::Equal)
            && matches!(self.a.total_cmp(&other.a), Ordering::Equal)
            && matches!(self.b.total_cmp(&other.b), Ordering::Equal)
    }
}

impl Eq for Oklab {}

impl fmt::Debug for Oklab {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:.2}%,{:.4},{:.4})", self.l*100.0, self.a, self.b)
    }
}

// TODO: make sure we're doing gamma the right way round.




#[derive(Debug, Copy, Clone, PartialEq)]
pub struct LinearSrgb { pub r: f32, pub g: f32, pub b: f32 }

impl LinearSrgb {
    fn delinearise(x: f32) -> f32 {
        if x >= 0.0031308 {
            1.055 * x.powf(1.0/2.4) - 0.055
        } else {
            12.92 * x
        }
    }

    pub fn to_gamma(&self) -> Srgb {
        Srgb {
            r: Self::delinearise(self.r),
            g: Self::delinearise(self.g),
            b: Self::delinearise(self.b),
        }
    }

    // https://bottosson.github.io/posts/oklab/#converting-from-linear-srgb-to-oklab
    pub fn to_oklab(&self) -> Oklab {
        // oklabify.
        let c = self;
        let l = 0.4122214708f32*c.r + 0.5363325363f32*c.g + 0.0514459929f32*c.b;
        let m = 0.2119034982f32*c.r + 0.6806995451f32*c.g + 0.1073969566f32*c.b;
        let s = 0.0883024619f32*c.r + 0.2817188376f32*c.g + 0.6299787005f32*c.b;

        let l_ = l.cbrt();
        let m_ = m.cbrt();
        let s_ = s.cbrt();

        return Oklab {
            l: 0.2104542553f32*l_ + 0.7936177850f32*m_ - 0.0040720468f32*s_,
            a: 1.9779984951f32*l_ - 2.4285922050f32*m_ + 0.4505937099f32*s_,
            b: 0.0259040371f32*l_ + 0.7827717662f32*m_ - 0.8086757660f32*s_,
        };
    }
}

/// Gamma-corrected sRGB.
// todo: derive eq where NaN is just somewhere?
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Srgb { pub r: f32, pub g: f32, pub b: f32 }

impl Srgb {
    fn linearise(x: f32) -> f32 {
        if x >= 0.04045 {
            ((x + 0.055)/(1.0 + 0.055)).powf(2.4)
        } else {
            x / 12.92
        }
    }

    pub fn to_linear(&self) -> LinearSrgb {
        LinearSrgb {
            r: Self::linearise(self.r),
            g: Self::linearise(self.g),
            b: Self::linearise(self.b),
        }
    }

    fn discretise_5bit(x: f32) -> u16 {
        (x * 31.0).clamp(0.0, 31.0).round() as u16
    }

    fn discretise_8bit(x: f32) -> u16 {
        (x * 255.0).clamp(0.0, 255.0).round() as u16
    }

    pub fn to_srgb555(&self) -> Srgb555 {
        let r_ = Self::discretise_5bit(self.r);
        let g_ = Self::discretise_5bit(self.g);
        let b_ = Self::discretise_5bit(self.b);
        let u = (r_ << 10) | (g_ << 5) | b_;
        Srgb555(u)
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct Srgb555(pub u16);

impl Srgb555 {
    pub fn to_srgb(self) -> Srgb {
        let r_ = (self.0 & 0b0_11111_00000_00000) >> 10;
        let g_ = (self.0 & 0b0_00000_11111_00000) >> 5;
        let b_ = self.0 & 0b0_00000_00000_11111;
        let r = r_ as f32 / 31.0;
        let g = g_ as f32 / 31.0;
        let b = b_ as f32 / 31.0;
        Srgb { r, g, b }
    }
}








// Bayer matrices. I've added 1/(2n^2) to each one so that the values
// are equally spaced within [0,1].

const BAYER_2: [[f32; 2]; 2] = [
    [0.125, 0.625],
    [0.875, 0.375],
];

const BAYER_4: [[f32; 4]; 4] = [
    [0.03125, 0.53125, 0.15625, 0.65625],
    [0.78125, 0.28125, 0.90625, 0.40625],
    [0.21875, 0.71875, 0.09375, 0.59375],
    [0.96875, 0.46875, 0.84375, 0.34375],
];

const BAYER_8: [[f32; 8]; 8] = [
    [0.0078125, 0.5078125, 0.1328125, 0.6328125, 0.0390625, 0.5390625, 0.1640625, 0.6640625],
    [0.7578125, 0.2578125, 0.8828125, 0.3828125, 0.7890625, 0.2890625, 0.9140625, 0.4140625],
    [0.1953125, 0.6953125, 0.0703125, 0.5703125, 0.2265625, 0.7265625, 0.1015625, 0.6015625],
    [0.9453125, 0.4453125, 0.8203125, 0.3203125, 0.9765625, 0.4765625, 0.8515625, 0.3515625],
    [0.0546875, 0.5546875, 0.1796875, 0.6796875, 0.0234375, 0.5234375, 0.1484375, 0.6484375],
    [0.8046875, 0.3046875, 0.9296875, 0.4296875, 0.7734375, 0.2734375, 0.8984375, 0.3984375],
    [0.2421875, 0.7421875, 0.1171875, 0.6171875, 0.2109375, 0.7109375, 0.0859375, 0.5859375],
    [0.9921875, 0.4921875, 0.8671875, 0.3671875, 0.9609375, 0.4609375, 0.8359375, 0.3359375],
];

const BAYER_16: [[f32; 16]; 16] = [
    [0.001953125, 0.501953125, 0.126953125, 0.626953125,
    0.033203125, 0.533203125, 0.158203125, 0.658203125,
    0.009765625, 0.509765625, 0.134765625, 0.634765625,
    0.041015625, 0.541015625, 0.166015625, 0.666015625],
    [0.751953125, 0.251953125, 0.876953125, 0.376953125,
    0.783203125, 0.283203125, 0.908203125, 0.408203125,
    0.759765625, 0.259765625, 0.884765625, 0.384765625,
    0.791015625, 0.291015625, 0.916015625, 0.416015625],
    [0.189453125, 0.689453125, 0.064453125, 0.564453125,
    0.220703125, 0.720703125, 0.095703125, 0.595703125,
    0.197265625, 0.697265625, 0.072265625, 0.572265625,
    0.228515625, 0.728515625, 0.103515625, 0.603515625],
    [0.939453125, 0.439453125, 0.814453125, 0.314453125,
    0.970703125, 0.470703125, 0.845703125, 0.345703125,
    0.947265625, 0.447265625, 0.822265625, 0.322265625,
    0.978515625, 0.478515625, 0.853515625, 0.353515625],
    [0.048828125, 0.548828125, 0.173828125, 0.673828125, 0.017578125, 0.517578125, 0.142578125, 0.642578125, 0.056640625, 0.556640625, 0.181640625, 0.681640625, 0.025390625, 0.525390625, 0.150390625, 0.650390625], [0.798828125, 0.298828125, 0.923828125, 0.423828125, 0.767578125, 0.267578125, 0.892578125, 0.392578125, 0.806640625, 0.306640625, 0.931640625, 0.431640625, 0.775390625, 0.275390625, 0.900390625, 0.400390625], [0.236328125, 0.736328125, 0.111328125, 0.611328125, 0.205078125, 0.705078125, 0.080078125, 0.580078125, 0.244140625, 0.744140625, 0.119140625, 0.619140625, 0.212890625, 0.712890625, 0.087890625, 0.587890625], [0.986328125, 0.486328125, 0.861328125, 0.361328125, 0.955078125, 0.455078125, 0.830078125, 0.330078125, 0.994140625, 0.494140625, 0.869140625, 0.369140625, 0.962890625, 0.462890625, 0.837890625, 0.337890625], [0.013671875, 0.513671875, 0.138671875, 0.638671875, 0.044921875, 0.544921875, 0.169921875, 0.669921875, 0.005859375, 0.505859375, 0.130859375, 0.630859375, 0.037109375, 0.537109375, 0.162109375, 0.662109375], [0.763671875, 0.263671875, 0.888671875, 0.388671875, 0.794921875, 0.294921875, 0.919921875, 0.419921875, 0.755859375, 0.255859375, 0.880859375, 0.380859375, 0.787109375, 0.287109375, 0.912109375, 0.412109375], [0.201171875, 0.701171875, 0.076171875, 0.576171875, 0.232421875, 0.732421875, 0.107421875, 0.607421875, 0.193359375, 0.693359375, 0.068359375, 0.568359375, 0.224609375, 0.724609375, 0.099609375, 0.599609375], [0.951171875, 0.451171875, 0.826171875, 0.326171875, 0.982421875, 0.482421875, 0.857421875, 0.357421875, 0.943359375, 0.443359375, 0.818359375, 0.318359375, 0.974609375, 0.474609375, 0.849609375, 0.349609375], [0.060546875, 0.560546875, 0.185546875, 0.685546875, 0.029296875, 0.529296875, 0.154296875, 0.654296875, 0.052734375, 0.552734375, 0.177734375, 0.677734375, 0.021484375, 0.521484375, 0.146484375, 0.646484375], [0.810546875, 0.310546875, 0.935546875, 0.435546875, 0.779296875, 0.279296875, 0.904296875, 0.404296875, 0.802734375, 0.302734375, 0.927734375, 0.427734375, 0.771484375, 0.271484375, 0.896484375, 0.396484375], [0.248046875, 0.748046875, 0.123046875, 0.623046875, 0.216796875, 0.716796875, 0.091796875, 0.591796875, 0.240234375, 0.740234375, 0.115234375, 0.615234375, 0.208984375, 0.708984375, 0.083984375, 0.583984375], [0.998046875, 0.498046875, 0.873046875, 0.373046875, 0.966796875, 0.466796875, 0.841796875, 0.341796875, 0.990234375, 0.490234375, 0.865234375, 0.365234375, 0.958984375, 0.458984375, 0.833984375, 0.333984375]];

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn oklab_impls() {
        let c1 = Oklab { l: 1.0, a: 0.5, b: 0.75 };
        let c2 = Oklab { l: 1.0, a: 0.5, b: 0.75 };
        let c3 = Oklab { l: 1.0, a: 0.6, b: 0.75 };
        assert!(c1 == c2);
        assert!(c2 <= c3);
        assert!(c2 < c3);
    }
}
