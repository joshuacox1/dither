use std::ops::{Add, AddAssign, Sub, SubAssign, Div, DivAssign, Mul, MulAssign};
use std::cmp::Ordering;
use std::fmt;
use std::fmt::Write;
use std::borrow::Cow;
use std::hash::Hash;
use std::collections::HashSet;

/// A colour in the Oklab colour space. Uses the revised lightness
/// estimate L_r described in <https://bottosson.github.io/posts/colorpicker/#intermission---a-new-lightness-estimate-for-oklab>,
/// which significantly improves lightness prediction.
/// TODO: rename Oklab and store l_r internally. Various discussions.
#[derive(Debug, Copy, Clone)]
pub struct Oklabr {
    pub l: f32,
    pub a: f32,
    pub b: f32,
}

// TODO: for good style, make add etc. traits take references.

const K1: f32 = 0.206;
const K2: f32 = 0.03;
const K3: f32 = (1.0 + K1) / (1.0 + K2);

impl Oklabr {
    /// The original `l` component of an `Oklab` colour. Not the "revised lightness
    /// estimate" described in [link].
    pub fn l(&self) -> f32 { Self::new_to_old_lightness(self.l) }

    /// Converts an `Oklabr` colour to an RGB888 colour, clamping
    /// each RGB channel if the colour is out of gamut. For best results
    /// it is recommended to call `srgb_gamut_map` on the colour first
    /// to obtain a suitably gamut-mapped `Oklabr` colour.
    pub fn to_rgb888(&self) -> Rgb888 {
        self.to_rgbf32().to_rgb888()
    }

    /// Convert to gamma-encoded RGB in the RgbF32 colour space.
    fn to_rgbf32(&self) -> RgbF32 {
        let old_l = Self::new_to_old_lightness(self.l);

        let l_ = old_l + 0.3963377774f32*self.a + 0.2158037573f32*self.b;
        let m_ = old_l - 0.1055613458f32*self.a - 0.0638541728f32*self.b;
        let s_ = old_l - 0.0894841775f32*self.a - 1.2914855480f32*self.b;

        let l = l_*l_*l_;
        let m = m_*m_*m_;
        let s = s_*s_*s_;

        let r_ = 4.0767416621f32*l - 3.3077115913f32*m + 0.2309699292f32*s;
        let g_ = -1.2684380046f32*l + 2.6097574011f32*m - 0.3413193965f32*s;
        let b_ = -0.0041960863f32*l - 0.7034186147f32*m + 1.7076147010f32*s;

        RgbF32 {
            r: RgbF32::delinearise(r_),
            g: RgbF32::delinearise(g_),
            b: RgbF32::delinearise(b_),
        }
    }

    fn old_to_new_lightness(l: f32) -> f32 {
        let x = K3*l - K1;
        0.5*(x + (x*x + 4.0*K2*K3*l).sqrt())
    }

    fn new_to_old_lightness(l_r: f32) -> f32 {
        let top = l_r * (l_r + K1);
        let bottom = K3 * (l_r + K2);
        top / bottom
    }

    /// Squared distance between this Oklabr colour and another.
    pub fn sq_dist(&self, other: &Self) -> f32 {
        let l_d = other.l - self.l;
        let a_d = other.a - self.a;
        let b_d = other.b - self.b;
        l_d*l_d + a_d*a_d + b_d*b_d
    }

    pub const WHITE: Oklabr = Oklabr { l: 1.0, a: 0.0, b: 0.0 };
    pub const BLACK: Oklabr = Oklabr { l: 0.0, a: 0.0, b: 0.0 };

    /// Returns whichever of white or black has the greater
    /// contrast with the given colour. One use of this is determining
    /// whether text placed over this colour should be black or white.
    pub fn white_or_black_most_contrast(&self) -> &'static Self {
        // Contrast is defined by the W3C as (L1 + 0.05) / (L2 + 0.05),
        // where L1 is the lightness of the lighter of the two colours
        // and L2 of the darker colour. Since white is lighter than
        // every colour and black darker than every colour, we want
        // to compare A = 0.05 / (l + 0.05) and B = (l + 0.05) / 1.05.
        // We have A <= B when 1.05/0.05 <= (l + 0.05)^2. Simplifying
        // we obtain the cut-off point of sqrt(0.0525)-0.05 = 0.17912...
        // Below this number the contrast with white is higher;
        // Above it, the contrast with black is higher.
        let cutoff = 0.0525f32.sqrt() - 0.05f32;
        if self.l < cutoff { &Self::WHITE } else { &Self::BLACK }
    }

    /// If the colour is within the RgbF32 gamut, returns a reference
    /// to it. Otherwise computes a suitably close colour within the
    /// RgbF32 gamut and returns that (owned). TODO: link the gamut
    /// mapping article. This is a highly complex operation.
    pub fn srgb_gamut_map(&self) -> Cow<'_, Self> {
        // todo: make this do something.
        Cow::Borrowed(self)
    }

    /// Returns whether an `Oklabr` colour is "well-formed", meaning
    /// none of the entries are infinities or NaN.
    pub fn is_well_formed(&self) -> bool {
        self.l.is_finite() && self.a.is_finite() && self.b.is_finite()
    }

    /// Returns whether an `Oklabr` colour is in the `Oklabr` gamut.
    /// This is true when all three components are non-infinite-or-NaN
    /// and `l` lies in `[0.0, 1.0]`.
    pub fn is_in_gamut(&self) -> bool {
        self.is_well_formed() && 0.0 <= self.l && self.l <= 1.0
    }
}

// TODO: hash?

impl Add for Oklabr {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            l: self.l + other.l,
            a: self.a + other.a,
            b: self.b + other.b,
        }
    }
}

impl Add<&Oklabr> for Oklabr {
    type Output = Self;

    fn add(self, other: &Oklabr) -> Self {
        Self {
            l: self.l + other.l,
            a: self.a + other.a,
            b: self.b + other.b,
        }
    }
}

impl Add<Oklabr> for &Oklabr {
    type Output = Oklabr;

    fn add(self, other: Oklabr) -> Oklabr {
        Oklabr {
            l: self.l + other.l,
            a: self.a + other.a,
            b: self.b + other.b,
        }
    }
}

impl Add for &Oklabr {
    type Output = Oklabr;

    fn add(self, other: Self) -> Oklabr {
        Oklabr {
            l: self.l + other.l,
            a: self.a + other.a,
            b: self.b + other.b,
        }
    }
}

impl AddAssign for Oklabr {
    fn add_assign(&mut self, other: Self) {
        self.l += other.l;
        self.a += other.a;
        self.b += other.b;
    }
}

impl AddAssign<&Oklabr> for Oklabr {
    fn add_assign(&mut self, other: &Oklabr) {
        self.l += other.l;
        self.a += other.a;
        self.b += other.b;
    }
}

impl Sub for Oklabr {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            l: self.l - other.l,
            a: self.a - other.a,
            b: self.b - other.b,
        }
    }
}

impl Sub<&Oklabr> for Oklabr {
    type Output = Self;

    fn sub(self, other: &Oklabr) -> Self {
        Self {
            l: self.l - other.l,
            a: self.a - other.a,
            b: self.b - other.b,
        }
    }
}

impl Sub<Oklabr> for &Oklabr {
    type Output = Oklabr;

    fn sub(self, other: Oklabr) -> Oklabr {
        Oklabr {
            l: self.l - other.l,
            a: self.a - other.a,
            b: self.b - other.b,
        }
    }
}

impl Sub for &Oklabr {
    type Output = Oklabr;

    fn sub(self, other: Self) -> Oklabr {
        Oklabr {
            l: self.l - other.l,
            a: self.a - other.a,
            b: self.b - other.b,
        }
    }
}

impl SubAssign for Oklabr {
    fn sub_assign(&mut self, other: Self) {
        self.l -= other.l;
        self.a -= other.a;
        self.b -= other.b;
    }
}

impl SubAssign<&Oklabr> for Oklabr {
    fn sub_assign(&mut self, other: &Oklabr) {
        self.l -= other.l;
        self.a -= other.a;
        self.b -= other.b;
    }
}

impl Mul<f32> for Oklabr {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self {
        Self {
            l: self.l * rhs,
            a: self.a * rhs,
            b: self.b * rhs,
        }
    }
}

impl Mul<&f32> for Oklabr {
    type Output = Self;

    fn mul(self, rhs: &f32) -> Self {
        Self {
            l: self.l * *rhs,
            a: self.a * *rhs,
            b: self.b * *rhs,
        }
    }
}

impl Mul<f32> for &Oklabr {
    type Output = Oklabr;

    fn mul(self, rhs: f32) -> Oklabr {
        Oklabr {
            l: self.l * rhs,
            a: self.a * rhs,
            b: self.b * rhs,
        }
    }
}

impl Mul<&f32> for &Oklabr {
    type Output = Oklabr;

    fn mul(self, rhs: &f32) -> Oklabr {
        Oklabr {
            l: self.l * *rhs,
            a: self.a * *rhs,
            b: self.b * *rhs,
        }
    }
}

impl MulAssign<f32> for Oklabr {
    fn mul_assign(&mut self, rhs: f32) {
        self.l *= rhs;
        self.a *= rhs;
        self.b *= rhs;
    }
}

impl MulAssign<&f32> for Oklabr {
    fn mul_assign(&mut self, rhs: &f32) {
        self.l *= *rhs;
        self.a *= *rhs;
        self.b *= *rhs;
    }
}

impl Div<f32> for Oklabr {
    type Output = Self;

    fn div(self, rhs: f32) -> Self {
        Self {
            l: self.l / rhs,
            a: self.a / rhs,
            b: self.b / rhs,
        }
    }
}

impl Div<&f32> for Oklabr {
    type Output = Self;

    fn div(self, rhs: &f32) -> Self {
        Self {
            l: self.l / *rhs,
            a: self.a / *rhs,
            b: self.b / *rhs,
        }
    }
}

impl Div<f32> for &Oklabr {
    type Output = Oklabr;

    fn div(self, rhs: f32) -> Oklabr {
        Oklabr {
            l: self.l / rhs,
            a: self.a / rhs,
            b: self.b / rhs,
        }
    }
}

impl Div<&f32> for &Oklabr {
    type Output = Oklabr;

    fn div(self, rhs: &f32) -> Oklabr {
        Oklabr {
            l: self.l / *rhs,
            a: self.a / *rhs,
            b: self.b / *rhs,
        }
    }
}

impl DivAssign<f32> for Oklabr {
    fn div_assign(&mut self, rhs: f32) {
        self.l /= rhs;
        self.a /= rhs;
        self.b /= rhs;
    }
}

impl DivAssign<&f32> for Oklabr {
    fn div_assign(&mut self, rhs: &f32) {
        self.l /= *rhs;
        self.a /= *rhs;
        self.b /= *rhs;
    }
}

impl Ord for Oklabr {
    fn cmp(&self, other: &Self) -> Ordering {
        self.l.total_cmp(&other.l)
            .then(self.a.total_cmp(&other.a))
            .then(self.b.total_cmp(&other.b))
    }
}

impl PartialOrd for Oklabr {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Oklabr {
    fn eq(&self, other: &Self) -> bool {
        matches!(self.l.total_cmp(&other.l), Ordering::Equal)
            && matches!(self.a.total_cmp(&other.a), Ordering::Equal)
            && matches!(self.b.total_cmp(&other.b), Ordering::Equal)
    }
}

impl Eq for Oklabr { }

impl fmt::Display for Oklabr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // TODO: should this show old or revised lightness
        write!(f, "({:.2}%,{:.4},{:.4})", self.l*100.0, self.a, self.b)
    }
}

// // Finds the maximum saturation possible for a given hue that fits in RgbF32
// // Saturation here is defined as S = C/L
// // a and b must be normalized so a^2 + b^2 == 1
// fn compute_max_saturation(a: f32, b: f32)
// {
//     // Max saturation will be when one of r, g or b goes below zero.

//     // Select different coefficients depending on which component goes below zero first
//     let k0, k1, k2, k3, k4, wl, wm, ws;

//     if (-1.88170328 * a - 0.80936493 * b > 1)
//     {
//         // Red component
//         k0 = 1.19086277; k1 = 1.76576728; k2 = 0.59662641;
//         k3 = 0.75515197; k4 = 0.56771245;
//         wl = 4.0767416621; wm = -3.3077115913; ws = 0.2309699292;
//     }
//     else if (1.81444104 * a - 1.19445276 * b > 1)
//     {
//         // Green component
//         k0 = 0.73956515; k1 = -0.45954404; k2 = 0.08285427;
//         k3 = 0.12541070; k4 = 0.14503204;
//         wl = -1.2684380046; wm = 2.6097574011; ws = -0.3413193965;
//     }
//     else
//     {
//         // Blue component
//         k0 = 1.35733652; k1 = -0.00915799; k2 = -1.15130210;
//         k3 = -0.50559606; k4 = 0.00692167;
//         wl = -0.0041960863; wm = -0.7034186147; ws = 1.7076147010;
//     }

//     // Approximate max saturation using a polynomial:
//     let S = k0 + k1 * a + k2 * b + k3 * a * a + k4 * a * b;

//     // Do one step Halley's method to get closer
//     // this gives an error less than 10e6, except for some blue hues where the dS/dh is close to infinite
//     // this should be sufficient for most applications, otherwise do two/three steps

//     let k_l = 0.3963377774 * a + 0.2158037573 * b;
//     let k_m = -0.1055613458 * a - 0.0638541728 * b;
//     let k_s = -0.0894841775 * a - 1.2914855480 * b;

//     {
//         let l_ = 1.0 + S * k_l;
//         let m_ = 1.0 + S * k_m;
//         let s_ = 1.0 + S * k_s;

//         let l = l_ * l_ * l_;
//         let m = m_ * m_ * m_;
//         let s = s_ * s_ * s_;

//         let l_dS = 3.0 * k_l * l_ * l_;
//         let m_dS = 3.0 * k_m * m_ * m_;
//         let s_dS = 3.0 * k_s * s_ * s_;

//         let l_dS2 = 6.0 * k_l * k_l * l_;
//         let m_dS2 = 6.0 * k_m * k_m * m_;
//         let s_dS2 = 6.0 * k_s * k_s * s_;

//         let f  = wl * l     + wm * m     + ws * s;
//         let f1 = wl * l_dS  + wm * m_dS  + ws * s_dS;
//         let f2 = wl * l_dS2 + wm * m_dS2 + ws * s_dS2;

//         S = S - f * f1 / (f1*f1 - 0.5 * f * f2);
//     }

//     return S;
// }

// // finds L_cusp and C_cusp for a given hue
// // a and b must be normalized so a^2 + b^2 == 1
// struct LC { float L; float C; };
// fn find_cusp(a: f32, b: f32) -> (f32, f32) {
//     // First, find the maximum saturation (saturation S = C/L)
//     let S_cusp = compute_max_saturation(a, b);

//     // Convert to linear RgbF32 to find the first point where at least one of r,g or b >= 1:
//     RGB rgb_at_max = oklab_to_linear_RgbF32(Oklabr { 1.0, S_cusp * a, S_cusp * b });
//     float L_cusp = cbrtf(1.0 / max(max(rgb_at_max.r, rgb_at_max.g), rgb_at_max.b));
//     float C_cusp = L_cusp * S_cusp;

//     return { L_cusp , C_cusp };
// }

// // Finds intersection of the line defined by
// // L = L0 * (1 - t) + t * L1;
// // C = t * C1;
// // a and b must be normalized so a^2 + b^2 == 1
// float find_gamut_intersection(float a, float b, float L1, float C1, float L0)
// {
//     // Find the cusp of the gamut triangle
//     LC cusp = find_cusp(a, b);

//     // Find the intersection for upper and lower half seprately
//     float t;
//     if (((L1 - L0) * cusp.C - (cusp.L - L0) * C1) <= 0.f)
//     {
//         // Lower half

//         t = cusp.C * L0 / (C1 * cusp.L + cusp.C * (L0 - L1));
//     }
//     else
//     {
//         // Upper half

//         // First intersect with triangle
//         t = cusp.C * (L0 - 1.f) / (C1 * (cusp.L - 1.f) + cusp.C * (L0 - L1));

//         // Then one step Halley's method
//         {
//             float dL = L1 - L0;
//             float dC = C1;

//             float k_l = +0.3963377774f * a + 0.2158037573f * b;
//             float k_m = -0.1055613458f * a - 0.0638541728f * b;
//             float k_s = -0.0894841775f * a - 1.2914855480f * b;

//             float l_dt = dL + dC * k_l;
//             float m_dt = dL + dC * k_m;
//             float s_dt = dL + dC * k_s;


//             // If higher accuracy is required, 2 or 3 iterations of the following block can be used:
//             {
//                 float L = L0 * (1.f - t) + t * L1;
//                 float C = t * C1;

//                 float l_ = L + C * k_l;
//                 float m_ = L + C * k_m;
//                 float s_ = L + C * k_s;

//                 float l = l_ * l_ * l_;
//                 float m = m_ * m_ * m_;
//                 float s = s_ * s_ * s_;

//                 float ldt = 3 * l_dt * l_ * l_;
//                 float mdt = 3 * m_dt * m_ * m_;
//                 float sdt = 3 * s_dt * s_ * s_;

//                 float ldt2 = 6 * l_dt * l_dt * l_;
//                 float mdt2 = 6 * m_dt * m_dt * m_;
//                 float sdt2 = 6 * s_dt * s_dt * s_;

//                 float r = 4.0767416621f * l - 3.3077115913f * m + 0.2309699292f * s - 1;
//                 float r1 = 4.0767416621f * ldt - 3.3077115913f * mdt + 0.2309699292f * sdt;
//                 float r2 = 4.0767416621f * ldt2 - 3.3077115913f * mdt2 + 0.2309699292f * sdt2;

//                 float u_r = r1 / (r1 * r1 - 0.5f * r * r2);
//                 float t_r = -r * u_r;

//                 float g = -1.2684380046f * l + 2.6097574011f * m - 0.3413193965f * s - 1;
//                 float g1 = -1.2684380046f * ldt + 2.6097574011f * mdt - 0.3413193965f * sdt;
//                 float g2 = -1.2684380046f * ldt2 + 2.6097574011f * mdt2 - 0.3413193965f * sdt2;

//                 float u_g = g1 / (g1 * g1 - 0.5f * g * g2);
//                 float t_g = -g * u_g;

//                 float b = -0.0041960863f * l - 0.7034186147f * m + 1.7076147010f * s - 1;
//                 float b1 = -0.0041960863f * ldt - 0.7034186147f * mdt + 1.7076147010f * sdt;
//                 float b2 = -0.0041960863f * ldt2 - 0.7034186147f * mdt2 + 1.7076147010f * sdt2;

//                 float u_b = b1 / (b1 * b1 - 0.5f * b * b2);
//                 float t_b = -b * u_b;

//                 t_r = u_r >= 0.f ? t_r : FLT_MAX;
//                 t_g = u_g >= 0.f ? t_g : FLT_MAX;
//                 t_b = u_b >= 0.f ? t_b : FLT_MAX;

//                 t += min(t_r, min(t_g, t_b));
//             }
//         }
//     }

//     return t;
// }





























/// Gamma-encoded RgbF32.
/// TODO: make this linear. Apply gamma later?
/// TODO: rename this RgbF32.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct RgbF32 { pub r: f32, pub g: f32, pub b: f32 }

impl RgbF32 {
    fn linearise(x: f32) -> f32 {
        if x >= 0.04045 {
            ((x + 0.055)/(1.0 + 0.055)).powf(2.4)
        } else {
            x / 12.92
        }
    }

    fn delinearise(x: f32) -> f32 {
        if x >= 0.0031308 {
            1.055 * x.powf(1.0/2.4) - 0.055
        } else {
            12.92 * x
        }
    }

    // https://bottosson.github.io/posts/oklab/#converting-from-linear-RgbF32-to-oklab
    pub fn to_oklabr(&self) -> Oklabr {
        // linearise
        let r = Self::linearise(self.r);
        let g = Self::linearise(self.g);
        let b = Self::linearise(self.b);

        // oklabify
        let l = 0.4122214708f32*r + 0.5363325363f32*g + 0.0514459929f32*b;
        let m = 0.2119034982f32*r + 0.6806995451f32*g + 0.1073969566f32*b;
        let s = 0.0883024619f32*r + 0.2817188376f32*g + 0.6299787005f32*b;

        let l_ = l.cbrt();
        let m_ = m.cbrt();
        let s_ = s.cbrt();

        let l = 0.2104542553f32*l_ + 0.7936177850f32*m_ - 0.0040720468f32*s_;
        let a = 1.9779984951f32*l_ - 2.4285922050f32*m_ + 0.4505937099f32*s_;
        let b = 0.0259040371f32*l_ + 0.7827717662f32*m_ - 0.8086757660f32*s_;

        // revised lightness estimate
        let l_r = Oklabr::old_to_new_lightness(l);

        Oklabr { l: l_r, a, b }
    }

    /// Converts to RGB 888, clamping within [0.0, 1.0] and rounding.
    pub fn to_rgb888(&self) -> Rgb888 {
        fn discretise_8bit(x: f32) -> u8 {
            (x * 255.0).clamp(0.0, 255.0).round() as u8
        }

        Rgb888 {
            r: discretise_8bit(self.r),
            g: discretise_8bit(self.g),
            b: discretise_8bit(self.b),
        }
    }
}

/// 24-bit RGB true colour. Intended to be interpreted in the RgbF32
/// colour space.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct Rgb888 {
    /// The red component.
    pub r: u8,
    /// The green component.
    pub g: u8,
    /// The blue component.
    pub b: u8,
}

impl Rgb888 {
    fn to_rgbf32(&self) -> RgbF32 {
        RgbF32 {
            r: self.r as f32 / 255.0,
            g: self.g as f32 / 255.0,
            b: self.b as f32 / 255.0,
        }
    }

    /// Maps an `Rgb888` interpreted in the RgbF32 colour space to its
    /// corresponding `Oklabr` colour.
    pub fn to_oklabr(&self) -> Oklabr {
        self.to_rgbf32().to_oklabr()
    }

    /// Converts an `Rgb888` into an #RRGGBB hex string.
    pub fn to_str(&self) -> String {
        let mut result = String::with_capacity(7);
        let Self { r, g, b } = self;
        write!(&mut result, "#{r:02x}{g:02x}{b:02x}").unwrap();
        result
    }

    /// Attempts to parse an #RRGGBB hex string into an `Rgb888`.
    pub fn from_str(input: &str) -> Result<Self, ()> {
        // 0xff = ERROR.
        fn from_char_opt(x: u8) -> u8 {
            match x {
                b'0' => 0, b'1' => 1, b'2' => 2, b'3' => 3, b'4' => 4,
                b'5' => 5, b'6' => 6, b'7' => 7, b'8' => 8, b'9' => 9,
                b'A' | b'a' => 0xa, b'B' | b'b' => 0xb, b'C' | b'c' => 0xc,
                b'D' | b'd' => 0xd, b'E' | b'e' => 0xe, b'F' | b'f' => 0xf,
                _ => 0xff,
            }
        }

        let input_b = input.as_bytes();
        if input_b.len() != 7 || input_b[0] != b'#' { return Err(()); }
        let r1 = from_char_opt(input_b[1]);
        let r2 = from_char_opt(input_b[2]);
        let g1 = from_char_opt(input_b[3]);
        let g2 = from_char_opt(input_b[4]);
        let b1 = from_char_opt(input_b[5]);
        let b2 = from_char_opt(input_b[6]);
        if r1 == 0xff || r2 == 0xff || g1 == 0xff || g2 == 0xff
                || b1 == 0xff || b2 == 0xff {
            Err(())
        } else {
            Ok(Self {
            r: 16*r1 + r2,
            g: 16*g1 + g2,
            b: 16*b1 + b2,
            })
        }
    }
}

/// A palette of some pixel type. Guaranteed to have length
/// between 1 and 256 inclusive.
pub struct Palette<PixelType>(Vec<PixelType>);

impl<PixelType> Palette<PixelType> {
    /// Must have between `1` and `256` colours. Otherwise will return
    /// an error.
    pub fn new(colours: Vec<PixelType>) -> Result<Self, ()> {
        let lp = colours.len();
        if 1 <= lp && lp <= 256 {
            Ok(Self(colours))
        } else {
            Err(())
        }
    }
}

impl<PixelType: Ord> Palette<PixelType> {
    /// Sorts the palette by the `PixelType` partial ordering.
    pub fn sort(&mut self) {
        self.0.sort_unstable();
    }
}

impl<PixelType: Copy + Eq + Hash> Palette<PixelType> {
    /// Creates a new palette with duplicates removed. Note
    /// that the order of the resulting palette is not guaranteed
    /// and is not deterministic.
    pub fn remove_duplicates(&self) -> Self {
        Self(self.0.iter().map(|&c| c)
            .collect::<HashSet<_>>()
            .into_iter()
            .collect::<Vec<_>>())
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn oklab_ord_eq_impls() {
        let c1 = Oklabr { l: 1.0, a: 0.5, b: 0.75 };
        let c2 = Oklabr { l: 1.0, a: 0.5, b: 0.75 };
        let c3 = Oklabr { l: 1.0, a: 0.6, b: 0.75 };
        assert!(c1 == c2);
        assert!(c2 <= c3);
        assert!(c2 < c3);
    }

    #[test]
    fn rgb888_from_str() {
        let tests = [
            ("#ab6519", Ok(Rgb888 { r: 0xab, g: 0x65, b: 0x19})),
            ("#00B8FF", Ok(Rgb888 { r: 0x00, g: 0xb8, b: 0xff})),
            ("#ab65g9", Err(())),
        ];

        for (s, expected_result) in tests {
            let rgb = Rgb888::from_str(s);
            assert_eq!(expected_result, rgb);
        }
    }

    #[test]
    fn rgb888_to_str() {
        let tests = [
            ("#ab6519", Rgb888 { r: 0xab, g: 0x65, b: 0x19}),
            ("#00b8ff", Rgb888 { r: 0x00, g: 0xb8, b: 0xff}),
        ];

        for (expected_result, rgb) in tests {
            let s = rgb.to_str();
            assert_eq!(expected_result, s);
        }
    }

    #[test]
    fn test_oklabrify() {
        let color = Rgb888::from_str("#ab6519").unwrap().to_oklabr();
        let expected = Oklabr { l: 0.50523514, a: 0.05841799, b: 0.10804674 };
        assert_eq!(expected, color);
    }

    fn assert_eq_eps(a: f32, b: f32, eps: f32) {
        assert!((a - b).abs() <= eps, "{a} and {b} differed by more than {eps}");
    }

    fn assert_eq_eps_oklabr(a: &Oklabr, b: &Oklabr, eps: f32) {
        assert_eq_eps(a.l, b.l, eps);
        assert_eq_eps(a.a, b.a, eps);
        assert_eq_eps(a.b, b.b, eps);
    }

    fn assert_eq_eps_rgbf32(a: &RgbF32, b: &RgbF32, eps: f32) {
        assert_eq_eps(a.r, b.r, eps);
        assert_eq_eps(a.g, b.g, eps);
        assert_eq_eps(a.b, b.b, eps);
    }

    #[test]
    fn test_rgb888_oklab_roundtrip() {
        let rgb888 = Rgb888::from_str("#ab6519").unwrap();
        let rgbf32 = rgb888.to_rgbf32();
        let oklabr = rgbf32.to_oklabr();
        let rgbf32_back = oklabr.to_rgbf32();
        let rgb888_back = rgbf32_back.to_rgb888();

        let expected_rgbf32 = RgbF32 {
            r: 0xab as f32 / 255.0,
            g: 0x65 as f32 / 255.0,
            b: 0x19 as f32 / 255.0,
        };
        let expected_oklabr = Oklabr {
            l: 0.50523514,
            a: 0.05841787,
            b: 0.10804674,
        };

        let eps = 1.0e-6;
        assert_eq_eps_rgbf32(&expected_rgbf32, &rgbf32, eps);
        assert_eq_eps_rgbf32(&expected_rgbf32, &rgbf32_back, eps);
        assert_eq_eps_oklabr(&expected_oklabr, &oklabr, eps);
        assert_eq!(rgb888_back, rgb888);
    }
}
