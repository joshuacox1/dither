use std::ops::{Add, AddAssign, Sub, SubAssign, DivAssign, MulAssign};
use std::cmp::Ordering;
use std::fmt;

/// A colour in the Oklab colour space.
#[derive(Debug, Copy, Clone)]
pub struct Oklab { pub l: f32, pub a: f32, pub b: f32 }

// TODO: for good style, make add etc. traits take references.

impl Oklab {
    /// Convert to gamma-encoded RGB in the sRGB colour space.
    pub fn to_srgb(&self) -> Srgb {
        let l_ = self.l + 0.3963377774f32*self.a + 0.2158037573f32*self.b;
        let m_ = self.l - 0.1055613458f32*self.a - 0.0638541728f32*self.b;
        let s_ = self.l - 0.0894841775f32*self.a - 1.2914855480f32*self.b;

        let l = l_*l_*l_;
        let m = m_*m_*m_;
        let s = s_*s_*s_;

        let r_ = 4.0767416621f32*l - 3.3077115913f32*m + 0.2309699292f32*s;
        let g_ = -1.2684380046f32*l + 2.6097574011f32*m - 0.3413193965f32*s;
        let b_ = -0.0041960863f32*l - 0.7034186147f32*m + 1.7076147010f32*s;

        Srgb {
            r: Srgb::delinearise(r_),
            g: Srgb::delinearise(g_),
            b: Srgb::delinearise(b_),
        }
    }

    /// Squared distance between colours.
    pub fn sq_dist(&self, other: &Self) -> f32 {
        let l_d = other.l - self.l;
        let a_d = other.a - self.a;
        let b_d = other.b - self.b;
        l_d*l_d + a_d*a_d + b_d*b_d
    }

    pub const WHITE: Oklab = Oklab { l: 1.0, a: 0.0, b: 0.0 };
    pub const BLACK: Oklab = Oklab { l: 0.0, a: 0.0, b: 0.0 };

    /// Panics if `rgb` is not of length 3.
    pub fn from_srgb_triple(rgb: &[f32]) -> Self {
        assert!(rgb.len() == 3);
        Srgb::to_oklab(&Srgb { r: rgb[0], g: rgb[1], b: rgb[2] })
    }

    /// From a string e.g. #ab6519 or #AB6519.
    pub fn from_srgb_888_str(hash_rrggbb: &str) -> Result<Self, ()> {
        Srgb::from_srgb_888_str(hash_rrggbb)
            .map(|c| Srgb::to_oklab(&c))
    }
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

impl fmt::Display for Oklab {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:.2}%,{:.4},{:.4})", self.l*100.0, self.a, self.b)
    }
}

/// Gamma-encoded sRGB.
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

    fn delinearise(x: f32) -> f32 {
        if x >= 0.0031308 {
            1.055 * x.powf(1.0/2.4) - 0.055
        } else {
            12.92 * x
        }
    }

    // https://bottosson.github.io/posts/oklab/#converting-from-linear-srgb-to-oklab
    fn to_oklab(&self) -> Oklab {
        // linearise
        let r = Self::linearise(self.r);
        let g = Self::linearise(self.g);
        let b = Self::linearise(self.b);

        // oklabify.
        let l = 0.4122214708f32*r + 0.5363325363f32*g + 0.0514459929f32*b;
        let m = 0.2119034982f32*r + 0.6806995451f32*g + 0.1073969566f32*b;
        let s = 0.0883024619f32*r + 0.2817188376f32*g + 0.6299787005f32*b;

        let l_ = l.cbrt();
        let m_ = m.cbrt();
        let s_ = s.cbrt();

        return Oklab {
            l: 0.2104542553f32*l_ + 0.7936177850f32*m_ - 0.0040720468f32*s_,
            a: 1.9779984951f32*l_ - 2.4285922050f32*m_ + 0.4505937099f32*s_,
            b: 0.0259040371f32*l_ + 0.7827717662f32*m_ - 0.8086757660f32*s_,
        };
    }

    pub fn to_srgb888(&self) -> [u8; 3] {
        fn discretise_8bit(x: f32) -> u8 {
            (x * 255.0).clamp(0.0, 255.0).round() as u8
        }

        [
            discretise_8bit(self.r),
            discretise_8bit(self.g),
            discretise_8bit(self.b),
        ]
    }

    pub fn to_srgb555(&self) -> u16 {
        fn discretise_5bit(x: f32) -> u16 {
            (x * 31.0).clamp(0.0, 31.0).round() as u16
        }
        let r_ = discretise_5bit(self.r);
        let g_ = discretise_5bit(self.g);
        let b_ = discretise_5bit(self.b);
        let u = (r_ << 10) | (g_ << 5) | b_;
        u
    }

    /// Parses a string of the form #xxyyzz where each pair is a valid
    /// hex number.
    fn from_srgb_888_str(input: &str) -> Result<Self, ()> {
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
            return Err(());
        }

        let r = (16*r1 + r2) as f32 / 255.0;
        let g = (16*g1 + g2) as f32 / 255.0;
        let b = (16*b1 + b2) as f32 / 255.0;

        Ok(Srgb { r, g, b })
    }
}

// TODO: delete?
// #[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
// pub struct Srgb555(pub u16);

// impl Srgb555 {
//     pub fn to_srgb(self) -> Srgb {
//         let r_ = (self.0 & 0b0_11111_00000_00000) >> 10;
//         let g_ = (self.0 & 0b0_00000_11111_00000) >> 5;
//         let b_ = self.0 & 0b0_00000_00000_11111;
//         let r = r_ as f32 / 31.0;
//         let g = g_ as f32 / 31.0;
//         let b = b_ as f32 / 31.0;
//         Srgb { r, g, b }
//     }
// }


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

    // TODO: make sure we're doing gamma the right way round.
    #[test]
    fn oklab_makesure() {
        let c1 = "#ab6519"; let c2 = "#AB6519";
        let ok1 = Oklab::from_srgb_888_str(c1).unwrap();
        let ok2 = Oklab::from_srgb_888_str(c2).unwrap();
        let expected = Oklab { l: 0.57339257, a: 0.05841787, b: 0.10804674 };

        assert_eq!(expected, ok1);
        assert_eq!(expected, ok2);
    }
}
