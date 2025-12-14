use super::{Image, Oklabr, Palette};

/// A trait for dithering algorithms.
pub trait Dither {
    /// Dithers an image with the provided palette.
    ///
    /// Preconditions:
    /// - `palette` has length between `1` and `256` inclusive. (It may
    ///   contain duplicate colours.)
    /// - All colours in `palette` satisfy `Oklab::is_valid`.
    ///
    /// Postconditions:
    /// - The returned image has the same dimensions as `image`.
    /// - Each `u8` index entry is a valid index into `palette`.
    fn dither(&self, image: &Image<Oklabr>, palette: &Palette<Oklabr>) -> Image<u8>;

    /// Dithers multiple images with the provided palette. The
    /// default implementation is to call `Dither::dither` on each
    /// image separately, but this may be overridden for efficiency
    /// (for example, if constructing a data structure for faster
    /// palette comparisons).
    fn dither_mult(
        &self,
        images: &[Image<Oklabr>],
        palette: &Palette<Oklabr>,
    ) -> Vec<Image<u8>> {
        images.iter()
            .map(|im| self.dither(im, palette))
            .collect::<Vec<_>>()
    }
}

/// No dithering at all. Just maps each pixel to its closest
/// palette entry in `Oklabr` space.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct NoDither;

impl Dither for NoDither {
    fn dither(
        &self,
        image: &Image<Oklabr>,
        palette: &Palette<Oklabr>,
    ) -> Image<u8> {
        // TODO: Consider implementing a kd-tree for faster nearest
        // neighbour comparisons.
        // TODO: par-iter everything. Or implement on the GPU?
        image.map_pixels(|px| {
            palette.iter().enumerate()
                .min_by(|(_,p),(_,q)| p.sq_dist(&px)
                    .total_cmp(&q.sq_dist(&px)))
                .unwrap().0 as u8
        }).unwrap()
    }
}

/// Implements the Knoll ordered dithering algorithm with a Bayer
/// matrix. TODO: further description
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct KnollDither {
    num_samples: u16,
    strength: f32,
}

impl KnollDither {
    /// Returns new Knoll dither settings. `num_samples` must be
    /// in the set `{1, 2, 4, 8, 16, 32, 64, 128, 256}`. `strength` must
    /// lie in `[0.0, 1.0]`.
    pub fn new(num_samples: u16, strength: f32) -> Result<Self, ()> {
        let samples_valid = matches!(
            num_samples,
            1 | 2 | 4 | 8 | 16 | 32 | 64 | 128 | 256,
        );
        let strength_valid = 0.0 <= strength && strength <= 1.0;
        if samples_valid && strength_valid {
            Ok(Self { num_samples, strength })
        } else {
            Err(())
        }
    }

    fn knoll_dither_inner(
        &self,
        palette: &Palette<Oklabr>,
        pixel: &Oklabr,
    ) -> Vec<u8> {
        let n = self.num_samples as usize;
        let mut candidates = Vec::with_capacity(n);
        let mut error = Oklabr { l: 0.0, a: 0.0, b: 0.0 };
        for it in 0..n {
            let attempt = pixel + error * self.strength;
            let (argclosest, closest) = palette.iter().enumerate()
                .min_by(|(_,p),(_,q)| p.sq_dist(&attempt)
                    .total_cmp(&q.sq_dist(&attempt)))
                .unwrap();
            error += pixel - closest;
            candidates.push(argclosest as u8);
        }

        // The default order on Oklabr colours is increasing by luminance,
        // which is what we need here.
        candidates.sort_unstable_by_key(|&i| palette[i as usize]);
        candidates
    }
}

impl Dither for KnollDither {
    fn dither(
        &self,
        image: &Image<Oklabr>,
        palette: &Palette<Oklabr>,
    ) -> Image<u8> {
        let bayer_divisor = (256 / self.num_samples) as usize;
        // select relevant bayer matrix?
        image.map_pixels_idx(|(i,j),p| {
            let candidates = self.knoll_dither_inner(palette, p);
            let bayer_index = BAYER_16[j % 16][i % 16] as usize / bayer_divisor;
            candidates[bayer_index]
        }).unwrap()
    }
}

// 16x16 Bayer matrix. By construction of the Bayer matrix, this
// will be suitable for smaller values of N: by indexing and flooring
// it will degrade into an 8x8, 4x4 or 2x2 Bayer matrix.
const BAYER_16: [[u8; 16]; 16] = [
    [0, 128, 32, 160, 8, 136, 40, 168, 2, 130, 34, 162, 10, 138, 42, 170],
    [192, 64, 224, 96, 200, 72, 232, 104, 194, 66, 226, 98, 202, 74, 234, 106],
    [48, 176, 16, 144, 56, 184, 24, 152, 50, 178, 18, 146, 58, 186, 26, 154],
    [240, 112, 208, 80, 248, 120, 216, 88, 242, 114, 210, 82, 250, 122, 218, 90],
    [12, 140, 44, 172, 4, 132, 36, 164, 14, 142, 46, 174, 6, 134, 38, 166],
    [204, 76, 236, 108, 196, 68, 228, 100, 206, 78, 238, 110, 198, 70, 230, 102],
    [60, 188, 28, 156, 52, 180, 20, 148, 62, 190, 30, 158, 54, 182, 22, 150],
    [252, 124, 220, 92, 244, 116, 212, 84, 254, 126, 222, 94, 246, 118, 214, 86],
    [3, 131, 35, 163, 11, 139, 43, 171, 1, 129, 33, 161, 9, 137, 41, 169],
    [195, 67, 227, 99, 203, 75, 235, 107, 193, 65, 225, 97, 201, 73, 233, 105],
    [51, 179, 19, 147, 59, 187, 27, 155, 49, 177, 17, 145, 57, 185, 25, 153],
    [243, 115, 211, 83, 251, 123, 219, 91, 241, 113, 209, 81, 249, 121, 217, 89],
    [15, 143, 47, 175, 7, 135, 39, 167, 13, 141, 45, 173, 5, 133, 37, 165],
    [207, 79, 239, 111, 199, 71, 231, 103, 205, 77, 237, 109, 197, 69, 229, 101],
    [63, 191, 31, 159, 55, 183, 23, 151, 61, 189, 29, 157, 53, 181, 21, 149],
    [255, 127, 223, 95, 247, 119, 215, 87, 253, 125, 221, 93, 245, 117, 213, 85],
];

