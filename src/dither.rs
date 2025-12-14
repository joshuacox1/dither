use super::{Image, Oklabr, Palette, Vec2D};

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
        image.map_pixels(|px| palette.nearest_neighbour(px) as u8)
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
            let argclosest = palette.nearest_neighbour(&attempt);
            error += pixel - palette[argclosest];
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
        })
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

/// Implements error-diffusion dithering.
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum ErrorDiffusionDither {
    /// blah.
    FloydSteinberg,
    /// blah.
    JarvisJudiceNinke,
}

impl Dither for ErrorDiffusionDither {
    fn dither(
        &self,
        image: &Image<Oklabr>,
        palette: &Palette<Oklabr>,
    ) -> Image<u8> {
        let diffusion_map = match self {
            ErrorDiffusionDither::FloydSteinberg => FLOYD_STEINBERG.as_slice(),
            ErrorDiffusionDither::JarvisJudiceNinke => JARVIS_JUDICE_NINKE.as_slice(),
        };

        image.map_frames(|frame| {
            let mut buf = frame.clone();
            let (w, h) = buf.dims();
            let mut result_1d = Vec::with_capacity(w*h);
            for j in 0..h {
                for i in 0..w {
                    let old_px = buf[(i,j)];
                    let new_px = palette.nearest_neighbour(&old_px);
                    result_1d.push(new_px as u8);
                    let error = old_px - palette[new_px]; // reverse error

                    for ((i_offset, j_offset), coeff) in diffusion_map.iter() {
                        let new_i = i.checked_add_signed(*i_offset);
                        let new_j = j.checked_add_signed(*j_offset);
                        match (new_i, new_j) {
                            (Some(new_i), Some(new_j)) => {
                                if let Some(px_) = buf.get_mut((new_i, new_j)) {
                                    *px_ += error * coeff;
                                }
                            },
                            _ => (),
                        };
                    }
                }
            }

            Vec2D::from_1d(result_1d, (w, h)).unwrap()
        }).unwrap()
    }
}


const FLOYD_STEINBERG: [((isize, isize), f32); 4] = [
    ((1,0), 7.0/16.0),
    ((-1,1), 3.0/16.0),
    ((0,1), 5.0/16.0),
    ((1,1), 1.0/16.0),
];

const JARVIS_JUDICE_NINKE: [((isize, isize), f32); 12] = [
    ((1,0), 7.0/48.0),
    ((2,0), 5.0/48.0),
    ((-2,1), 3.0/48.0),
    ((-1,1), 5.0/48.0),
    ((0,1), 7.0/48.0),
    ((1,1), 5.0/48.0),
    ((2,1), 3.0/48.0),
    ((-2,2), 1.0/48.0),
    ((-1,2), 3.0/48.0),
    ((0,2), 5.0/48.0),
    ((1,2), 3.0/48.0),
    ((2,2), 1.0/48.0),
];
