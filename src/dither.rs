use super::{Image, Oklabr};

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
    fn dither(&self, image: &Image<Oklabr>, palette: &[Oklabr]) -> Image<u8>;

    /// Dithers multiple images with the provided palette. The
    /// default implementation is to call `Dither::dither` on each
    /// image separately, but this may be overridden for efficiency
    /// (for example, if constructing a data structure for faster
    /// palette comparisons).
    fn dither_mult(
        &self,
        images: &[Image<Oklabr>],
        palette: &[Oklabr],
    ) -> Vec<Image<u8>> {
        images.iter()
            .map(|im| self.dither(im, palette))
            .collect::<Vec<_>>()
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
}

impl Dither for KnollDither {
    fn dither(
        &self,
        image: &Image<Oklabr>,
        palette: &[Oklabr],
    ) -> Image<u8> {
        todo!()
    }
}

/// No dithering at all. Just maps each pixel to its closest
/// palette entry in Oklabr space.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct NoDither;

impl Dither for NoDither {
    fn dither(
        &self,
        image: &Image<Oklabr>,
        palette: &[Oklabr],
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


