use super::{Image, Vec2D, Oklabr};

/// A trait for dithering algorithms.
pub trait Dither {
    /// Dithers an image with the provided palette.
    ///
    /// Preconditions:
    ///
    /// - `image` has both width and height at least `1`.
    /// - `palette` has length between `1` and `256` inclusive and is
    ///   sorted according to the `Oklabr` `Ord` implementation. (It may
    ///   contain duplicate colours.)
    /// - All `Oklabr` colours in `image` and `palette` are _valid_
    ///   (i.e. `l_r`, `a`, and `b` are all finite).
    ///
    /// Postconditions:
    ///
    /// - The returned image has the same dimensions as `image`.
    /// - Each `u8` index entry is a valid index into `palette`.
    fn dither_mult(
        &self,
        images: &[Image<Oklabr>],
        palette: &[Oklabr],
    ) -> Vec<Image<u8>>;

    /// Dither a single image.
    /// TODO: what?
    /// Make it take an impl?
    fn dither(&self, image: &Image<Oklabr>, palette: &[Oklabr]) -> Image<u8> {
        self.dither_mult(std::slice::from_ref(image), palette)
            .into_iter().next().unwrap()
    }
}



// todo: impl Eq or something
/// The Knoll ordered dithering algorithm.
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
    fn dither_mult(
        &self,
        images: &[Image<Oklabr>],
        palette: &[Oklabr],
    ) -> Vec<Image<u8>> {
        todo!()
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct NearestNeighbour;

impl Dither for NearestNeighbour {
    fn dither_mult(
        &self,
        images: &[Image<Oklabr>],
        palette: &[Oklabr],
    ) -> Vec<Image<u8>> {
        todo!()
        // TODO: Consider implementing a kd-tree for faster nearest
        // neighbour comparisons.
        // TODO: par-iter everything. Or implement on the GPU?
        // image.map_frames(|img| {
        //     todo!()
        //     // map the below over every pixel
        //     // let (argclosest, closest) = palette.iter().enumerate()
        //     //     .min_by(|(_,p),(_,q)| p.sq_dist(&attempt)
        //     //         .total_cmp(&q.sq_dist(&attempt)))
        //     //     .unwrap();

        // })
    }
}


