mod color;
mod pximage;

pub use color::{Oklabr, Rgb888};
pub use pximage::{PxImage, write_indexed_png};

// TODO: delete
// #[derive(Debug, Clone)]
// pub enum PaletteOptions {
//     /// An explicit list of (unique) colours.
//     Explicit(Vec<Oklabr>),
//     /// Automatically selects a palette based on the input image.
//     /// (Uses k-means clustering with kmeans++ initialisation.)
//     Auto {
//         /// The number of colors in the image. Zero is an invalid value.
//         num_colors: u8,
//         /// A seed for the random number generator used in palette
//         /// generation.
//         random_seed: u64,
//     },
// }

#[derive(Debug, Clone)]
pub struct DitherOptions {
    /// Number of samples taken of the dither matrix.
    /// I believe this number should lie between B^2 / 4 and B^2,
    /// where B = matrix size. Consdier
    pub num_samples: u8,
    /// Dither matrix side length. Valid values are 2, 4, 8 and 16.
    pub matrix_size: u8,
    /// Threshold. TODO: find a
    /// geometric intuition for any threshold. And offer a limited
    /// number of thresholds (e.g. 16).
    pub strength: f32,
}

// #[derive(Debug, Clone)]
// pub struct Problem {
//     pub input_path: String,
//     pub palette: PaletteOptions,
//     pub dither: DitherOptions,
//     pub output_path: String,
// }


/// A two-dimensional array of fixed size determined at runtime.
pub struct Array2D<T> {
    data: Vec<T>,
    dims: (usize, usize),
}


/// A trait for colour quantisation algorithms.
pub trait Quantise {
    /// Given an image (whose width and height will both be positive),
    /// returns a limited palette of size between `1` and `n+1` inclusive
    /// that represents the colours of the image well. TODO rename `n`
    fn quantise(&self, image: Array2D<Oklabr>, n: u8) -> Vec<Oklabr>;
}

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
    fn dither(
        &self,
        image: Array2D<Oklabr>,
        palette: &[Oklabr],
    ) -> Array2D<u8>;
}

/// Quantisation method using the k-means clustering algorithm.
pub struct KmeansClustering {
    random_seed: u64,
}

impl KmeansClustering {
    /// All random seed values are valid and reasonable.
    pub fn new(random_seed: u64) -> Self {
        Ok(Self { random_seed })
    }
}

impl Quantise for KmeansClustering {
    fn quantise(&self, image: Array2D<Oklabr>, n: u8) -> Vec<Oklabr> {
        todo!()
    }
}

/// The trivial quantisation method of providing an explicit palette.
pub struct ExplicitPalette {
    palette: Vec<Oklabr>,
}

impl ExplicitPalette {
    /// `colours` must have between `1` and `256` members inclusive.
    pub fn new(palette: Vec<Oklabr>) -> Result<Self, ()> {
        let l_p = palette.len();
        if 1 <= l_p && l_p <= 256 {
            Ok(Self { palette })
        } else {
            Err(())
        }
    }
}

impl Quantise for ExplicitPalette {
    fn quantise(&self, image: Array2D<Oklabr>, n: u8) -> Vec<Oklabr> {
        self.palette.clone()
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
    fn dither(
        &self,
        image: Array2D<Oklabr>,
        palette: &[Oklabr],
    ) -> Array2D<u8> {
        todo!()
    }
}








use std::io::Write;
use std::fs::File;
use std::io::BufWriter;

/// Writes an image in the PNG format using indexed colour (packed in
/// as few bits as possible). TODO: preconditions.
pub fn write_indexed_png_2<W: Write>(
    writer: W,
    indexed_image: Array2D<u8>,
    palette: &[Rgb888],
) -> Result<(), ()> {
    todo!()
}

/// Convenience method TODO
pub fn write_indexed_png_to_path(
    path: &str,
    indexed_image: Array2D<u8>,
    palette: &[Rgb888],
) -> Result<(), ()> {
    let file = File::create(path).unwrap();
    let mut w = BufWriter::new(file);
    write_indexed_png_2(&mut w, indexed_image, palette)
}
