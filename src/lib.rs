mod color;
mod dither;
mod img;
mod pximage;
//mod quantise;
mod vec2d;

pub use color::{Oklabr, Rgb888};
pub use pximage::{PxImage, write_indexed_png};
//pub use dither::{Dither};
//pub use quantise::{Quantise};
pub use vec2d::Vec2D;
pub use img::{Image, ImageError, ImageErrorInfo};

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







use std::io::Write;
use std::fs::File;
use std::io::BufWriter;

/// Writes an image in the PNG format using indexed colour (packed in
/// as few bits as possible). TODO: preconditions.
pub fn write_indexed_png_2<W: Write>(
    writer: W,
    indexed_image: Vec2D<u8>,
    palette: &[Rgb888],
) -> Result<(), ()> {
    todo!()
}

/// Convenience method TODO
pub fn write_indexed_png_to_path(
    path: &str,
    indexed_image: Vec2D<u8>,
    palette: &[Rgb888],
) -> Result<(), ()> {
    let file = File::create(path).unwrap();
    let mut w = BufWriter::new(file);
    write_indexed_png_2(&mut w, indexed_image, palette)
}
