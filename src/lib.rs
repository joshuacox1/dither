mod color;
mod pximage;

pub use color::Oklab;
pub use pximage::PxImage;

#[derive(Debug, Clone)]
pub enum PaletteOptions {
    /// An explicit list of (unique) colours.
    Explicit(Vec<Oklab>),
    /// Automatically selects a palette based on the input image.
    /// (Uses k-means clustering with kmeans++ initialisation.)
    Auto {
        /// The number of colors in the image. Zero is an invalid value.
        num_colors: u8,
        /// A seed for the random number generator used in palette
        /// generation.
        random_seed: u64,
    },
}

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

#[derive(Debug, Clone)]
pub struct Problem {
    pub input_path: String,
    pub palette: PaletteOptions,
    pub dither: DitherOptions,
    pub output_path: String,
}
