use std::borrow::Cow;

use super::{Oklabr, Image, Vec2D};

/// A trait for extracting a representative limited palette from
/// a list of images.
pub trait Quantise {
    /// Given a slice of images, return a limited palette for use in
    /// image quantisation. The limited palette must have between
    /// `1` and `256` colours (though no need to guarantee lack
    /// of duplicates).
    fn quantise(
        &self,
        images: &[Image<Oklabr>]
    ) -> Cow<'_, [Oklabr]>;
}

/// The trivial quantisation method of providing an explicit palette.
pub struct ExplicitPalette {
    palette: Vec<Oklabr>,
}

impl ExplicitPalette {
    /// `colours` must have between `1` and `256` members inclusive.
    pub fn new(palette: Vec<Oklabr>) -> Result<Self, ()> {
        let l_p = palette.len();
        if 1 <= l_p && l_p <= 256 && palette.iter().all(|p| p.is_well_formed()) {
            Ok(Self { palette })
        } else {
            Err(())
        }
    }
}

impl Quantise for ExplicitPalette {
    fn quantise(&self, images: &[Image<Oklabr>]) -> Cow<'_, [Oklabr]> {
        Cow::from(&self.palette)
    }
}

/// Quantisation method using the k-means clustering algorithm.
pub struct KmeansClustering {
    k: u16,
    random_seed: u64,
}

impl KmeansClustering {
    /// `k` must be between `1` and `256`.
    /// All `u64` values are valid and reasonable random seeds.
    pub fn new(k: u16, random_seed: u64) -> Result<Self, ()> {
        if 1 <= k && k <= 256 {
            Ok(Self { k, random_seed })
        } else {
            Err(())
        }
    }
}

impl Quantise for KmeansClustering {
    fn quantise(&self, images: &[Image<Oklabr>]) -> Cow<'_, [Oklabr]> {
        todo!()
    }
}