use std::borrow::Cow;

use super::{Oklabr, Image, Vec2D, Palette};

/// A trait for extracting a representative limited palette from
/// a list of images.
pub trait Quantise {
    /// Given a slice of images, return a limited palette for use in
    /// image quantisation.
    fn quantise(
        &self,
        images: &[Image<Oklabr>]
    ) -> Cow<'_, Palette<Oklabr>>;
}

impl Quantise for Palette<Oklabr> {
    fn quantise(&self, images: &[Image<Oklabr>]) -> Cow<'_, Palette<Oklabr>> {
        Cow::Borrowed(&self)
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
    fn quantise(&self, images: &[Image<Oklabr>]) -> Cow<'_, Palette<Oklabr>> {
        todo!()
    }
}