use super::{Oklabr, Image, Vec2D};

/// A trait for image colour quantisation.
pub trait Quantise {
    /// Given a list of images (whose width and height will both be positive),
    /// returns a limited palette of size between `1` and `n+1` inclusive
    /// that represents the colours of the image well. TODO rename `n`
    /// TODO: return a Cow?
    fn quantise(&self, image: Image<Oklabr>) -> Vec<Oklabr>;
}


/// Quantisation method using the k-means clustering algorithm.
pub struct KmeansClustering {
    k: u16,
    random_seed: u64,
}

impl KmeansClustering {
    /// All random seed values are valid and reasonable.
    pub fn new(k: u16, random_seed: u64) -> Result<Self, ()> {
        if 1 <= k && k <= 256 {
            Ok(Self { k, random_seed })
        } else {
            Err(())
        }
    }
}

impl Quantise for KmeansClustering {
    fn quantise(&self, images: &[Vec2D<Oklabr>]) -> Vec<Oklabr> {
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
    fn quantise(&self, image: Vec2D<Oklabr>) -> Vec<Oklabr> {
        self.palette.clone()
    }
}
