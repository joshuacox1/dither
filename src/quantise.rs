use std::borrow::Cow;

use super::{Oklabr, Image, Vec2D, Palette};

use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;

/// A trait for extracting a representative limited palette from
/// a list of images.
pub trait Quantise {
    /// Given a slice of images, return a limited palette for use in
    /// image quantisation.
    fn quantise_mult(
        &self,
        images: &[Image<Oklabr>]
    ) -> Cow<'_, Palette<Oklabr>>;

    /// Quantises a single image. Calls `quantise_mult` by default internally.
    fn quantise(&self, image: &Image<Oklabr>) -> Cow<'_, Palette<Oklabr>> {
        self.quantise_mult(std::slice::from_ref(image))
    }
}

impl Quantise for Palette<Oklabr> {
    fn quantise_mult(&self, images: &[Image<Oklabr>]) -> Cow<'_, Palette<Oklabr>> {
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

    // Naive k-means on Oklabr colours.
    // Precondition: points is non-empty.
    fn naive_k_means(
        &self,
        points: &[Oklabr],
        rng: &mut impl Rng,
    ) -> Vec<Oklabr> {
        let k = self.k as usize;
        let mut centroids = self.kmeans_initial(points, rng);
        if centroids.len() < k {
            return centroids;
        }

        let mut centroid_idxes = vec![usize::MAX; points.len()];
        // TODO: use doubles? are we worried about f32 accuracy here?
        let mut cluster_counts = vec![0; k];
        let mut converged;

        let mut counter = 0;
        loop {
            converged = true;

            // Assign each point to the cluster with the nearest centroid.
            // TODO: check we can't oscillate...
            for (i_p, p) in points.iter().enumerate() {
                let min_centroid_idx = centroids.iter().enumerate()
                    .min_by(|(_,m1),(_,m2)| p.sq_dist(m1).total_cmp(&p.sq_dist(m2)))
                    .unwrap().0;

                if centroid_idxes[i_p] != min_centroid_idx {
                    converged = false;
                }
                centroid_idxes[i_p] = min_centroid_idx;
                cluster_counts[min_centroid_idx] += 1;
            }

            if converged {
                break;
            }

            // Recalculate the centroids of each cluster.
            for c in centroids.iter_mut() {
                // todo: just modify directly
                *c = Oklabr { l: 0.0, a: 0.0, b: 0.0 };
            }
            for (i_p, p) in points.iter().enumerate() {
                let centroid_idx = centroid_idxes[i_p];
                centroids[centroid_idx] += *p;
            }
            for (i_c, c) in centroids.iter_mut().enumerate() {
                let cluster_size = cluster_counts[i_c];
                // TODO: delete this assertion. It may not always be true.
                // make it robust
                assert!(cluster_size > 0);
                *c /= cluster_size as f32;
            }

            for c in cluster_counts.iter_mut() { *c = 0; }
            counter += 1;
        }

        centroids
    }

    /// Uses the kmeans++ method to generate a well-behaved initial
    /// set of centroids for the kmeans clustering algorithm.
    /// Returns a Vec of at MOST length k. This function should only
    /// return a Vec of length less than k if there are fewer than k
    /// unique points. In this case we can skip the whole kmeans
    /// algorithm as these points are a perfect cluster.
    /// Will always return at least one element.
    fn kmeans_initial(
        &self,
        points: &[Oklabr],
        rng: &mut impl Rng,
    ) -> Vec<Oklabr> {
        let k = self.k as usize;
        // Initialise list of centroids with a random point.
        let mut centroids = Vec::with_capacity(k);
        let idx0 = rng.random_range(..points.len());
        centroids.push(points[idx0]);

        let mut sq_dists = Vec::with_capacity(points.len());
        for _ in 1..k {
            // For each point, compute D^2, the squared distance to
            // the nearest centroid to it in the current list.
            let mut total = 0.0f32;
            for point in points {
                let min_sq_dist = centroids.iter()
                    .map(|c| point.sq_dist(c))
                    .min_by(f32::total_cmp)
                    .unwrap();
                sq_dists.push(min_sq_dist);
                total += min_sq_dist;
            }

            // If `total` is zero, then there are fewer than `k`
            // unique points. In this case return early. This is
            // then a perfect palette and we can skip kmeans.
            if total == 0.0 {
                return centroids;
            }

            // Pick a random point with probability proportional
            // to D^2.
            let threshold = rng.random::<f32>() * total;
            let mut cumulative = 0.0f32;
            let mut point_picked = false;
            for (i, point) in points.iter().enumerate() {
                cumulative += sq_dists[i];
                if cumulative >= threshold {
                    centroids.push(*point);
                    point_picked = true;
                    break;
                }
            }
            if !point_picked {
                // todo: add a warning log around this.
                // this should only occur if there is floating
                // point wonkiness.
                centroids.push(*points.last().unwrap());
            }

            sq_dists.clear();
        }

        centroids
    }
}

impl Quantise for KmeansClustering {
    fn quantise_mult(&self, images: &[Image<Oklabr>]) -> Cow<'_, Palette<Oklabr>> {
        let capacity = images.iter().map(|img| {
            let (w, h) = img.dims();
            w*h
        }).sum();
        let mut all_pixels = Vec::with_capacity(capacity);
        for image in images {
            for frame in image.frames() {
                for pixel in frame.data_1d() {
                    all_pixels.push(*pixel);
                }
            }
        }

        let mut rng = Xoshiro256StarStar::seed_from_u64(self.random_seed);
        let centroids = self.naive_k_means(&all_pixels, &mut rng);
        let palette = Palette::new(centroids).unwrap();
        Cow::Owned(palette)
    }
}