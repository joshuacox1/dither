//! Wrappers to convert to and from PNG.
//! It would be nice to support APNG at some point. But for now, just PNG.

pub use super::types::{Oklab, Srgb};
use image::{ImageReader, Rgb, ImageBuffer};

use std::ops::{Index, IndexMut};

use rand::{Rng, SeedableRng};
use rand_xoshiro::{Xoshiro256StarStar};

// /// A 2D array of stuffs.
pub struct OkImage {
    data: Vec<Oklab>,
    dims: (usize, usize),
}


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

pub struct DitherOptions {
    /// Dither matrix side length. Valid values are 2, 4, 8 and 16.
    matrix_size: u8,
    /// Number of samples taken of the dither matrix.
    /// I believe this number should lie between B^2 / 4 and B^2,
    /// where B = matrix size. Consdier
    num_samples: u8,
}

pub struct Problem {
    input_path: String,
    palette: PaletteOptions,
    dither: DitherOptions,
    output_path: String,
}


pub fn image_from_path() {

}








impl OkImage {
    pub fn from_path(s: &str) -> Self {
        // todo: remove these unwraps lol.
        let img = ImageReader::open(s)
            .unwrap().decode().unwrap()
            .into_rgb32f();
        let w = img.width() as usize;
        let h = img.height() as usize;
        let data = img.pixels()
            .map(|rgb| {
                let srgb = Srgb { r: rgb[0], g: rgb[1], b: rgb[2] };
                srgb.to_linear().to_oklab()
            })
            .collect::<Vec<_>>();
        Self::new(data, (w, h))
    }

    pub fn save_as_srgb_png(&self, s: &str) {
        // rayon?!
        let img = ImageBuffer::from_fn(
            self.width() as u32,
            self.height() as u32,
            |i,j| {
                let (i,j) = (i as usize, j as usize);
                let px = self[(i,j)];
                let srgb = px.to_linear_srgb().to_gamma();
                Rgb::<u8>([
                    (srgb.r * 255.0).clamp(0.0, 255.0) as u8,
                    (srgb.g * 255.0).clamp(0.0, 255.0) as u8,
                    (srgb.b * 255.0).clamp(0.0, 255.0) as u8,
                ])
            });
        img.save(s).unwrap();
    }

    pub fn new(data: Vec<Oklab>, dims: (usize, usize)) -> Self {
        assert!(data.len() == dims.0*dims.1);
        assert!(data.len() > 0);
        Self { data, dims }
    }

    pub fn width(&self) -> usize { self.dims.0 }
    pub fn height(&self) -> usize { self.dims.1 }
    pub fn dims(&self) -> (usize, usize) { self.dims }

    pub fn round(&mut self, palette: &[Oklab]) {
        for d in self.data.iter_mut() {
            let rounded = palette.iter()
                .min_by(|p,q| p.sq_dist(d).total_cmp(&q.sq_dist(d)))
                .unwrap();
            *d = *rounded;
        }
    }

    pub fn knoll_dither(&mut self,
        n: usize,
        bayer_size: usize,
        palette: &[Oklab],
    ) {
        let (w, h) = self.dims;
        let mut candidates = vec![Oklab::BLACK; n];
        let threshold = 0.3f32;
        let bayer_interped = match bayer_size {
            4 => interp_bayer(BAYER_4, n),
            8 => interp_bayer(BAYER_8, n),
            16 => interp_bayer(BAYER_16, n),
            _ => panic!("Invalid Bayer matrix size {bayer_size}."),
        };
        for j in 0..h {
            for i in 0..w {
                let pixel = self[(i,j)];
                let mut goal = pixel;
                for it in 0..n {
                    // closest colour to goal. todo: put this out
                    let closest = palette.iter()
                        .min_by(|p,q| p.sq_dist(&goal).total_cmp(&q.sq_dist(&goal)))
                        .unwrap();
                    let mut diff = pixel - *closest; diff *= threshold;
                    if i == 0 && j == 0 {
                        println!("{pixel:?}, {goal:?}");
                    }

                    goal += diff;
                    candidates[it] = *closest;
                }
                // sort by luminance
                candidates.sort_unstable_by(|c1,c2|
                    c1.l.total_cmp(&c2.l));


                let bayer_idx = (j % bayer_size)*bayer_size + (i % bayer_size);
                let bayer_interp = bayer_interped[bayer_idx];
                let chosen = candidates[bayer_interp];

                if i == 0 && j == 0 {
                    println!("{candidates:?}, {bayer_idx}, {bayer_interp}");
                }

                //println!("{i},{j},{candidates:?},{bayer_idx},{bayer_interp}");
                self[(i,j)] = chosen;
            }
        }
    }

    // // Overwrites candidates.
    // fn knoll_dither_inner(
    //     candidates: &mut [Oklab],
    //     palette: &[Oklab],
    //     pixel: &Oklab,
    //     n: usize,
    //     threshold: f32)
    // {
    //     let mut goal = pixel;
    //     for it in 0..n {
    //         // closest colour to goal. todo: put this out
    //         let closest = palette.iter()
    //             .min_by(|p,q| p.sq_dist(&goal).total_cmp(&q.sq_dist(&goal)))
    //             .unwrap();
    //         let mut diff = pixel - *closest; diff *= threshold;
    //         goal += diff;
    //         candidates[it] = *closest;
    //     }

    //     // The default order on Oklab colours is increasing by luminance,
    //     // which is what we need here.
    //     candidates.sort_unstable();
    // }

    /// Obtain a representative palette from the image of at most k colours.
    /// k must be positive. The palette will generally consist of
    /// exactly k colours unless the original image has fewer colours.
    /// Performs kmeans clustering with kmeans++ initialisation.
    /// (todo: should I transpose from Oklab back into eg srgb555,
    /// then back again?)
    /// The palette will have at least one colour in it.
    pub fn palette(&self, k: usize, random_seed: u64) -> Vec<Oklab> {
        let mut rng = Xoshiro256StarStar::seed_from_u64(random_seed);

        let mut centroids = Self::naive_k_means(&self.data, k, &mut rng);
        // TODO: consider ensuring here that the palette colours
        // are valid eg SRGB555 or RGB888 colours by mapping there,
        // rounding and mapping back to Oklab. This will mean our
        // palette is correct for the dithering computations.
        // Probably doesn't make a difference but we should probably
        // round somewhere.
        centroids.sort_unstable();
        // for c in centroids {

        // }
        centroids
    }

    // Naive k-means on Oklab colours.
    // Precondition: points is non-empty and k > 0.
    fn naive_k_means(
        points: &[Oklab],
        k: usize,
        rng: &mut impl Rng,
    ) -> Vec<Oklab> {
        let mut centroids = Self::kmeans_initial(points, k, rng);
        if centroids.len() < k {
            return centroids;
        }

        let mut centroid_idxes = vec![usize::MAX; points.len()];
        // TODO: use doubles? are we worried about f32 accuracy here?
        let mut cluster_counts = vec![0; k];
        let mut converged = false;

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
            //println!("kmeans step. part: {centroid_idxes:?}\ncounts:{cluster_counts:?}");

            if converged {
                break;
            }

            // Recalculate the centroids of each cluster.
            for c in centroids.iter_mut() {
                *c = Oklab { l: 0.0, a: 0.0, b: 0.0 };
            }
            for (i_p, p) in points.iter().enumerate() {
                let centroid_idx = centroid_idxes[i_p];
                centroids[centroid_idx] += *p;
            }
            for (i_c, c) in centroids.iter_mut().enumerate() {
                let cluster_size = cluster_counts[i_c];
                // TODO: delete this assertion.
                assert!(cluster_size > 0);
                *c /= cluster_size as f32;
            }

            for c in cluster_counts.iter_mut() { *c = 0; }
        }

        centroids
    }

    /// Uses the kmeans++ method to generate a well-behaved initial
    /// set of centroids for the kmeans clustering algorithm.
    /// We must have 0 < k and 0 < points.len().
    /// Returns a Vec of at MOST length k. This function will only
    /// return a Vec of length less than k if there are fewer than k
    /// unique points. In this case we can skip the whole kmeans
    /// algorithm as these points are a perfect cluster.
    /// Will always return at least one element.
    fn kmeans_initial(
        points: &[Oklab],
        k: usize,
        rng: &mut impl Rng,
    ) -> Vec<Oklab> {
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

            // TODO: Floating point error could maybe make this strange?
            // Since cumulative and threshold do the exact same calc
            // (adding all the numbers in the same order) you'd think
            // this is guaranteed. In practice you could imagine floating
            // point wonkiness. So really should accept this may not
            // be true and select the last point if so.
            assert!(point_picked);
            // if !point_picked {
            //     centroids.push(*points.last().unwrap());
            // }
            sq_dists.clear();
        }

        centroids
    }
}

impl Index<(usize, usize)> for OkImage {
    type Output = Oklab;

    fn index(&self, ij: (usize, usize)) -> &Self::Output {
        let (i, j) = ij; let (width, height) = self.dims;
        if i < width && j < height {
            &self.data[width*j + i]
        } else {
            panic!("Failed to index img with dims {:?} with indices {ij:?}",
                self.dims);
        }
    }
}

impl IndexMut<(usize, usize)> for OkImage {
    fn index_mut(&mut self, ij: (usize, usize)) -> &mut Self::Output {
        let (i, j) = ij; let (width, height) = self.dims;
        if i < width && j < height {
            &mut self.data[width*j + i]
        } else {
            panic!("Failed to index img with dims {:?} with indices {ij:?}",
                self.dims);
        }
    }
}

const BAYER_4: [[u8; 4]; 4] = [
    [0, 8, 2, 10],
    [12, 4, 14, 6],
    [3, 11, 1, 9],
    [15, 7, 13, 5],
];

const BAYER_8: [[u8; 8]; 8] = [
    [0, 32, 8, 40, 2, 34, 10, 42],
    [48, 16, 56, 24, 50, 18, 58, 26],
    [12, 44, 4, 36, 14, 46, 6, 38],
    [60, 28, 52, 20, 62, 30, 54, 22],
    [3, 35, 11, 43, 1, 33, 9, 41],
    [51, 19, 59, 27, 49, 17, 57, 25],
    [15, 47, 7, 39, 13, 45, 5, 37],
    [63, 31, 55, 23, 61, 29, 53, 21],
];

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

// returns a 2D array encoded in a 1D one.
fn interp_bayer<const B: usize>(bayer_mat: [[u8; B]; B], knoll_n: usize)
-> Vec<usize> {
    let mut result = Vec::with_capacity(B*B);
    for j in 0..B {
        for i in 0..B {
            let b = bayer_mat[j][i];
            let b_interp = ((b as f32) * (knoll_n as f32 / (B*B) as f32))
                .floor() as usize;
            result.push(b_interp);
        }
    }
    result
}




#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_kmeans_rect_full() {
        let img = OkImage {
            data: vec![
                Oklab { l: 0.5, a: 0.0, b: 0.0 },
                Oklab { l: 0.5, a: 0.0, b: 0.5 },
                Oklab { l: 0.5, a: 1.0, b: 0.0 },
                Oklab { l: 0.5, a: 1.0, b: 0.5 },
            ],
            dims: (2,2),
        };
        let k = 2usize;
        let seed = 0xdeadbeef;
        let q = img.palette(k, seed);
        assert_eq!(q, vec![
            Oklab { l: 0.5, a: 1.0, b: 0.25 },
            Oklab { l: 0.5, a: 0.0, b: 0.25 },
        ]);
    }

    #[test]
    fn test_kmeans_rect() {
        let points = [
            Oklab { l: 0.5, a: 0.0, b: 0.0 },
            Oklab { l: 0.5, a: 0.0, b: 0.5 },
            Oklab { l: 0.5, a: 1.0, b: 0.0 },
            Oklab { l: 0.5, a: 1.0, b: 0.5 },
        ];
        let k = 2usize;
        let mut rng = Xoshiro256StarStar::seed_from_u64(0xdeadbeef);
        let initial = OkImage::kmeans_initial(&points, k, &mut rng);
        assert_eq!(vec![
            Oklab { l: 0.5, a: 1.0, b: 0.5 },
            Oklab { l: 0.5, a: 0.0, b: 0.0 },
        ], initial);
    }

    #[test]
    fn test_kmeans_all_the_same() {
        let points = [
            Oklab { l: 0.5, a: 0.0, b: 0.5 },
            Oklab { l: 0.5, a: 0.0, b: 0.5 },
            Oklab { l: 0.5, a: 0.0, b: 0.5 },
            Oklab { l: 0.5, a: 0.0, b: 0.5 },
        ];
        let k = 2usize;
        let mut rng = Xoshiro256StarStar::seed_from_u64(0xdeadbeef);
        let initial = OkImage::kmeans_initial(&points, k, &mut rng);
        assert_eq!(vec![
            Oklab { l: 0.5, a: 0.0, b: 0.5 },
        ], initial);
    }

    #[test]
    fn test_kmeans_k_larger_than_points() {
        let points = [
            Oklab { l: 0.5, a: 0.0, b: 0.5 },
            Oklab { l: 0.5, a: 1.0, b: 0.5 },
        ];
        let k = 3usize;
        let mut rng = Xoshiro256StarStar::seed_from_u64(0xdeadbeef);
        let initial = OkImage::kmeans_initial(&points, k, &mut rng);
        assert_eq!(vec![
            Oklab { l: 0.5, a: 1.0, b: 0.5 },
            Oklab { l: 0.5, a: 0.0, b: 0.5 },
        ], initial);
    }
}
