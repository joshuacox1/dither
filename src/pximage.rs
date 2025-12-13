//! Wrappers to convert to and from PNG.
//! It would be nice to support APNG at some point. But for now, just PNG.

use super::{Oklabr, Rgb888, DitherOptions};
use image::{ImageReader, Rgb, ImageBuffer};

use std::ops::{Index, IndexMut};
use std::collections::HashSet;

use rand::{Rng, SeedableRng};
use rand_xoshiro::{Xoshiro256StarStar};

// /// A 2D array of stuffs.
pub struct PxImage {
    data: Vec<Oklabr>,
    dims: (u32, u32),
}

impl PxImage {
    pub fn from_path(s: &str) -> Self {
        // todo: remove these unwraps lol.
        let img = ImageReader::open(s)
            .unwrap().decode().unwrap()
            .into_rgb8();
        let w = img.width();
        let h = img.height();
        let data = img.pixels()
            .map(|rgb| Rgb888 { r: rgb[0], g: rgb[1], b: rgb[2] }.to_oklabr())
            .collect::<Vec<_>>();
        Self::new(data, (w, h))
    }

    // TODO: delete
    // pub fn save_as_srgb_png(&self, s: &str) {
    //     // rayon?!
    //     // TODO: also, indexed colour?!
    //     let (w,h) = self.dims;
    //     let img = ImageBuffer::from_fn(
    //         w, h,
    //         |i,j| {
    //             let (i,j) = (i as usize, j as usize);
    //             let px = self[(i,j)];
    //             // TODO: make this use Oklabr methods.
    //             let srgb = px.to_srgb();
    //             Rgb::<u8>([
    //                 (srgb.r * 255.0).clamp(0.0, 255.0) as u8,
    //                 (srgb.g * 255.0).clamp(0.0, 255.0) as u8,
    //                 (srgb.b * 255.0).clamp(0.0, 255.0) as u8,
    //             ])
    //         });
    //     img.save(s).unwrap();
    // }

    fn new(data: Vec<Oklabr>, dims: (u32, u32)) -> Self {
        assert!(data.len() == dims.0 as usize * dims.1 as usize);
        assert!(data.len() > 0);
        Self { data, dims }
    }

    /// The dimensions of the image (width, height). Both will be
    /// at least `1`.
    pub fn dims(&self) -> (u32, u32) { self.dims }

    /// Perform Knoll dithering on the image.
    /// For now: 8-bit indexed colour only. Consider other bits later.
    pub fn knoll_dither(&self,
        options: &DitherOptions,
        palette: &[([u8; 3], Oklabr)],
    ) -> Vec<u8> {
        let w = self.dims.0 as usize; let h = self.dims.1 as usize;
        let n = options.num_samples as usize;
        let b = options.matrix_size as usize;
        // todo: delete this
        let palette = palette.iter().map(|(_,c)| *c).collect::<Vec<_>>();

        let mut result = Vec::with_capacity(w*h); // divide by 8/bit depth
        // todo: maybeuninit
        let mut candidates = vec![usize::MAX; n];
        let bayer_interped = bayer_1d(b, n);
        for j in 0..h {
            for i in 0..w {
                let pixel = self[(i,j)];
                Self::knoll_dither_inner(&mut candidates, &palette,
                    &pixel, n, options.strength);

                // todo: reconsider bayer interping.
                let bayer_idx = bayer_interped[(j % b)*b + (i % b)];
                let palette_idx = candidates[bayer_idx];
                result.push(palette_idx as u8);
            }
        }

        result
    }

    // Overwrites candidates.
    pub fn knoll_dither_inner(
        candidates: &mut [usize],
        palette: &[Oklabr],
        pixel: &Oklabr,
        n: usize,
        strength: f32)
    {
        let mut error = Oklabr { l: 0.0, a: 0.0, b: 0.0 };
        for it in 0..n {
            let attempt = *pixel + error * strength;
            let (argclosest, closest) = palette.iter().enumerate()
                .min_by(|(i,p),(j,q)| p.sq_dist(&attempt)
                    .total_cmp(&q.sq_dist(&attempt)))
                .unwrap();
            error += *pixel - *closest;
            candidates[it] = argclosest;
        }

        // The default order on Oklabr colours is increasing by luminance,
        // which is what we need here.
        candidates.sort_unstable_by_key(|&i| palette[i]);
    }

    /// Obtain a representative palette from the image of at most k colours.
    /// k must be positive. The palette will generally consist of
    /// exactly k colours unless the original image has fewer colours.
    /// Performs kmeans clustering with kmeans++ initialisation.
    /// (todo: should I transpose from Oklabr back into eg rgb888,
    /// then back again?)
    /// The palette will have at least one colour in it.
    pub fn palette(&self, k: usize, random_seed: u64) -> Vec<(Rgb888, Oklabr)> {
        let mut rng = Xoshiro256StarStar::seed_from_u64(random_seed);

        let centroids = Self::naive_k_means(&self.data, k, &mut rng);

        // Discretize to RGB888. Ensure uniqueness (though in practice
        // it should almost never cause centroids to be merged this way).
        let centroids_rgb = centroids.into_iter()
            .map(|c| c.to_rgb888())
            .collect::<HashSet<_>>();

        // Map back to Oklabr colours and sort. (The Oklabr ordering
        // sort by luminance first, which is what we are mainly
        // interested in.)
        let mut centroids = centroids_rgb.into_iter()
            .map(|c| (c, c.to_oklabr()))
            .collect::<Vec<_>>();
        centroids.sort_unstable_by_key(|(_,c)| *c);

        centroids
    }

    // Naive k-means on Oklabr colours.
    // Precondition: points is non-empty and k > 0.
    fn naive_k_means(
        points: &[Oklabr],
        k: usize,
        rng: &mut impl Rng,
    ) -> Vec<Oklabr> {
        let mut centroids = Self::kmeans_initial(points, k, rng);
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
            println!("kmeans step {counter}...");

            if converged {
                break;
            }

            // Recalculate the centroids of each cluster.
            for c in centroids.iter_mut() {
                *c = Oklabr { l: 0.0, a: 0.0, b: 0.0 };
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
            counter += 1;
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
        points: &[Oklabr],
        k: usize,
        rng: &mut impl Rng,
    ) -> Vec<Oklabr> {
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

impl Index<(usize, usize)> for PxImage {
    type Output = Oklabr;

    fn index(&self, ij: (usize, usize)) -> &Self::Output {
        let (i, j) = ij;
        let width = self.dims.0 as usize; let height = self.dims.1 as usize;
        if i < width && j < height {
            &self.data[width*j + i]
        } else {
            panic!("Failed to index img with dims {:?} with indices {ij:?}",
                self.dims);
        }
    }
}

impl IndexMut<(usize, usize)> for PxImage {
    fn index_mut(&mut self, ij: (usize, usize)) -> &mut Self::Output {
        let (i, j) = ij;
        let width = self.dims.0 as usize; let height = self.dims.1 as usize;
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






fn bayer_1d(matrix_size: usize, num_samples: usize) -> Vec<usize> {
    match matrix_size {
        4 => interp_bayer(BAYER_4, num_samples),
        8 => interp_bayer(BAYER_8, num_samples),
        16 => interp_bayer(BAYER_16, num_samples),
        _ => panic!("Invalid Bayer matrix size {matrix_size}."),
    }
}

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


/// Optimal sorting network (in both comparisons and layers) for 8
/// elements.
fn sort_8<T: PartialOrd>(arr: &mut [T]) {
    let mut cswap = |i,j| if arr[i] > arr[j] { arr.swap(i,j); };
    cswap(0,2); cswap(1,3); cswap(4,6); cswap(5,7);
    cswap(0,4); cswap(1,5); cswap(2,6); cswap(3,7);
    cswap(0,1); cswap(2,3); cswap(4,5); cswap(6,7);
    cswap(2,4); cswap(3,5);
    cswap(1,4); cswap(3,6);
    cswap(1,2); cswap(3,4); cswap(5,6);
}


use std::fs::File;
use std::io::BufWriter;
use png::{BitDepth, SrgbRenderingIntent, ColorType, DeflateCompression};

pub fn write_indexed_png(
    path: &str,
    dims: (u32, u32),
    palette: &[([u8; 3], Oklabr)],
    data: &[u8],
) {
    let file = File::create(path).unwrap();
    let ref mut w = BufWriter::new(file);

    let mut e = png::Encoder::new(w, dims.0, dims.1);
    e.set_source_srgb(SrgbRenderingIntent::RelativeColorimetric);
    e.set_color(ColorType::Indexed);
    e.set_deflate_compression(DeflateCompression::Level(9));

    let lp = palette.len();
    assert!(1 <= lp && lp <= 256);
    let bit_depth = if lp <= 2 {
        BitDepth::One
    } else if lp <= 4 {
        BitDepth::Two
    } else if lp <= 16 {
        BitDepth::Four
    } else {
        BitDepth::Eight
    };
    e.set_depth(bit_depth);

    // By the PNG spec, the palette chunk does not need to contain padding
    // entries up to 256 or whatever. It is fine as long as all the pixels
    // in the image have a valid palette index.
    // The palette must be 8-bit RGB.
    let mut palette_bytes = Vec::<u8>::with_capacity(3*lp);
    for p in palette {
        palette_bytes.push(p.0[0]);
        palette_bytes.push(p.0[1]);
        palette_bytes.push(p.0[2]);
    }
    e.set_palette(palette_bytes);
    let mut writer = e.write_header().unwrap();

    let d = data.len();
    let w = dims.0 as usize; let h = dims.1 as usize;

    // Per the PNG spec, for < 8-bit depth images, rows always start
    // at a byte boundary (so there is padding of unspecified value
    // if the width is not a clean multiple). Within a byte, pixels are
    // laid out where left to right = highest bits to lowest bits.
    match bit_depth {
        BitDepth::One => {
            let w_bytes = w.div_ceil(8);
            let data_size = 8*w_bytes * h;
            let mut z = Vec::with_capacity(data_size);
            for j in 0..h {
                for i in 0..w_bytes {
                    let idx = j*w + 8*i;
                    let mut b = data[idx] << 7;
                    if let Some(p) = data.get(idx+1) { b |= p << 6; }
                    if let Some(p) = data.get(idx+2) { b |= p << 5; }
                    if let Some(p) = data.get(idx+3) { b |= p << 4; }
                    if let Some(p) = data.get(idx+4) { b |= p << 3; }
                    if let Some(p) = data.get(idx+5) { b |= p << 2; }
                    if let Some(p) = data.get(idx+6) { b |= p << 1; }
                    if let Some(p) = data.get(idx+7) { b |= p; }
                    z.push(b);
                }
            }
            writer.write_image_data(&z).unwrap();
        },
        BitDepth::Two => {
            let w_bytes = w.div_ceil(4);
            let data_size = 4*w_bytes * h;
            let mut z = Vec::with_capacity(data_size);
            for j in 0..h {
                for i in 0..w_bytes {
                    let idx = j*w + 4*i;
                    let mut b = data[idx] << 6;
                    if let Some(p) = data.get(idx+1) { b |= p << 4; }
                    if let Some(p) = data.get(idx+2) { b |= p << 2; }
                    if let Some(p) = data.get(idx+3) { b |= p; }
                    z.push(b);
                }
            }
            writer.write_image_data(&z).unwrap();
        },
        BitDepth::Four => {
            let w_bytes = w.div_ceil(2);
            let data_size = 2*w_bytes * h;
            let mut z = Vec::with_capacity(data_size);
            for j in 0..h {
                for i in 0..w_bytes {
                    let idx = j*w + 2*i;
                    let mut b = data[idx] << 4;
                    if let Some(p) = data.get(idx+1) { b |= p; }
                    z.push(b);
                }
            }
            writer.write_image_data(&z).unwrap();
        },
        BitDepth::Eight => {
            writer.write_image_data(data).unwrap();
        },
        BitDepth::Sixteen => unsafe { std::hint::unreachable_unchecked() },
    }

    // TODO: some of these unwraps should not be unwraps.
    // only the write_header().unwrap is legitimate as that should
    // never fail as I've inputted valid values
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_kmeans_rect_full() {
        let img = PxImage {
            data: vec![
                Oklabr { l: 0.5, a: 0.0, b: 0.0 },
                Oklabr { l: 0.5, a: 0.0, b: 0.5 },
                Oklabr { l: 0.5, a: 1.0, b: 0.0 },
                Oklabr { l: 0.5, a: 1.0, b: 0.5 },
            ],
            dims: (2,2),
        };
        let k = 2usize;
        let seed = 0xdeadbeef;
        let q = img.palette(k, seed).into_iter().map(|(_,x)| x).collect::<Vec<_>>();
        assert_eq!(q, vec![
            Oklabr { l: 0.5, a: 1.0, b: 0.25 },
            Oklabr { l: 0.5, a: 0.0, b: 0.25 },
        ]);
    }

    #[test]
    fn test_kmeans_rect() {
        let points = [
            Oklabr { l: 0.5, a: 0.0, b: 0.0 },
            Oklabr { l: 0.5, a: 0.0, b: 0.5 },
            Oklabr { l: 0.5, a: 1.0, b: 0.0 },
            Oklabr { l: 0.5, a: 1.0, b: 0.5 },
        ];
        let k = 2usize;
        let mut rng = Xoshiro256StarStar::seed_from_u64(0xdeadbeef);
        let initial = PxImage::kmeans_initial(&points, k, &mut rng);
        assert_eq!(vec![
            Oklabr { l: 0.5, a: 1.0, b: 0.5 },
            Oklabr { l: 0.5, a: 0.0, b: 0.0 },
        ], initial);
    }

    #[test]
    fn test_kmeans_all_the_same() {
        let points = [
            Oklabr { l: 0.5, a: 0.0, b: 0.5 },
            Oklabr { l: 0.5, a: 0.0, b: 0.5 },
            Oklabr { l: 0.5, a: 0.0, b: 0.5 },
            Oklabr { l: 0.5, a: 0.0, b: 0.5 },
        ];
        let k = 2usize;
        let mut rng = Xoshiro256StarStar::seed_from_u64(0xdeadbeef);
        let initial = PxImage::kmeans_initial(&points, k, &mut rng);
        assert_eq!(vec![
            Oklabr { l: 0.5, a: 0.0, b: 0.5 },
        ], initial);
    }

    #[test]
    fn test_kmeans_k_larger_than_points() {
        let points = [
            Oklabr { l: 0.5, a: 0.0, b: 0.5 },
            Oklabr { l: 0.5, a: 1.0, b: 0.5 },
        ];
        let k = 3usize;
        let mut rng = Xoshiro256StarStar::seed_from_u64(0xdeadbeef);
        let initial = PxImage::kmeans_initial(&points, k, &mut rng);
        assert_eq!(vec![
            Oklabr { l: 0.5, a: 1.0, b: 0.5 },
            Oklabr { l: 0.5, a: 0.0, b: 0.5 },
        ], initial);
    }

    // #[test]
    // fn test_knoll_dither_inner() {
    //     let n = 8;
    //     let mut candidates = [0; 8];
    //     let pixel = Oklabr { l: 0.0, a: 0.0, b: 0.0 };
    //     let palette = [
    //         Oklabr { l: -1.0, a: 0.0, b: 0.0 },
    //         Oklabr { l: 2.0, a: 0.0, b: 0.0 },
    //     ];
    //     let strength = 0.5;

    //     PxImage::knoll_dither_inner(&mut candidates, &palette, &pixel,
    //         n, strength);

    //     assert_eq!(vec![Oklabr::BLACK], candidates);
    // }

    // #[test]
    // fn test_sort_8() {
    //     for i1 in 0..8 {
    //         for i2 in 0..8 {
    //             for i3 in 0..8 {
    //                 for i4 in 0..8 {
    //                     for i5 in 0..8 {
    //                         for i6 in 0..8 {
    //                             for i7 in 0..8 {
    //                                 for i8 in 0..8 {
    //                                     let mut v = [i1,i2,i3,i4,i5,i6,i7,i8];
    //                                     sort_8(&mut v);
    //                                     for j in 0..7 {
    //                                         assert!(!(v[j] > v[j+1]));
    //                                     }
    //                                 }
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
}
