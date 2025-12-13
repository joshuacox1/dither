use std::ops::{Index, IndexMut};

/// A two-dimensional non-resizable array with size determined at runtime.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Vec2D<T> {
    data: Vec<T>,
    dims: (usize, usize),
}

// todo: better errors
// implement indexing, map, map where the index is a part, etc.
// todo: par_iter, or better yet, GPU-stuff
impl<T> Vec2D<T> {
    /// Creates a `Vec2D` from a 1-dimensional buffer of increasing rows.
    fn from_1d(data: Vec<T>, dims: (usize, usize)) -> Result<Self, ()> {
        if dims.0 * dims.1 == data.len() {
            Ok(Self { data, dims })
        } else {
            Err(())
        }
    }

    fn valid_coords(&self, coords: &(usize, usize)) -> bool {
        let (i,j) = coords; let (w,h) = self.dims;
        *i < w && *j < h
    }

    fn coords_2d_to_1d(&self, coords: &(usize, usize)) -> usize {
        let (i,j) = *coords;
        self.dims.0 * j + i
    }

    /// Maps over the elements of a `Vec2D`, preserving dimensions.
    pub fn map<F, U>(&self, f: F) -> Vec2D<U>
    where F: Fn(&T) -> U {
        Vec2D {
            data: self.data.iter().map(f).collect::<Vec<_>>(),
            dims: self.dims,
        }
    }

    /// Maps over the elements of a `Vec2D` in-place.
    pub fn map_inplace<F>(&mut self, f: F)
    where F: Fn(&mut T) {
        for v in self.data.iter_mut() {
            f(v);
        }
    }

    /// The dimensions of the `Vec2D`, i.e. (row count, column count).
    pub fn dims(&self) -> (usize, usize) { self.dims }
}

impl<T> Index<(usize, usize)> for Vec2D<T> {
    type Output = T;

    fn index(&self, coords: (usize, usize)) -> &Self::Output {
        if self.valid_coords(&coords) {
            &self.data[self.coords_2d_to_1d(&coords)]
        } else {
            panic!("Invalid indices {coords:?} for dimensions {:?}", self.dims)
        }
    }
}

impl<T> IndexMut<(usize, usize)> for Vec2D<T> {
    fn index_mut(&mut self, coords: (usize, usize)) -> &mut Self::Output {
        if self.valid_coords(&coords) {
            let idx = self.coords_2d_to_1d(&coords);
            &mut self.data[idx]
        } else {
            // todo combine
            panic!("Invalid indices {coords:?} for dimensions {:?}", self.dims)
        }
    }
}


#[cfg(test)]
mod test {
    use super::*;

    fn test_vec2d() -> Vec2D<u8> {
        let data = (1u8..=12).collect::<Vec<_>>();
        let dims = (4,3);
        Vec2D::from_1d(data, dims).unwrap()
    }

    #[test]
    fn test_indexing() {
        let mut vec2d = test_vec2d();
        assert_eq!(vec2d[(2,1)], 7);

        vec2d[(1,2)] = 18;
        assert_eq!(vec2d[(1,2)], 18);
    }

    #[test] #[should_panic]
    fn test_index_fail() {
        let vec2d = test_vec2d();
        let _ = vec2d[(1,5)];
    }

    #[test]
    fn test_map() {
        let vec2d = test_vec2d();
        let vec2dnew = vec2d.map(|&u| (u as i16) * 10);

        let expected: [[i16; 4]; 3] = [
            [10, 20, 30, 40],
            [50, 60, 70, 80],
            [90, 100, 110, 120],
        ];

        for j in 0..3 {
            for i in 0..4 {
                let idx = (i,j);
                assert_eq!(vec2dnew[idx], expected[j][i]);
            }
        }
    }

    #[test]
    fn test_map_inplace() {
        let mut vec2d = test_vec2d();
        vec2d.map_inplace(|u| { *u *= 10; });

        let expected: [[u8; 4]; 3] = [
            [10, 20, 30, 40],
            [50, 60, 70, 80],
            [90, 100, 110, 120],
        ];

        for j in 0..3 {
            for i in 0..4 {
                let idx = (i,j);
                assert_eq!(vec2d[idx], expected[j][i]);
            }
        }
    }
}
