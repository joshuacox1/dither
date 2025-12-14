use std::borrow::Cow;
use std::io::{Read, Write};
use std::io;

use super::{Vec2D, Rgb888};

/// An image consisting of one or more frames, each a 2D grid of pixels
/// (width and height both positive).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Image<PixelType> {
    frames: Vec<Vec2D<PixelType>>,
}

/// Information associated with an `ImageError`.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct ImageErrorInfo {
    /// The index of the offending frame.
    pub frame_idx: usize,
    /// The dimensions of the offending frame.
    pub dims: (usize, usize),
}

/// An error encountered when creating or manipulating an `Image`.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum ImageError {
    /// No image frames were given. An image must have at least one frame.
    NoFrames,
    /// One of the image frames had zero width or height.
    ZeroDims(ImageErrorInfo),
    /// Two image frames had different dimensions.
    InconsistentDims(ImageErrorInfo, ImageErrorInfo),
}

// TODO: custom imagecreationerror type
// TODO: error for empty vec width or height


impl<PixelType> Image<PixelType> {
    /// Creates an image from a single frame (hence not animated).
    pub fn single(image: Vec2D<PixelType>) -> Result<Self, ImageError> {
        let (w,h) = image.dims();
        if w == 0 || h == 0 {
            Err(ImageError::ZeroDims(ImageErrorInfo {
                frame_idx: 0,
                dims: (w, h),
            }))
        } else {
            Ok(Self { frames: vec![image] })
        }
    }

    /// Creates an image from multiple frames. Fails if there
    /// are no frames or if each frame has different dimensions.
    pub fn from_frames(frames: Vec<Vec2D<PixelType>>) -> Result<Self, ImageError> {
        if frames.len() == 0 {
            return Err(ImageError::NoFrames);
        }

        let (w, h) = frames[0].dims();
        if w == 0 || h == 0 {
            return Err(ImageError::ZeroDims(ImageErrorInfo {
                frame_idx: 0,
                dims: (w, h),
            }));
        }
        for i in 1..frames.len() {
            let img = &frames[i];
            let (w2, h2) = img.dims();
            if w != w2 || h != h2 {
                let info1 = ImageErrorInfo { frame_idx: 0, dims: (w, h) };
                let info2 = ImageErrorInfo { frame_idx: i, dims: (w2, h2) };
                return Err(ImageError::InconsistentDims(info1, info2));
            }
        }

        Ok(Self { frames })
    }

    /// The number of frames in the image. Will be at least one.
    pub fn num_frames(&self) -> usize { self.frames.len() }

    /// Gets the frame at the given zero-based index, if it exists.
    pub fn get_frame(&self, idx: usize) -> Option<&Vec2D<PixelType>> {
        self.frames.get(idx)
    }

    /// Maps an image transformation function across all image frames.
    pub fn map_frames<F, P>(&self, f: F) -> Result<Image<P>, ImageError>
    where F: Fn(&Vec2D<PixelType>) -> Vec2D<P> {
        Image::from_frames(self.frames.iter().map(f).collect::<Vec<_>>())
    }

    /// Maps an image transformation function across all pixels in all
    /// frames.
    pub fn map_pixels<F, P>(&self, f: F) -> Result<Image<P>, ImageError>
    where F: Fn(&PixelType) -> P {
        self.map_frames(|g| g.map(|p| f(p)))
    }

    /// Maps an image transformation function across all pixels in all
    /// frames.
    pub fn map_pixels_idx<F, P>(&self, f: F) -> Result<Image<P>, ImageError>
    where F: Fn((usize, usize), &PixelType) -> P {
        self.map_frames(|g| g.map_idx(|xy, p| f(xy, p)))
    }

    /// The dimensions of the image.
    pub fn dims(&self) -> (usize, usize) {
        self.frames[0].dims()
    }
}

impl Image<Rgb888> {
    /// Writes an image with pixel type `Rgb888` to the PNG format.
    pub fn write_png<W: Write>(writer: W) -> io::Result<()> {
        todo!()
    }

    /// Reads an image from the given reader
    pub fn read_png<R: Read>(reader: R) -> io::Result<Self> {
        todo!()
    }
}

impl Image<u8> {
    /// Writes an indexed image with the given palette to the PNG format.
    /// The PNG will have indexed colour type with as few bits as possible
    /// (1-bit if <= 2 colours, 2-bit if <= 4 colours, 4-bit if <= 16 colours,
    /// 8-bit if <= 256 colours). The palette size must be between 1 and 256.
    /// Throws an error if a pixel is OOB? Or just writes an invalid PNG?
    pub fn write_png<W: Write>(writer: W, palette: &[Rgb888]) -> io::Result<()> {
        todo!()
    }
}



