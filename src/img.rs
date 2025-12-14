use std::borrow::Cow;
use std::io::{Read, Write};
use std::io;

use super::{Vec2D, Rgb888};

use png::{Encoder, EncodingError,
    ColorType, SrgbRenderingIntent, DeflateCompression,
    BitDepth};

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
    ///
    /// NOTE: This currently only writes out the first frame. APNG is not supported.
    /// Further work is to write a valid APNG with all the frames.
    /// Further further work is to optimise that APNG with tricks.
    pub fn write_png<W: Write>(&self, writer: W, palette: &[Rgb888]) -> Result<(), EncodingError> {
        // todo fix
        // pub fn write_indexed_png(
        //     path: &str,
        //     dims: (u32, u32),
        //     palette: &[([u8; 3], Oklabr)],
        //     data: &[u8],
        // ) {
        // let file = File::create(path).unwrap();
        // let ref mut w = BufWriter::new(file);
        let (w, h) = self.dims();

        let mut e = png::Encoder::new(writer, w as u32, h as u32);
        e.set_source_srgb(SrgbRenderingIntent::RelativeColorimetric);
        e.set_color(ColorType::Indexed);
        // the highest level of compression
        e.set_deflate_compression(DeflateCompression::Level(9));

        let lp = palette.len();
        let bit_depth = match lp {
            _ if lp <= 2 => BitDepth::One,
            _ if lp <= 4 => BitDepth::Two,
            _ if lp <= 16 => BitDepth::Four,
            _ if lp <= 256 => BitDepth::Eight,
            _ => panic!("Palette with {lp} > 256 colour entries was passed in"),
        };
        e.set_depth(bit_depth);

        // By the PNG spec, the palette chunk does not need to contain padding
        // entries up to 256 or whatever. It is fine as long as all the pixels
        // in the image have a valid palette index.
        // The palette must be 8-bit RGB.
        let mut palette_bytes = Vec::<u8>::with_capacity(3*lp);
        for p in palette {
            palette_bytes.push(p.r);
            palette_bytes.push(p.g);
            palette_bytes.push(p.b);
        }
        e.set_palette(palette_bytes);
        let mut writer = e.write_header()?;

        let data = &self.frames[0];

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
                        let idx = 8*i;
                        let mut b = data[(idx,j)] << 7;
                        if let Some(p) = data.get((idx+1,j)) { b |= p << 6; }
                        if let Some(p) = data.get((idx+2,j)) { b |= p << 5; }
                        if let Some(p) = data.get((idx+3,j)) { b |= p << 4; }
                        if let Some(p) = data.get((idx+4,j)) { b |= p << 3; }
                        if let Some(p) = data.get((idx+5,j)) { b |= p << 2; }
                        if let Some(p) = data.get((idx+6,j)) { b |= p << 1; }
                        if let Some(p) = data.get((idx+7,j)) { b |= p; }
                        z.push(b);
                    }
                }
                writer.write_image_data(&z)?;
            },
            BitDepth::Two => {
                let w_bytes = w.div_ceil(4);
                let data_size = 4*w_bytes * h;
                let mut z = Vec::with_capacity(data_size);
                for j in 0..h {
                    for i in 0..w_bytes {
                        let idx = 4*i;
                        let mut b = data[(idx,j)] << 6;
                        if let Some(p) = data.get((idx+1,j)) { b |= p << 4; }
                        if let Some(p) = data.get((idx+2,j)) { b |= p << 2; }
                        if let Some(p) = data.get((idx+3,j)) { b |= p; }
                        z.push(b);
                    }
                }
                writer.write_image_data(&z)?;
            },
            BitDepth::Four => {
                let w_bytes = w.div_ceil(2);
                let data_size = 2*w_bytes * h;
                let mut z = Vec::with_capacity(data_size);
                for j in 0..h {
                    for i in 0..w_bytes {
                        let idx = 2*i;
                        let mut b = data[(idx,j)] << 4;
                        if let Some(p) = data.get((idx+1,j)) { b |= p; }
                        z.push(b);
                    }
                }
                writer.write_image_data(&z)?;
            },
            BitDepth::Eight => {
                writer.write_image_data(data.data_1d())?;
            },
            BitDepth::Sixteen => unsafe { std::hint::unreachable_unchecked() },
        }

        Ok(())
    }
}






