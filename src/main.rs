use std::io::{Write, BufWriter};
use std::fs::File;
use image::ImageReader;
use std::collections::HashSet;

use dither::{Oklabr, Image, Rgb888, NoDither, KnollDither,
    Palette, Dither, ErrorDiffusion, Quantise, KmeansClustering};

fn main() {
    // let imgpath = "/Users/joshua/OneDrive/projects/dither/cat_undithered.png";
    // let catto = Image::<Rgb888>::read_img_sloppy(imgpath);
    // let pixels = catto.get_frame(0).unwrap().data_1d().iter()
    //     .map(|&rgb| rgb).collect::<HashSet<_>>();
    // println!("{}", pixels.len());
    // for p in pixels {
    //     println!("{}", p.to_str());
    // }

    // let cat_palette = Palette::new([
    //     "#ffff99",
    //     "#996666",
    //     "#999966",
    //     "#ffcccc",
    //     "#333300",
    //     "#cccccc",
    //     "#993333",
    //     "#666666",
    //     "#336666",
    //     "#ffcc99",
    //     "#999999",
    //     "#006633",
    //     "#660000",
    //     "#cccc99",
    //     "#333333",
    //     "#99cc66",
    //     "#339999",
    //     "#99cccc",
    //     "#cc9999",
    //     "#000000",
    //     "#663300",
    //     "#336633",
    //     "#663333",
    //     "#006666",
    //     "#339966",
    //     "#ccffff",
    //     "#cc9966",
    //     "#669999",
    //     "#cccc66",
    //     "#996633",
    //     "#99cc99",
    //     "#999933",
    //     "#666633",
    //     "#ffffff",
    //     "#330000",
    //     "#ffffcc",
    //     "#669966",
    //     "#1e1612",
    //     "#2e2522",
    //     "#3c3330",
    //     "#49423c",
    //     "#564f49",
    //     "#645d56",
    //     "#716b65",
    //     "#417f63",
    //     "#837d78",
    //     "#91908c",
    //     "#a0a29d",
    //     "#b2b3ae",
    //     "#c4c4bf",
    //     "#cfd4d3",
    //     "#e0e3e3",
    //     "#f6f8f8",
    // ].into_iter()
    //     .map(|s| Rgb888::from_str(s).unwrap())
    //     .collect::<Vec<_>>()).unwrap().map(|rgb| rgb.to_oklabr());

    let imgpath = "/Users/joshua/OneDrive/projects/dither/blossom.jpg";
    let img = Image::<Rgb888>::read_img_sloppy(imgpath)
        .map_pixels(|rgb| rgb.to_oklabr());

    // let quantise_impl = cat_palette;
    let quantise_impl = KmeansClustering::new(16, 0xdeadbeef).unwrap();
    let palette = quantise_impl.quantise(&img);
    let mut palette = palette.snap_to_unique_srgb888();
    palette.sort();


    // let palette_rgb888 = Palette::new([
    //     // "#000000", "#ffffff",
    //     "#080000", "#201A0B", "#432817", "#492910",
    //     "#234309", "#5D4F1E", "#9c6b20", "#a9220f",
    //     "#2b347c", "#2b7409", "#d0ca40", "#e8a077",
    //     "#6a94ab", "#d5c4b3", "#fce76e", "#fcfae2",
    // ].into_iter()
    //     .map(|s| Rgb888::from_str(s).unwrap())
    //     .collect::<Vec<_>>()
    // ).unwrap();
    // let palette_oklabr = palette_rgb888.map(|rgb| rgb.to_oklabr());

    // let dither_impl = NoDither;
    // let dither_impl = KnollDither::new(8, 1.0/3.0).unwrap();
    let dither_impl = ErrorDiffusion::FloydSteinberg;
    let dithered_img = dither_impl.dither(&img, &palette);

    let file = File::create("output.png").unwrap();
    let ref mut w = BufWriter::new(file);
    let palette_rgb888 = palette.map(|c| c.to_rgb888());
    dithered_img.write_png(w, &palette_rgb888).unwrap();


    for rgb in palette_rgb888.iter() {
        println!("\"{}\",", rgb.to_str());
    }
}
