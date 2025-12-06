use std::io;
use image::ImageReader;

use dither::{Oklabr, PxImage, DitherOptions, PaletteOptions, Problem};

fn main() {
    // //"C:\\Users\\Joshua\\OneDrive\\pictures\\floor_goban.jpg",
    let mut img = PxImage::from_path(
        //"C:\\Users\\Joshua\\OneDrive\\projects\\dither\\test4x4.png");
         "C:\\Users\\Joshua\\OneDrive\\projects\\dither\\david.png");
        //"C:\\Users\\Joshua\\OneDrive\\pictures\\account\\ken_amada_winter_small.png");
        //"C:\\Users\\Joshua\\OneDrive\\pictures\\floor_goban.jpg");
        //"C:\\Users\\Joshua\\OneDrive\\pictures\\the_music_lesson_hq.jpg");
    println!("Image loaded in");

    fn discretize(x: f32) -> u8 {
        (x * 255.0).clamp(0.0, 255.0).round() as u8
    }


    let blue = Oklabr::from_srgb_888_str("#0000ff").unwrap();
    let cyan = Oklabr::from_srgb_888_str("#00ffff").unwrap();
    println!("{blue:?}, {cyan:?}");
    let halfway = (blue + cyan) / 2.0;
    println!("halfway: {halfway:?}, {:?}", halfway.to_srgb());
    // return;


    println!("Computing palette");
    // let pal = img.palette(2, 0xbeefdead);
    let pal = [Oklabr::WHITE, Oklabr::BLACK];

    // let pal = [
    //     Oklabr { l: 0.24, a: 0.01, b: -0.04 },
    //     Oklabr { l: 0.72, a: -0.18, b: 0.12 },
    //     Oklabr { l: 0.68, a: 0.08, b: 0.13 },
    //     Oklabr { l: 0.56, a: 0.19, b: 0.08 },
    //     Oklabr { l: 0.53, a: 0.08, b: -0.12 },
    //     Oklabr { l: 0.64, a: -0.06, b: -0.13 },
    //     Oklabr { l: 0.99, a: 0.0, b: 0.01 },
    // ];

    // let pal = [
    //     Oklabr { l: 0.084, a: 0.03, b: 0.02 },
    //     Oklabr { l: 0.22, a: 0.0, b: 0.03 },
    //     Oklabr { l: 0.31, a: 0.03, b: 0.04 },
    //     Oklabr { l: 0.32, a: 0.03, b: 0.05 },

    //     Oklabr { l: 0.35, a: -0.06, b: 0.07 },
    //     Oklabr { l: 0.43, a: 0.0, b: 0.07 },
    //     Oklabr { l: 0.57, a: 0.03, b: 0.1 },
    //     Oklabr { l: 0.48, a: 0.15, b: 0.09 },

    //     Oklabr { l: 0.36, a: 0.01, b: -0.12 },
    //     Oklabr { l: 0.49, a: -0.11, b: 0.1 },
    //     Oklabr { l: 0.82, a: -0.04, b: 0.15 },
    //     Oklabr { l: 0.77, a: 0.07, b: 0.08 },

    //     Oklabr { l: 0.64, a: -0.03, b: -0.05 },
    //     Oklabr { l: 0.83, a: 0.01, b: 0.03 },
    //     Oklabr { l: 0.92, a: -0.02, b: 0.14 },
    //     Oklabr { l: 0.98, a: -0.01, b: 0.03 },
    // ];
    let pal = pal.into_iter().map(|c| (c.to_srgb().to_srgb888(), c)).collect::<Vec<_>>();


    for (crgb, c) in &pal {
        let [r,g,b] = crgb;
        println!("#{r:02x}{g:02x}{b:02x} | {c}");
    }


    println!("Palette computed.");
    // for i in 0..100 {
    //     let pal = img.palette(4, 0xbeefdead + i);
    //     for c in &pal {
    //         let q = c.to_linear_srgb().to_gamma();
    //         let r = discretize(q.r);
    //         let g = discretize(q.g);
    //         let b = discretize(q.b);
    //         print!("{i:02}: #{r:02x}{g:02x}{b:02x} | ");
    //     }
    //     // print!("{i:02}: {:?} {:?}", pal[0], pal[1]);
    //     println!("");
    // }


    let dither_options = DitherOptions {
        num_samples: 8,
        matrix_size: 4,
        strength: 0.3,
    };

    let palette_indices = img.knoll_dither(&dither_options, &pal);
    println!("Dithering done");
    dither::write_indexed_png("result.png", img.dims(), &pal, &palette_indices);
    // img.save_as_srgb_png("meep2.png");
}

/* inputs:
- aa





*/
