// //! The RNG. Simple 64-bit xorshift.
// //! TODO: do something slightly better? Like the 192 bit thing.

// use std::time::SystemTime;

// pub struct Rng {
//     state: u64,
// }

// impl Rng {
//     /// Creates a seed using the system clock.
//     pub fn seed_from_system_clock() -> u64 {
//         // The size of the Monster group mod 2**64
//         let random_constant = 10101925807912910848u64;
//         let seed = (SystemTime::now()
//             .duration_since(SystemTime::UNIX_EPOCH).unwrap()
//             .as_nanos()
//             as u64) // truncates
//             .wrapping_add(random_constant)
//             .max(1);
//         seed
//     }

//     /// Creates the RNG with a seed. All seed values aside from zero
//     /// are equally good (since this is a xorshift RNG which is
//     /// simply cycling through [1, ..., 2^64-1]). Zero is shifted
//     /// to one if passed in.
//     pub fn new(seed: u64) -> Self {
//         Self { state: seed.max(1) }
//     }

//     /// Advances the RNG and returns the next value.
//     pub fn poll(&mut self) -> u64 {
//         let mut x = self.state;
//         x ^= x << 13;
//         x ^= x >> 7;
//         x ^= x << 17;
//         self.state = x;
//         x
//     }
// }
