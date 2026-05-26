// Stylistic clippy lints with broad disagreement — left as-is here:
//   - needless_range_loop: many sites index into multiple parallel arrays,
//     where `for i in 0..n` is clearer than `.iter().enumerate().take(n)`.
//   - too_many_arguments: a handful of pipeline functions take >7 args by
//     design; refactoring into a struct would just shuffle complexity.
//   - type_complexity: one-off Vec<(String, usize, Vec<...>)> result types
//     that don't merit a named alias.
#![allow(clippy::needless_range_loop)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::type_complexity)]

pub mod monomer;
pub mod kmer;
pub mod io;
pub mod cmd;
