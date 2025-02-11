//! Random utility functions for Spacerocks.

use strsim::damerau_levenshtein;

pub fn find_closest_match<'a>(input: &'a str, valid_options: &'a [&'a str]) -> Option<&'a str> {
    valid_options
        .iter()
        .min_by_key(|&&a| {
            damerau_levenshtein(&input.to_uppercase(), &a)
        })
        .copied()
}