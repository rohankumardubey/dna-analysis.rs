//! Slightly modified and commented version of exercise from https://www.sotr.blog/articles/dna-analysis

use std::collections::HashMap;

/// Constant number of nucleotides in each codon
const NUCLEOTIDES_IN_CODON: usize = 3;

// Traits depend on each other. Directly needed traits:
// * Hash to serve as keys for a HashMap
// * Ord to sort the vector of codons
// * Copy to build up vector of codons using copies
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum Nucleotide {
    A,
    C,
    G,
    T,
}

/// Tuple of three nucleotides representing a codon
/// TODO: Should depend on `NUCLEOTIDES_IN_CODON`
type Codon = (Nucleotide, Nucleotide, Nucleotide);
/// Vector of codons representing a gene
type Gene = Vec<Codon>;

/// Maps a string to a vector of nucleotides
/// Panics if an invalid character occurs.
/// Valid characters are: A, C, G, T.
fn str_to_nucleotides(s: &str) -> Vec<Nucleotide> {
    s.chars()
        .map(|c| match c {
            'A' => Nucleotide::A,
            'C' => Nucleotide::C,
            'G' => Nucleotide::G,
            'T' => Nucleotide::T,
            _ => panic!("Not a nucleotide"),
        })
        .collect()
}

/// Returns a histogram of occurrences per nucleotide
fn nucleotide_frequency(dna: &[Nucleotide]) -> HashMap<Nucleotide, u32> {
    let mut frequencies: HashMap<Nucleotide, u32> = HashMap::from([
        (Nucleotide::A, 0),
        (Nucleotide::T, 0),
        (Nucleotide::G, 0),
        (Nucleotide::C, 0),
    ]);

    for nucleotide in dna {
        increment_nucleotide_count(&mut frequencies, nucleotide);
    }

    frequencies
}

fn increment_nucleotide_count(frequencies: &mut HashMap<Nucleotide, u32>, nucleotide: &Nucleotide) {
    // `unwrap` is acceptable here since the map is
    // populated with all possible nucleotide keys
    *frequencies.get_mut(nucleotide).unwrap() += 1;
}

/// Takes a string and converts it to a vector of codons.
/// Will ignore the remainder if string length is not multiple of `NUCLEOTIDES_IN_CODON`.
fn str_to_gene(s: &str) -> Gene {
    let nucleotides = str_to_nucleotides(s);
    let num_nucleotides_in_codons = nucleotides.len() - (nucleotides.len() % NUCLEOTIDES_IN_CODON);

    (0..num_nucleotides_in_codons)
        .step_by(NUCLEOTIDES_IN_CODON)
        .map(|i| (nucleotides[i], nucleotides[i + 1], nucleotides[i + 2]))
        .collect()
}

/// Check if codon is contained in a given gene
fn binary_search_for_codon(sorted_gene: &Gene, target_codon: &Codon) -> bool {
    let mut low = 0;
    let mut high = sorted_gene.len() - 1;

    while low <= high {
        let middle_codon_index = (low + high) / 2;
        let middle_codon = &sorted_gene[middle_codon_index];

        match middle_codon.cmp(target_codon) {
            std::cmp::Ordering::Less => low = middle_codon_index + 1,
            std::cmp::Ordering::Greater => high = middle_codon_index - 1,
            std::cmp::Ordering::Equal => return true,
        }
    }

    false
}

/// Naive pattern matching algorithm using brute-force.
/// For every possible pattern location, check if pattern is present.
/// Returns vector of indices of given sequence in dna.
fn naive_match(dna: &[Nucleotide], target_sequence: &[Nucleotide]) -> Vec<usize> {
    assert!(dna.len() >= target_sequence.len());

    let mut target_occurance_at = Vec::new();

    let possible_starting_positions = 0..(dna.len() - target_sequence.len() + 1);

    for i in possible_starting_positions {
        for j in 0..target_sequence.len() {
            if dna[i + j] != target_sequence[j] {
                break;
            }
            if j == target_sequence.len() - 1 {
                target_occurance_at.push(i);
            }
        }
    }

    target_occurance_at
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn nucleotide_frequency_test() {
        let dna_str = "ATATCTTAGAGGGAG";
        let freq = nucleotide_frequency(&str_to_nucleotides(dna_str));

        assert_eq!(freq.get(&Nucleotide::A), Some(&5));
        assert_eq!(freq.get(&Nucleotide::T), Some(&4));
        assert_eq!(freq.get(&Nucleotide::C), Some(&1));
        assert_eq!(freq.get(&Nucleotide::G), Some(&5));
    }

    #[test]
    fn binary_search() {
        let gene_str = "ATATCTTAGAGGGAGGGCTGAGGGTTTGAAGTCC";
        let mut gene = str_to_gene(gene_str);
        gene.sort();

        let ata: Codon = (Nucleotide::A, Nucleotide::T, Nucleotide::A);
        let atc: Codon = (Nucleotide::A, Nucleotide::T, Nucleotide::C);
        let agg: Codon = (Nucleotide::A, Nucleotide::G, Nucleotide::G);
        let tcc: Codon = (Nucleotide::T, Nucleotide::C, Nucleotide::C);

        assert!(binary_search_for_codon(&gene, &ata));
        assert!(!binary_search_for_codon(&gene, &atc));
        assert!(binary_search_for_codon(&gene, &agg));
        assert!(!binary_search_for_codon(&gene, &tcc));
    }

    #[test]
    fn naive_match_test() {
        let dna = str_to_nucleotides("ATATCTTAGAGGGAGGGAGG");
        let target_sequence = str_to_nucleotides("AGG");

        let match_indices = naive_match(&dna, &target_sequence);

        assert_eq!(match_indices.len(), 3);
        assert_eq!(match_indices[0], 9);
        assert_eq!(match_indices[1], 13);
        assert_eq!(match_indices[2], 17);
    }
}
