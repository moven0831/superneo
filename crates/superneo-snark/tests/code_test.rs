//! Foldable-code invariants: roots of unity have the right order, the NTT agrees with
//! naive evaluation, the multilinear basis transforms invert, and — the load-bearing
//! invariant for FRI — folding the codeword agrees with folding the coefficients.

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_snark::code::{
    coeffs_to_evals, encode, evals_to_coeffs, fold_codeword, fold_coeffs, lift, log2, root_of_unity,
};

fn rng() -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(0x5151_2626_3737_4848)
}

fn rand_fp(r: &mut ChaCha8Rng) -> Fp {
    Fp::new(r.next_u64())
}

#[test]
fn roots_of_unity_have_exact_order() {
    for bits in 1..=20 {
        let w = root_of_unity(bits);
        // w^{2^bits} = 1 and w^{2^{bits-1}} = −1 (primitive order exactly 2^bits).
        assert_eq!(w.pow(1u64 << bits), Fp::ONE, "order divides 2^{bits}");
        assert_eq!(
            w.pow(1u64 << (bits - 1)),
            -Fp::ONE,
            "order is exactly 2^{bits}"
        );
    }
    // root(b)^2 = root(b-1): the domains nest by squaring.
    for bits in 2..=20 {
        assert_eq!(root_of_unity(bits).square(), root_of_unity(bits - 1));
    }
}

/// Naive evaluation of a coefficient polynomial on the size-`n` subgroup.
fn naive_encode(coeffs: &[Fp], n: usize) -> Vec<Fp> {
    let w = root_of_unity(log2(n));
    (0..n)
        .map(|t| {
            let x = w.pow(t as u64);
            let mut acc = Fp::ZERO;
            let mut xp = Fp::ONE;
            for &c in coeffs {
                acc += c * xp;
                xp *= x;
            }
            acc
        })
        .collect()
}

#[test]
fn ntt_matches_naive_evaluation() {
    let mut r = rng();
    for log_n in 1..=8 {
        let n = 1 << log_n;
        let coeffs: Vec<Fp> = (0..n).map(|_| rand_fp(&mut r)).collect();
        // encode with blow-up 0 = plain NTT on the size-n domain.
        let got = encode(&coeffs, 0);
        let want = naive_encode(&coeffs, n);
        assert_eq!(got, want, "NTT disagrees with naive eval at log_n={log_n}");
    }
}

#[test]
fn basis_transforms_invert() {
    let mut r = rng();
    for log_n in 0..=10 {
        let n = 1 << log_n;
        let evals: Vec<Fp> = (0..n).map(|_| rand_fp(&mut r)).collect();
        let coeffs = evals_to_coeffs(&evals);
        assert_eq!(coeffs_to_evals(&coeffs), evals, "möbius/zeta not inverse");
    }
}

#[test]
fn fold_codeword_matches_fold_coeffs() {
    // encode(fold_coeffs(a, α)) == fold_codeword(encode(a), α). Use α ∈ F so the
    // comparison stays in the base field (the K case is exercised in the PCS tests).
    let mut r = rng();
    for log_n in 1..=8 {
        let n = 1 << log_n;
        let log_blowup = 2;
        let coeffs: Vec<Fp> = (0..n).map(|_| rand_fp(&mut r)).collect();
        let alpha = rand_fp(&mut r);

        let folded_coeffs = fold_coeffs(&lift(&coeffs), Ext2::from_base(alpha));
        // Re-encode the folded coefficients on the (halved) domain.
        let folded_coeffs_fp: Vec<Fp> = folded_coeffs.iter().map(|c| c.c0).collect();
        let lhs = encode(&folded_coeffs_fp, log_blowup);

        let codeword = lift(&encode(&coeffs, log_blowup));
        let rhs = fold_codeword(&codeword, Ext2::from_base(alpha));
        let rhs_fp: Vec<Fp> = rhs.iter().map(|c| c.c0).collect();

        assert_eq!(lhs, rhs_fp, "fold commutation fails at log_n={log_n}");
        assert!(rhs.iter().all(|c| c.c1 == Fp::ZERO), "α∈F ⇒ result in F");
    }
}
