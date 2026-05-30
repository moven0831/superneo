//! The foldable Reed–Solomon code underlying the BaseFold PCS (M6).
//!
//! A multilinear polynomial in `ν` variables is represented by either its evaluation
//! table over `{0,1}^ν` (the form our witnesses arrive in) or its coefficient vector
//! in the monomial basis `Π_k X_k^{i_k}`. The two are related by the subset-sum (zeta)
//! and Möbius transforms. BaseFold Reed–Solomon-encodes the **coefficient** vector over
//! a smooth multiplicative subgroup of `F_q` and folds it FRI-style; the key invariant
//! is that folding the codeword agrees with folding the coefficients
//! (`encode(fold_coeffs(a, α)) = fold_codeword(encode(a), α)`), which ties the FRI
//! folding to binding a sum-check variable.
//!
//! Goldilocks has 2-adicity 32 (`q − 1 = 2^32·(2^32 − 1)`), and `7` is a multiplicative
//! generator, so a primitive `2^b`-th root of unity is `7^{(q−1) ≫ b}`.

use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_field::Q;

/// The Goldilocks multiplicative generator.
const GENERATOR: u64 = 7;

/// The 2-adicity of `F_q^*` (`q − 1 = 2^32 · (2^32 − 1)`).
pub const TWO_ADICITY: usize = 32;

/// A primitive `2^{log_n}`-th root of unity in `F_q`.
///
/// `root(log_n) = 7^{(q−1) ≫ log_n}`, so `root(log_n)^2 = root(log_n − 1)` and the
/// evaluation domains `D_i = ⟨root(log_n − i)⟩` are nested by squaring — exactly the
/// structure FRI folding needs.
pub fn root_of_unity(log_n: usize) -> Fp {
    assert!(log_n <= TWO_ADICITY, "no 2^{log_n}-th root of unity in F_q");
    Fp::new(GENERATOR).pow((Q - 1) >> log_n)
}

/// `log2` of a power of two.
pub fn log2(n: usize) -> usize {
    assert!(n.is_power_of_two(), "{n} is not a power of two");
    n.trailing_zeros() as usize
}

// ---- Multilinear basis transforms --------------------------------------------

/// Evaluation table → monomial coefficients (Möbius transform).
///
/// `a[x] = Σ_{i ⊆ x} (−1)^{|x|−|i|} e[i]`, computed by the in-place butterfly
/// `a[hi] -= a[lo]` over each bit. Inverse of [`coeffs_to_evals`].
pub fn evals_to_coeffs(evals: &[Fp]) -> Vec<Fp> {
    let n = evals.len();
    assert!(n.is_power_of_two());
    let mut a = evals.to_vec();
    let mut step = 1;
    while step < n {
        let mut base = 0;
        while base < n {
            for j in base..base + step {
                let lo = a[j];
                a[j + step] -= lo;
            }
            base += 2 * step;
        }
        step <<= 1;
    }
    a
}

/// Monomial coefficients → evaluation table (subset-sum / zeta transform).
///
/// `e[x] = Σ_{i ⊆ x} a[i]`. Inverse of [`evals_to_coeffs`].
pub fn coeffs_to_evals(coeffs: &[Fp]) -> Vec<Fp> {
    let n = coeffs.len();
    assert!(n.is_power_of_two());
    let mut e = coeffs.to_vec();
    let mut step = 1;
    while step < n {
        let mut base = 0;
        while base < n {
            for j in base..base + step {
                let lo = e[j];
                e[j + step] += lo;
            }
            base += 2 * step;
        }
        step <<= 1;
    }
    e
}

// ---- NTT-based Reed–Solomon encoding -----------------------------------------

/// In-place radix-2 decimation-in-time NTT of length `n = a.len()` (a power of two),
/// using the supplied primitive `n`-th root of unity. Maps coefficients to
/// `[P(root^0), …, P(root^{n−1})]`.
fn ntt_in_place(a: &mut [Fp], root: Fp) {
    let n = a.len();
    if n <= 1 {
        return;
    }
    // Bit-reversal permutation.
    let log_n = log2(n);
    for i in 0..n {
        let j = (i as u32).reverse_bits() >> (32 - log_n);
        if (j as usize) > i {
            a.swap(i, j as usize);
        }
    }
    // Butterflies. The length-`len` stage uses the primitive `len`-th root.
    let mut len = 2;
    while len <= n {
        // w_len = root^{n/len} is a primitive len-th root of unity.
        let w_len = root.pow((n / len) as u64);
        let half = len / 2;
        let mut base = 0;
        while base < n {
            let mut w = Fp::ONE;
            for k in 0..half {
                let u = a[base + k];
                let v = a[base + k + half] * w;
                a[base + k] = u + v;
                a[base + k + half] = u - v;
                w *= w_len;
            }
            base += len;
        }
        len <<= 1;
    }
}

/// Reed–Solomon encode a coefficient vector: evaluate the degree-`<N` polynomial with
/// coefficients `coeffs` (length `N`, a power of two) on the size-`N·2^{log_blowup}`
/// subgroup `⟨root(log N + log_blowup)⟩`.
pub fn encode(coeffs: &[Fp], log_blowup: usize) -> Vec<Fp> {
    let n = coeffs.len();
    assert!(n.is_power_of_two());
    let n0 = n << log_blowup;
    let mut buf = vec![Fp::ZERO; n0];
    buf[..n].copy_from_slice(coeffs);
    ntt_in_place(&mut buf, root_of_unity(log2(n0)));
    buf
}

// ---- Folding -----------------------------------------------------------------

/// Fold a coefficient vector by binding the least-significant variable to `α`:
/// `a'[j] = a[2j] + α·a[2j+1]`.
pub fn fold_coeffs(coeffs: &[Ext2], alpha: Ext2) -> Vec<Ext2> {
    coeffs
        .chunks_exact(2)
        .map(|c| c[0] + alpha * c[1])
        .collect()
}

/// Fold a codeword over the size-`n` domain `⟨root(log n)⟩` to the size-`n/2` domain
/// `⟨root(log n)²⟩`, FRI-style:
/// `C'[t] = (C[t] + C[t+n/2])/2 + α·(C[t] − C[t+n/2])/(2·root(log n)^t)`.
///
/// This is the evaluation-domain image of [`fold_coeffs`] (see the module-level
/// commutation invariant).
pub fn fold_codeword(codeword: &[Ext2], alpha: Ext2) -> Vec<Ext2> {
    let n = codeword.len();
    assert!(n.is_power_of_two() && n >= 2);
    let half = n / 2;
    let g = root_of_unity(log2(n));
    let two_inv = Ext2::from_base(Fp::new(2).inv().expect("2 ≠ 0 in F_q"));
    let mut out = Vec::with_capacity(half);
    let mut s = Fp::ONE; // root^t
    for t in 0..half {
        let lo = codeword[t];
        let hi = codeword[t + half];
        let s_inv = Ext2::from_base(s.inv().expect("domain points are nonzero"));
        let even = (lo + hi) * two_inv;
        let odd = (lo - hi) * two_inv * s_inv;
        out.push(even + alpha * odd);
        s *= g;
    }
    out
}

/// Lift a base-field codeword to the extension field.
pub fn lift(codeword: &[Fp]) -> Vec<Ext2> {
    codeword.iter().map(|&x| Ext2::from_base(x)).collect()
}
