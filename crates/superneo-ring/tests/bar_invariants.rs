//! The bar / inner-product transform: Theorem 5 (`ct(bar(a)·b) = ⟨a,b⟩`),
//! Theorem 6 (`Mz = ct(M̄z)`), `det(G) = ±1`, and bar linearity.

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_field::fp::Fp;
use superneo_ring::bar::{bar_block, bar_ring, determinant, gram};
use superneo_ring::maps::{bar_matrix_row, embed_vector, mbar_row_dot_z};
use superneo_ring::ring::{RingElem, D};

fn rng() -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(0xB5_AD_4E_CE_DA_1C_E2_A9)
}

fn rand_block(r: &mut ChaCha8Rng) -> [Fp; D] {
    let mut c = [Fp::ZERO; D];
    for x in c.iter_mut() {
        *x = Fp::new(r.next_u64());
    }
    c
}

fn dot(a: &[Fp; D], b: &[Fp; D]) -> Fp {
    let mut acc = Fp::ZERO;
    for i in 0..D {
        acc += a[i] * b[i];
    }
    acc
}

#[test]
fn gram_determinant_is_plus_or_minus_one() {
    let det = determinant(&gram());
    assert!(
        det == Fp::ONE || det == -Fp::ONE,
        "det(G) must be ±1 (Theorem 5), got {det}"
    );
}

#[test]
fn theorem5_on_basis_vectors() {
    // ct(bar(e_i)·e_j) = ⟨e_i, e_j⟩ = δ_ij — a complete algebraic check of the transform.
    for i in 0..D {
        let mut ei = [Fp::ZERO; D];
        ei[i] = Fp::ONE;
        let bar_ei = bar_ring(&ei);
        for j in 0..D {
            let mut ej = [Fp::ZERO; D];
            ej[j] = Fp::ONE;
            let got = (bar_ei * RingElem::from_coeffs(ej)).ct();
            let want = if i == j { Fp::ONE } else { Fp::ZERO };
            assert_eq!(got, want, "Theorem 5 failed at (e{i}, e{j})");
        }
    }
}

#[test]
fn theorem5_on_random_vectors() {
    // ct( cf⁻¹(bar(a)) · cf⁻¹(b) ) = ⟨a, b⟩.
    let mut r = rng();
    for _ in 0..20_000 {
        let a = rand_block(&mut r);
        let b = rand_block(&mut r);
        let got = (bar_ring(&a) * RingElem::from_coeffs(b)).ct();
        assert_eq!(got, dot(&a, &b));
    }
}

#[test]
fn theorem6_matrix_vector_transform() {
    // Mz = ct(M̄z): the constant terms of the ring matrix-vector product recover the
    // field matrix-vector product.
    let mut r = rng();
    let m = 6usize; // rows
    let n_r = 3usize; // ring-vector length
    let n_f = D * n_r;
    for _ in 0..200 {
        let z: Vec<Fp> = (0..n_f).map(|_| Fp::new(r.next_u64())).collect();
        let z_ring = embed_vector(&z).unwrap();
        for _row in 0..m {
            let row: Vec<Fp> = (0..n_f).map(|_| Fp::new(r.next_u64())).collect();
            // direct field matrix-vector entry
            let mut mz = Fp::ZERO;
            for k in 0..n_f {
                mz += row[k] * z[k];
            }
            // ct of the ring matrix-vector entry
            let bar_row = bar_matrix_row(&row).unwrap();
            let got = mbar_row_dot_z(&bar_row, &z_ring).ct();
            assert_eq!(got, mz, "Theorem 6 failed");
        }
    }
}

#[test]
fn bar_is_linear() {
    let mut r = rng();
    for _ in 0..10_000 {
        let a = rand_block(&mut r);
        let b = rand_block(&mut r);
        let s = Fp::new(r.next_u64());

        let sum: [Fp; D] = std::array::from_fn(|i| a[i] + b[i]);
        let ba = bar_block(&a);
        let bb = bar_block(&b);
        let bsum = bar_block(&sum);
        for i in 0..D {
            assert_eq!(bsum[i], ba[i] + bb[i], "bar not additive");
        }

        let scaled: [Fp; D] = std::array::from_fn(|i| a[i] * s);
        let bscaled = bar_block(&scaled);
        for i in 0..D {
            assert_eq!(bscaled[i], ba[i] * s, "bar not scalar-linear");
        }
    }
}
