//! Ajtai public parameters: the random matrix `A ∈ R_F^{κ×n_cols}` and the
//! commitment map `Commit(A, z) = A·z` (Definitions 4 & 18).
//!
//! The matrix is expanded deterministically from a 32-byte seed via ChaCha8, so
//! `setup_seeded` is reproducible (the basis for differential and determinism tests).

use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use superneo_field::fp::Fp;
use superneo_field::Q;
use superneo_ring::{RingElem, D};

use crate::commit::Commitment;
use crate::error::CommitError;

/// A uniform field element in `[0, q)` by rejection sampling.
fn uniform_fp<R: RngCore>(rng: &mut R) -> Fp {
    loop {
        let x = rng.next_u64();
        if x < Q {
            return Fp::new(x);
        }
    }
}

/// A uniform ring element (each coefficient uniform in `F`).
fn uniform_ring<R: RngCore>(rng: &mut R) -> RingElem {
    let mut c = [Fp::ZERO; D];
    for ci in c.iter_mut() {
        *ci = uniform_fp(rng);
    }
    RingElem::from_coeffs(c)
}

/// Ajtai public parameters: a `κ × n_cols` matrix over `R_F`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PublicParams {
    /// Module rank `κ`.
    pub kappa: usize,
    /// Number of columns (the ring-vector length `n_R` of committed witnesses).
    pub n_cols: usize,
    /// Row-major matrix `A`, `a[i][j] ∈ R_F`.
    a: Vec<Vec<RingElem>>,
}

impl PublicParams {
    /// Expand `A` deterministically from `seed` (Setup, Definition 18).
    pub fn setup_seeded(seed: [u8; 32], kappa: usize, n_cols: usize) -> Self {
        let mut rng = ChaCha8Rng::from_seed(seed);
        let a = (0..kappa)
            .map(|_| (0..n_cols).map(|_| uniform_ring(&mut rng)).collect())
            .collect();
        PublicParams { kappa, n_cols, a }
    }

    /// Borrow the matrix entry `a[i][j]`.
    pub fn entry(&self, i: usize, j: usize) -> &RingElem {
        &self.a[i][j]
    }

    /// Commit to a ring vector: `Commit(A, z) = A·z ∈ R_F^κ`.
    ///
    /// Errors with [`CommitError::DimMismatch`] if `z.len() != n_cols`.
    pub fn commit(&self, z: &[RingElem]) -> Result<Commitment, CommitError> {
        if z.len() != self.n_cols {
            return Err(CommitError::DimMismatch(format!(
                "witness length {} != n_cols {}",
                z.len(),
                self.n_cols
            )));
        }
        let mut out = Vec::with_capacity(self.kappa);
        for row in &self.a {
            let mut acc = RingElem::ZERO;
            for (aij, zj) in row.iter().zip(z.iter()) {
                acc = acc + (*aij * *zj);
            }
            out.push(acc);
        }
        Ok(Commitment(out))
    }
}
