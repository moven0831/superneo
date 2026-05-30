//! A minimal R1CS constraint-system builder that finalizes into a fold-compatible
//! [`CcsStructure`] (M7).
//!
//! Variables index a witness vector (`z[0] = 1` is the constant "one" wire). Constraints
//! are rank-1: `⟨a, z⟩ · ⟨b, z⟩ = ⟨c, z⟩`, with `a, b, c` linear combinations. The
//! builder both records the constraints and evaluates them against the live witness, so
//! [`ConstraintSystem::is_satisfied`] is the in-memory satisfaction oracle the recursive
//! tests use. [`ConstraintSystem::finalize`] emits the `t = 4` SuperNeo CCS encoding
//! `(M_1=I, M_2=A, M_3=B, M_4=C)`, `f = X_2·X_3 − X_4` (Remark 3), padded to a square
//! `m × m` with `m` a multiple of the ring degree `d` — exactly the shape `fold` consumes.

use superneo_field::fp::Fp;
use superneo_fold::types::field_matvec;
use superneo_fold::{CcsStructure, SparsePoly};
use superneo_ring::D;

/// A witness-vector index.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Var(pub usize);

/// A linear combination `Σ coeff_i · z[var_i]` over the witness.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Lc(pub Vec<(usize, Fp)>);

impl Lc {
    /// The zero combination.
    pub fn zero() -> Lc {
        Lc(Vec::new())
    }

    /// A single variable with coefficient one.
    pub fn from_var(v: Var) -> Lc {
        Lc(vec![(v.0, Fp::ONE)])
    }

    /// `self + other`.
    pub fn add(&self, other: &Lc) -> Lc {
        let mut out = self.0.clone();
        out.extend_from_slice(&other.0);
        Lc(out)
    }

    /// `self − other`.
    pub fn sub(&self, other: &Lc) -> Lc {
        let mut out = self.0.clone();
        out.extend(other.0.iter().map(|&(v, c)| (v, -c)));
        Lc(out)
    }

    /// `s · self`.
    pub fn scale(&self, s: Fp) -> Lc {
        Lc(self.0.iter().map(|&(v, c)| (v, c * s)).collect())
    }
}

/// An R1CS constraint system that finalizes into a SuperNeo CCS instance.
#[derive(Clone, Debug)]
pub struct ConstraintSystem {
    witness: Vec<Fp>,
    a: Vec<Lc>,
    b: Vec<Lc>,
    c: Vec<Lc>,
}

impl Default for ConstraintSystem {
    fn default() -> Self {
        Self::new()
    }
}

impl ConstraintSystem {
    /// A fresh system whose only variable is the constant-one wire `z[0] = 1`.
    pub fn new() -> Self {
        ConstraintSystem {
            witness: vec![Fp::ONE],
            a: Vec::new(),
            b: Vec::new(),
            c: Vec::new(),
        }
    }

    /// The constant-one wire.
    pub fn one(&self) -> Var {
        Var(0)
    }

    /// A constant as a linear combination (`value · one`).
    pub fn constant(&self, value: Fp) -> Lc {
        Lc(vec![(0, value)])
    }

    /// Allocate a witness variable holding `value`.
    pub fn alloc(&mut self, value: Fp) -> Var {
        let idx = self.witness.len();
        self.witness.push(value);
        Var(idx)
    }

    /// Evaluate a linear combination against the live witness.
    pub fn eval(&self, lc: &Lc) -> Fp {
        let mut acc = Fp::ZERO;
        for &(v, c) in &lc.0 {
            acc += c * self.witness[v];
        }
        acc
    }

    /// Impose `⟨a, z⟩ · ⟨b, z⟩ = ⟨c, z⟩`.
    ///
    /// The constraint is *recorded*, not checked here — synthesizing from tampered
    /// advice yields a system that [`is_satisfied`](Self::is_satisfied) reports as
    /// unsatisfied (rather than panicking), which is what the recursive tests rely on.
    pub fn enforce(&mut self, a: Lc, b: Lc, c: Lc) {
        self.a.push(a);
        self.b.push(b);
        self.c.push(c);
    }

    /// Allocate `p = ⟨a, z⟩ · ⟨b, z⟩` and return it.
    pub fn mul(&mut self, a: &Lc, b: &Lc) -> Var {
        let value = self.eval(a) * self.eval(b);
        let p = self.alloc(value);
        self.enforce(a.clone(), b.clone(), Lc::from_var(p));
        p
    }

    /// Impose `⟨a, z⟩ = ⟨b, z⟩` (as `(a − b)·1 = 0`).
    pub fn assert_eq(&mut self, a: &Lc, b: &Lc) {
        self.enforce(a.sub(b), Lc::from_var(self.one()), Lc::zero());
    }

    /// Number of variables (witness length).
    pub fn num_vars(&self) -> usize {
        self.witness.len()
    }

    /// Number of rank-1 constraints.
    pub fn num_constraints(&self) -> usize {
        self.a.len()
    }

    /// Whether every constraint holds under the live witness.
    pub fn is_satisfied(&self) -> bool {
        self.a
            .iter()
            .zip(&self.b)
            .zip(&self.c)
            .all(|((a, b), c)| self.eval(a) * self.eval(b) == self.eval(c))
    }

    /// Densify a linear-combination list into a length-`m` row.
    fn dense_row(lc: &Lc, m: usize) -> Vec<Fp> {
        let mut row = vec![Fp::ZERO; m];
        for &(v, c) in &lc.0 {
            row[v] += c;
        }
        row
    }

    /// Finalize into a SuperNeo CCS structure and its padded witness vector.
    ///
    /// Returns `(structure, z)` where `z` (length `m`) satisfies the CCS relation iff the
    /// system [`is_satisfied`](Self::is_satisfied). `m` is the smallest multiple of `d`
    /// at least `max(num_vars, num_constraints)`.
    pub fn finalize(&self) -> (CcsStructure, Vec<Fp>) {
        let need = self.num_vars().max(self.num_constraints()).max(1);
        let m = need.div_ceil(D) * D;

        let mut id = vec![vec![Fp::ZERO; m]; m];
        for (i, row) in id.iter_mut().enumerate() {
            row[i] = Fp::ONE;
        }
        let mut mat_a = vec![vec![Fp::ZERO; m]; m];
        let mut mat_b = vec![vec![Fp::ZERO; m]; m];
        let mut mat_c = vec![vec![Fp::ZERO; m]; m];
        for (r, ((a, b), c)) in self.a.iter().zip(&self.b).zip(&self.c).enumerate() {
            mat_a[r] = Self::dense_row(a, m);
            mat_b[r] = Self::dense_row(b, m);
            mat_c[r] = Self::dense_row(c, m);
        }

        let f = SparsePoly::new(
            4,
            vec![
                (Fp::ONE, vec![0, 1, 1, 0]),  // +X_2·X_3
                (-Fp::ONE, vec![0, 0, 0, 1]), // −X_4
            ],
        );
        let structure = CcsStructure::new(vec![id, mat_a, mat_b, mat_c], f);

        let mut z = vec![Fp::ZERO; m];
        z[..self.witness.len()].copy_from_slice(&self.witness);
        (structure, z)
    }

    /// Loop-closure check: the finalized CCS instance's witness satisfies `A z ∘ B z = C z`
    /// (the CCS satisfaction predicate the folding scheme's `F`-path enforces).
    pub fn finalized_relation_holds(&self) -> bool {
        let (s, z) = self.finalize();
        let az = field_matvec(&s.matrices[1], &z);
        let bz = field_matvec(&s.matrices[2], &z);
        let cz = field_matvec(&s.matrices[3], &z);
        az.iter().zip(&bz).zip(&cz).all(|((&a, &b), &c)| a * b == c)
    }
}
