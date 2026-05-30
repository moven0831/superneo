//! **Phase 2: the recursive verifier circuit (M7).**
//!
//! Expresses the folding verifier's work as CCS constraints, so that a satisfying
//! assignment certifies that the verifier accepts — the step that closes a true
//! recursive IVC loop. This PoC implements the dominant component, the **sum-check
//! verifier** ([`sumcheck_verifier`]), built from a reusable R1CS constraint-system
//! builder ([`cs`]) and extension-field gadgets ([`gadgets`], Karatsuba `K`-mul). A
//! synthesized system [`finalize`](cs::ConstraintSystem::finalize)s into the same
//! `t = 4` CCS shape `fold` consumes, demonstrating loop closure at the relation level.
//!
//! Scope (PoC, documented residuals): the Fiat–Shamir challenges are supplied as advice
//! (Blake3-in-circuit is out of range); the Π_CCS final-`Q(r′)` reconstruction and the
//! Π_RLC / Π_DEC verifier checks are not yet synthesized; and folding the resulting CCS
//! instance through the *norm-bounded* scheme needs the augmented-witness decomposition.

pub mod cs;
pub mod gadgets;
pub mod sumcheck_verifier;

pub use cs::{ConstraintSystem, Lc, Var};
pub use gadgets::KVar;
pub use sumcheck_verifier::synthesize_sumcheck_verifier;

use superneo_field::ext2::Ext2;

/// Build a complete recursive sum-check-verifier circuit: synthesize the verifier chain
/// and bind its final reduced claim to `expected_final` (the value a full verifier would
/// match against `Q(r′)`). Returns the finished constraint system.
pub fn build_sumcheck_verifier_circuit(
    init_claim: Ext2,
    round_evals: &[Vec<Ext2>],
    challenges: &[Ext2],
    expected_final: Ext2,
) -> ConstraintSystem {
    let mut cs = ConstraintSystem::new();
    let claim = synthesize_sumcheck_verifier(&mut cs, init_claim, round_evals, challenges);
    let expected = KVar::alloc(&mut cs, expected_final);
    claim.assert_eq(&mut cs, &expected);
    cs
}
