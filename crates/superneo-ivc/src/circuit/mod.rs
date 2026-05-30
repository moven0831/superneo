//! **Phase 2: the recursive verifier circuit (M7).**
//!
//! Expresses the folding verifier's work as CCS constraints, so that a satisfying
//! assignment certifies that the verifier accepts — the step that closes a true recursive
//! IVC loop. This PoC synthesizes the **complete `verify_fold` chain** over advice
//! challenges — `Π_CCS` then `Π_RLC` then `Π_DEC`:
//!
//! * the **sum-check verifier** ([`sumcheck_verifier`]) and the **Π_CCS `Q(r')`
//!   reconstruction** ([`pi_ccs_verifier`]) tying the reduced claim to the claimed evals;
//! * the **Π_RLC ring linear combination** ([`pi_rlc_verifier`], on the [`ring_gadgets`]
//!   `R_F`/`R_K` multiplication gadgets);
//! * the **Π_DEC decomposition additivity** check ([`pi_dec_verifier`]).
//!
//! All are built from a reusable R1CS constraint-system builder ([`cs`]) and the
//! extension-field gadget ([`gadgets`], Karatsuba `K`-mul). A synthesized system
//! [`finalize`](cs::ConstraintSystem::finalize)s into the same `t = 4` CCS shape `fold`
//! consumes, demonstrating loop closure at the relation level.
//!
//! Scope (PoC, documented residuals): the Fiat–Shamir challenges are supplied as advice
//! (Blake3-in-circuit is out of range — a later phase introduces an arithmetization-
//! friendly transcript so the challenges are derived in-circuit); and folding the
//! resulting CCS instance back through the *norm-bounded* scheme needs the augmented-
//! witness decomposition. The ring-combination cost (`D²` base mults per product) means
//! the Π_RLC step runs at scaled-down `κ, K, n_R` here.

pub mod cs;
pub mod gadgets;
pub mod pi_ccs_verifier;
pub mod pi_dec_verifier;
pub mod pi_rlc_verifier;
pub mod ring_gadgets;
pub mod sumcheck_verifier;

pub use cs::{ConstraintSystem, Lc, Var};
pub use gadgets::KVar;
pub use pi_ccs_verifier::build_pi_ccs_verifier_circuit;
pub use pi_dec_verifier::build_pi_dec_verifier_circuit;
pub use pi_rlc_verifier::build_pi_rlc_verifier_circuit;
pub use ring_gadgets::{RFVar, RKVar};
pub use sumcheck_verifier::synthesize_sumcheck_verifier;

use superneo_field::ext2::Ext2;

/// Build a standalone recursive sum-check-verifier circuit: synthesize the verifier chain
/// and bind its final reduced claim to `expected_final` (the value a full verifier would
/// match against `Q(r')`). Both endpoints `init_claim` and `expected_final` are public
/// values of the outer relation, so they are pinned as constant wires — the circuit is
/// satisfiable iff the (advice) transcript reduces the public `init_claim` to the public
/// `expected_final`. Returns the finished constraint system.
///
/// For a complete Π_CCS verifier that reconstructs `Q(r')` in-circuit (rather than taking
/// `expected_final` as advice), see [`build_pi_ccs_verifier_circuit`].
pub fn build_sumcheck_verifier_circuit(
    init_claim: Ext2,
    round_evals: &[Vec<Ext2>],
    challenges: &[Ext2],
    expected_final: Ext2,
) -> ConstraintSystem {
    let mut cs = ConstraintSystem::new();
    let (claim, _r_wires) =
        synthesize_sumcheck_verifier(&mut cs, init_claim, round_evals, challenges);
    let expected = KVar::constant(&cs, expected_final);
    claim.assert_eq(&mut cs, &expected);
    cs
}
