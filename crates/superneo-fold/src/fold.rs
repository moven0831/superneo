//! The SuperNeo folding scheme: `Π_SuperNeo = Π_DEC ∘ Π_RLC ∘ Π_CCS` (Theorem 3).
//!
//! Folds `K` fresh CCS instances together with the `k` running CE instances of the
//! accumulator into a fresh accumulator of `k` CE instances of norm bound `b`.

use superneo_commit::PublicParams;

use crate::error::FoldError;
use crate::pi_ccs::{pi_ccs_prove, pi_ccs_verify};
use crate::pi_dec::{pi_dec_prove, pi_dec_verify};
use crate::pi_rlc::{pi_rlc_prove, pi_rlc_verify};
use crate::transcript::Transcript;
use crate::types::{
    CcsInstance, CcsStructure, CcsWitness, CeInstance, CeWitness, FoldProof, GlobalParams,
};

/// Prove one folding step. Returns the new accumulator (`k` CE instances), its
/// witnesses, and the fold proof.
#[allow(clippy::too_many_arguments)]
pub fn fold(
    tr: &mut Transcript,
    gp: &GlobalParams,
    pp: &PublicParams,
    s: &CcsStructure,
    fresh: &[CcsInstance],
    fresh_wit: &[CcsWitness],
    carried: &[CeInstance],
    carried_wit: &[CeWitness],
) -> Result<(Vec<CeInstance>, Vec<CeWitness>, FoldProof), FoldError> {
    let (mid_i, mid_w, ccs) = pi_ccs_prove(tr, gp, s, fresh, fresh_wit, carried, carried_wit);
    let (rlc_i, rlc_w, rlc) = pi_rlc_prove(tr, gp, &mid_i, &mid_w);
    let (dec_i, dec_w, dec) = pi_dec_prove(tr, gp, pp, s, &rlc_i, &rlc_w)?;
    Ok((dec_i, dec_w, FoldProof { ccs, rlc, dec }))
}

/// Verify one folding step. Returns the new accumulator (`k` CE instances).
pub fn verify_fold(
    tr: &mut Transcript,
    gp: &GlobalParams,
    s: &CcsStructure,
    fresh: &[CcsInstance],
    carried: &[CeInstance],
    proof: &FoldProof,
) -> Result<Vec<CeInstance>, FoldError> {
    let mid = pi_ccs_verify(tr, gp, s, fresh, carried, &proof.ccs)?;
    let rlc = pi_rlc_verify(tr, gp, &mid, &proof.rlc)?;
    let dec = pi_dec_verify(tr, gp, &rlc, &proof.dec)?;
    Ok(dec)
}
