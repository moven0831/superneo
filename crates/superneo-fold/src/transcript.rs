//! Fiat–Shamir transcript (clean-room, Blake3-based).
//!
//! A length-framed Blake3 sponge: `absorb` folds `(label, data)` into a 32-byte
//! state; `squeeze` derives output from `(state, label, counter)` without consuming the
//! state, resetting the counter on the next absorb. All protocol randomness
//! (α, γ, sum-check challenges, the Π_RLC ρ's) is derived here, so prover and
//! verifier must absorb identical bytes in identical order (plan, Risk R3).
//!
//! Blake3 is the clean-room choice over Poseidon2; the type is the single transcript
//! used across the fold, IVC, and SNARK layers.

use superneo_commit::{Commitment, PublicParams};
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_ring::{RingElem, D};

use crate::types::CcsStructure;

/// A Blake3 Fiat–Shamir transcript.
#[derive(Clone)]
pub struct Transcript {
    state: [u8; 32],
    squeeze_ctr: u64,
}

impl Transcript {
    /// Start a transcript bound to a domain-separation label.
    pub fn new(domain: &'static [u8]) -> Self {
        let mut t = Transcript {
            state: [0u8; 32],
            squeeze_ctr: 0,
        };
        t.absorb_bytes(b"dom", domain);
        t
    }

    /// Absorb raw bytes under a label (length-framed to prevent ambiguity).
    pub fn absorb_bytes(&mut self, label: &'static [u8], data: &[u8]) {
        let mut h = blake3::Hasher::new();
        h.update(&self.state);
        h.update(&(label.len() as u64).to_le_bytes());
        h.update(label);
        h.update(&(data.len() as u64).to_le_bytes());
        h.update(data);
        self.state = *h.finalize().as_bytes();
        self.squeeze_ctr = 0;
    }

    /// Absorb a field element.
    pub fn absorb_fp(&mut self, label: &'static [u8], x: Fp) {
        self.absorb_bytes(label, &x.to_u64().to_le_bytes());
    }

    /// Absorb an extension element.
    pub fn absorb_ext(&mut self, label: &'static [u8], x: Ext2) {
        let mut buf = [0u8; 16];
        buf[..8].copy_from_slice(&x.c0.to_u64().to_le_bytes());
        buf[8..].copy_from_slice(&x.c1.to_u64().to_le_bytes());
        self.absorb_bytes(label, &buf);
    }

    /// Absorb a ring element (all `d` coefficients).
    pub fn absorb_ring(&mut self, label: &'static [u8], a: &RingElem) {
        let mut buf = [0u8; 8 * D];
        for (i, c) in a.coeffs().iter().enumerate() {
            buf[8 * i..8 * i + 8].copy_from_slice(&c.to_u64().to_le_bytes());
        }
        self.absorb_bytes(label, &buf);
    }

    /// Absorb a commitment (all `κ` ring slots).
    pub fn absorb_commitment(&mut self, label: &'static [u8], c: &Commitment) {
        for slot in &c.0 {
            self.absorb_ring(label, slot);
        }
    }

    /// A 32-byte digest of the Ajtai matrix `A` (κ × n_cols ring elements).
    fn pp_digest(pp: &PublicParams) -> [u8; 32] {
        let mut h = blake3::Hasher::new();
        h.update(b"superneo/vk/pp/v1");
        h.update(&(pp.kappa as u64).to_le_bytes());
        h.update(&(pp.n_cols as u64).to_le_bytes());
        for i in 0..pp.kappa {
            for j in 0..pp.n_cols {
                for c in pp.entry(i, j).coeffs() {
                    h.update(&c.to_u64().to_le_bytes());
                }
            }
        }
        *h.finalize().as_bytes()
    }

    /// A 32-byte digest of the CCS structure (the `t` matrices and the polynomial `f`).
    fn structure_digest(s: &CcsStructure) -> [u8; 32] {
        let mut h = blake3::Hasher::new();
        h.update(b"superneo/vk/s/v1");
        h.update(&(s.t() as u64).to_le_bytes());
        h.update(&(s.m as u64).to_le_bytes());
        for mat in &s.matrices {
            for row in mat {
                for v in row {
                    h.update(&v.to_u64().to_le_bytes());
                }
            }
        }
        h.update(&(s.f.t as u64).to_le_bytes());
        h.update(&(s.f.terms.len() as u64).to_le_bytes());
        for (coeff, exps) in &s.f.terms {
            h.update(&coeff.to_u64().to_le_bytes());
            h.update(&(exps.len() as u64).to_le_bytes());
            for &e in exps {
                h.update(&(e as u64).to_le_bytes());
            }
        }
        *h.finalize().as_bytes()
    }

    /// Bind the verifier key — the commitment key `pp` (Ajtai `A`) and the relation `s` —
    /// into the transcript. MUST be called once at the transcript origin, before any
    /// challenge is squeezed, so Fiat–Shamir binds the public parameters and the relation
    /// (otherwise an adversary free to choose `A`/`s` is not committed to them).
    pub fn absorb_vk(&mut self, pp: &PublicParams, s: &CcsStructure) {
        let pp_d = Self::pp_digest(pp);
        let s_d = Self::structure_digest(s);
        let mut h = blake3::Hasher::new();
        h.update(b"superneo/vk/v1");
        h.update(&pp_d);
        h.update(&s_d);
        self.absorb_bytes(b"vk", h.finalize().as_bytes());
    }

    /// Bind only the relation structure `s` (used inside Π_CCS so standalone folds, which
    /// build their own transcript without [`absorb_vk`], still bind the relation).
    pub fn absorb_structure(&mut self, s: &CcsStructure) {
        let s_d = Self::structure_digest(s);
        self.absorb_bytes(b"structure", &s_d);
    }

    /// Derive output from `(state, label, counter)`. The label is folded in (length-framed)
    /// so distinct challenge types are domain-separated even at the same state; the counter
    /// keeps successive squeezes under one label distinct. Does not mutate `state`.
    fn squeeze_bytes(&mut self, label: &'static [u8], n: usize) -> Vec<u8> {
        let mut h = blake3::Hasher::new();
        h.update(&self.state);
        h.update(b"squeeze");
        h.update(&(label.len() as u64).to_le_bytes());
        h.update(label);
        h.update(&self.squeeze_ctr.to_le_bytes());
        self.squeeze_ctr += 1;
        let mut reader = h.finalize_xof();
        let mut out = vec![0u8; n];
        reader.fill(&mut out);
        out
    }

    /// Squeeze a uniform field element (16 bytes → `u128` → reduce mod q).
    pub fn challenge_fp(&mut self, label: &'static [u8]) -> Fp {
        let b = self.squeeze_bytes(label, 16);
        let mut x = [0u8; 16];
        x.copy_from_slice(&b);
        Fp::reduce128(u128::from_le_bytes(x))
    }

    /// Squeeze a uniform extension element (two field elements).
    pub fn challenge_ext(&mut self, label: &'static [u8]) -> Ext2 {
        let c0 = self.challenge_fp(label);
        let c1 = self.challenge_fp(label);
        Ext2::new(c0, c1)
    }

    /// Squeeze `n` extension elements.
    pub fn challenge_ext_vec(&mut self, label: &'static [u8], n: usize) -> Vec<Ext2> {
        (0..n).map(|_| self.challenge_ext(label)).collect()
    }

    /// Squeeze `n` challenges from a strong sampling set `C`: ring elements whose
    /// coefficients are uniform in `{−bound, …, bound}` (Definition 17).
    pub fn challenge_ring_set(
        &mut self,
        label: &'static [u8],
        n: usize,
        bound: i64,
    ) -> Vec<RingElem> {
        let span = (2 * bound + 1) as u64;
        let bytes = self.squeeze_bytes(label, 8 * D * n);
        (0..n)
            .map(|i| {
                let mut c = [Fp::ZERO; D];
                for (l, slot) in c.iter_mut().enumerate() {
                    let off = 8 * (i * D + l);
                    let mut w = [0u8; 8];
                    w.copy_from_slice(&bytes[off..off + 8]);
                    let digit = (u64::from_le_bytes(w) % span) as i64 - bound;
                    *slot = Fp::from_i64(digit);
                }
                RingElem::from_coeffs(c)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn squeeze_is_label_separated() {
        // FS-3: distinct challenge labels at the same transcript state must yield distinct
        // challenges (pre-fix the label was ignored, so these were equal).
        let mut t1 = Transcript::new(b"dom");
        let mut t2 = Transcript::new(b"dom");
        assert_ne!(t1.challenge_fp(b"alpha"), t2.challenge_fp(b"beta"));
    }
}
