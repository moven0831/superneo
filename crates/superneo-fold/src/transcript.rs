//! Fiat–Shamir transcript (clean-room, Blake3-based).
//!
//! A length-framed Blake3 sponge: `absorb` folds `(label, data)` into a 32-byte
//! state; `squeeze` derives output from `(state, counter)` without consuming the
//! state, resetting the counter on the next absorb. All protocol randomness
//! (α, γ, sum-check challenges, the Π_RLC ρ's) is derived here, so prover and
//! verifier must absorb identical bytes in identical order (plan, Risk R3).
//!
//! Blake3 is the clean-room choice over Poseidon2; the type is the single transcript
//! used across the fold, IVC, and SNARK layers.

use superneo_commit::Commitment;
use superneo_field::ext2::Ext2;
use superneo_field::fp::Fp;
use superneo_ring::{RingElem, D};

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

    fn squeeze_bytes(&mut self, n: usize) -> Vec<u8> {
        let mut h = blake3::Hasher::new();
        h.update(&self.state);
        h.update(b"squeeze");
        h.update(&self.squeeze_ctr.to_le_bytes());
        self.squeeze_ctr += 1;
        let mut reader = h.finalize_xof();
        let mut out = vec![0u8; n];
        reader.fill(&mut out);
        out
    }

    /// Squeeze a uniform field element (16 bytes → `u128` → reduce mod q).
    pub fn challenge_fp(&mut self, _label: &'static [u8]) -> Fp {
        let b = self.squeeze_bytes(16);
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
        _label: &'static [u8],
        n: usize,
        bound: i64,
    ) -> Vec<RingElem> {
        let span = (2 * bound + 1) as u64;
        let bytes = self.squeeze_bytes(8 * D * n);
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
