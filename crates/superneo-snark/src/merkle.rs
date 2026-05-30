//! A binary Merkle tree over extension-field codeword leaves (Blake3), used to
//! commit each FRI layer of the BaseFold PCS.
//!
//! Leaves are the codeword entries (`Ext2`), serialized to 16 bytes and domain-tagged.
//! `commit` returns the root; `open` returns a leaf with its authentication path;
//! `verify_path` recomputes the root from a leaf + path. Corrupting any leaf or path
//! node changes the recomputed root, which is how the FRI query phase detects tampering.

use superneo_field::ext2::Ext2;

/// Domain tags keep leaf and internal-node hashes in disjoint ranges.
const LEAF_TAG: &[u8] = b"superneo/merkle/leaf";
const NODE_TAG: &[u8] = b"superneo/merkle/node";

/// Hash one codeword entry into a 32-byte leaf digest.
fn hash_leaf(x: Ext2) -> [u8; 32] {
    let mut h = blake3::Hasher::new();
    h.update(LEAF_TAG);
    h.update(&x.c0.to_u64().to_le_bytes());
    h.update(&x.c1.to_u64().to_le_bytes());
    *h.finalize().as_bytes()
}

/// Hash two child digests into their parent.
fn hash_node(l: &[u8; 32], r: &[u8; 32]) -> [u8; 32] {
    let mut h = blake3::Hasher::new();
    h.update(NODE_TAG);
    h.update(l);
    h.update(r);
    *h.finalize().as_bytes()
}

/// A Merkle tree over a power-of-two number of leaves; stores every level so paths
/// are cheap to produce.
#[derive(Clone, Debug)]
pub struct MerkleTree {
    /// `levels[0]` are the leaf digests; `levels.last()` is the single root.
    levels: Vec<Vec<[u8; 32]>>,
}

/// An authentication path: sibling digests from the leaf up to (but excluding) the root.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MerklePath {
    /// Sibling digests, leaf level first.
    pub siblings: Vec<[u8; 32]>,
}

impl MerkleTree {
    /// Build the tree over a codeword (length a power of two).
    pub fn commit(codeword: &[Ext2]) -> MerkleTree {
        let n = codeword.len();
        assert!(
            n.is_power_of_two() && n >= 1,
            "leaf count must be a power of two"
        );
        let leaves: Vec<[u8; 32]> = codeword.iter().map(|&x| hash_leaf(x)).collect();
        let mut levels = vec![leaves];
        while levels.last().unwrap().len() > 1 {
            let cur = levels.last().unwrap();
            let next: Vec<[u8; 32]> = cur
                .chunks_exact(2)
                .map(|pair| hash_node(&pair[0], &pair[1]))
                .collect();
            levels.push(next);
        }
        MerkleTree { levels }
    }

    /// The Merkle root.
    pub fn root(&self) -> [u8; 32] {
        self.levels.last().unwrap()[0]
    }

    /// The authentication path for leaf `index`.
    pub fn open(&self, index: usize) -> MerklePath {
        let mut idx = index;
        let mut siblings = Vec::with_capacity(self.levels.len() - 1);
        for level in &self.levels[..self.levels.len() - 1] {
            siblings.push(level[idx ^ 1]);
            idx >>= 1;
        }
        MerklePath { siblings }
    }
}

/// Recompute the root implied by `leaf` at position `index` with authentication
/// `path`, for a tree of `n` leaves. Returns `true` iff it matches `root`.
pub fn verify_path(root: &[u8; 32], n: usize, index: usize, leaf: Ext2, path: &MerklePath) -> bool {
    if !n.is_power_of_two() || index >= n {
        return false;
    }
    // `n` is a power of two (checked above), so its tree has `log2(n)` levels.
    if path.siblings.len() != n.trailing_zeros() as usize {
        return false;
    }
    let mut acc = hash_leaf(leaf);
    let mut idx = index;
    for sib in &path.siblings {
        acc = if idx & 1 == 0 {
            hash_node(&acc, sib)
        } else {
            hash_node(sib, &acc)
        };
        idx >>= 1;
    }
    &acc == root
}
