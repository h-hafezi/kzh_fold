use ark_ec::CurveGroup;
use ark_ec::VariableBaseMSM;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use rand::thread_rng;

#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct MultiCommitGens<G>
where
    G: CurveGroup,
{
    pub n: usize,
    pub G: Vec<G>,
    pub h: G,
}

impl<G: CurveGroup> MultiCommitGens<G> {
    pub fn new(n: usize, _label: &[u8]) -> Self {
        let mut gens: Vec<G> = Vec::new();
        for _ in 0..n + 1 {
            gens.push(G::rand(&mut thread_rng()));
        }

        MultiCommitGens {
            n,
            G: gens[..n].to_vec(),
            h: gens[n],
        }
    }

    pub fn clone(&self) -> Self {
        MultiCommitGens {
            n: self.n,
            h: self.h,
            G: self.G.clone(),
        }
    }

    pub fn split_at(&self, mid: usize) -> (Self, Self) {
        let (G1, G2) = self.G.split_at(mid);

        (
            MultiCommitGens {
                n: G1.len(),
                G: G1.to_vec(),
                h: self.h,
            },
            MultiCommitGens {
                n: G2.len(),
                G: G2.to_vec(),
                h: self.h,
            },
        )
    }
}

pub trait Commitments<G: CurveGroup>: Sized {
    fn commit(&self, blind: &G::ScalarField, gens_n: &MultiCommitGens<G>) -> G;
    fn batch_commit(inputs: &[Self], blind: &G::ScalarField, gens_n: &MultiCommitGens<G>) -> G;
}

impl<G: CurveGroup> Commitments<G> for G::ScalarField {
    fn commit(&self, blind: &G::ScalarField, gens_n: &MultiCommitGens<G>) -> G {
        assert_eq!(gens_n.n, 1);

        gens_n.G[0] * self + gens_n.h * blind
    }

    fn batch_commit(inputs: &[Self], blind: &G::ScalarField, gens_n: &MultiCommitGens<G>) -> G {
        assert_eq!(gens_n.n, inputs.len());

        let mut bases = CurveGroup::normalize_batch(gens_n.G.as_ref());
        let mut scalars = inputs.to_vec();
        bases.push(gens_n.h.into_affine());
        scalars.push(*blind);

        VariableBaseMSM::msm(bases.as_ref(), scalars.as_ref()).unwrap()
    }
}