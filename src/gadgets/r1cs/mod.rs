use ark_ec::short_weierstrass::Projective;

pub mod r1cs;
pub mod ova;
pub mod conversion;
pub mod r1cs_var;

// the following are easier interface for using these types since we don't have to restate Projective
pub type R1CSShape<G> = r1cs::R1CSShape<Projective<G>>;
pub type R1CSInstance<G, C> = r1cs::R1CSInstance<Projective<G>, C>;
pub type R1CSWitness<G> = r1cs::R1CSWitness<Projective<G>>;
pub type RelaxedR1CSInstance<G, C> = r1cs::RelaxedR1CSInstance<Projective<G>, C>;
pub type RelaxedR1CSWitness<G> = r1cs::RelaxedR1CSWitness<Projective<G>>;
pub type OvaInstance<G, C> = ova::OvaInstance<Projective<G>, C>;
pub type OvaWitness<G> = ova::OvaWitness<Projective<G>>;
pub type RelaxedOvaInstance<G, C> = ova::RelaxedOvaInstance<Projective<G>, C>;
pub type RelaxedOvaWitness<G> = ova::RelaxedOvaWitness<Projective<G>>;

