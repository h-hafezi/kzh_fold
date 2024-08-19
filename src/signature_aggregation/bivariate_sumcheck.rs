use ark_ec::pairing::Pairing;
use transcript::IOPTranscript;

use crate::{polynomial::bivariate_poly::BivariatePolynomial, polynomial::univariate_poly::UnivariatePolynomial};



pub struct SumcheckProof<E: Pairing> {
    r_poly: UnivariatePolynomial<E::ScalarField>, // first round polynomial
    s_poly: UnivariatePolynomial<E::ScalarField>, // second round polynomial
    // more more more
}

#[allow(unused_variables)]  // XXX remove
pub fn bivariate_sumcheck<E: Pairing>(transcript: &mut IOPTranscript<E::ScalarField>, f_poly: &BivariatePolynomial<E::ScalarField>) -> SumcheckProof<E> {
    let r_poly = f_poly.sum_partial_evaluations_in_domain();

    // XXX Do more stuff..
    // XXX If you need to compute a quotient poly use arkworks' divide_by_vanishing_poly()
    //     for example see https://gist.github.com/shuklaayush/b92e6b53b0ff8571c0e73d42b504f7e3

    // Squeeze alpha challenge
    transcript.append_serializable_element(b"r_poly", &[r_poly]).unwrap();
    let alpha = transcript.get_and_append_challenge(b"alpha").unwrap();

    let s_poly = f_poly.partially_evaluate_at_x(&alpha);

    unimplemented!();
}
