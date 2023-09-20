#![allow(non_snake_case)]

mod pedersen;
mod r1cs;
pub mod traits;
mod transcript;

use crate::pedersen::Pedersen;
use crate::traits::CurveLike;
use halo2curves::group::ff::Field;
use r1cs::R1CS;
use rand;
use std::marker::PhantomData;
use transcript::PoseidonTranscript;

#[derive(Clone, Debug)]
pub struct NARKProof<C: CurveLike> {
    c_A: C::CurveExt,
    c_B: C::CurveExt,
    c_C: C::CurveExt,
    c_A_prime: C::CurveExt,
    c_B_prime: C::CurveExt,
    c_C_prime: C::CurveExt,
    c_1: C::CurveExt,
    c_2: C::CurveExt,
    s: Vec<C::ScalarExt>,
    sigma_A: C::ScalarExt,
    sigma_B: C::ScalarExt,
    sigma_C: C::ScalarExt,
    sigma_o: C::ScalarExt,
}

pub struct NARK<C: CurveLike> {
    _marker: PhantomData<C>,
}

// Implements the Sigma protocol in Figure 3 of https://eprint.iacr.org/2020/1618.pdf
// We use the notations from the paper for the variables.
impl<C: CurveLike> NARK<C> {
    pub fn new() -> Self {
        Self {
            _marker: PhantomData,
        }
    }

    fn hadamard(a: &[C::ScalarExt], b: &[C::ScalarExt]) -> Vec<C::ScalarExt> {
        assert_eq!(a.len(), b.len());

        (0..a.len())
            .map(|i| a[i] * b[i])
            .collect::<Vec<C::ScalarExt>>()
    }

    fn add_vecs(a: &[C::ScalarExt], b: &[C::ScalarExt]) -> Vec<C::ScalarExt> {
        assert_eq!(a.len(), b.len());

        (0..a.len())
            .map(|i| a[i] + b[i])
            .collect::<Vec<C::ScalarExt>>()
    }

    fn scale_vec(a: &[C::ScalarExt], b: &C::ScalarExt) -> Vec<C::ScalarExt> {
        (0..a.len())
            .map(|i| a[i] * b)
            .collect::<Vec<C::ScalarExt>>()
    }

    pub fn prove(
        r1cs: &R1CS<C::ScalarExt>,
        witness: &[C::ScalarExt],
        pub_input: &[C::ScalarExt],
        transcript: &mut PoseidonTranscript<C>,
    ) -> NARKProof<C> {
        let n = r1cs.z_len();

        let z = R1CS::construct_z(witness, pub_input);

        // Compute the matrix-vector product of all three matrices and the vector Z
        let z_A = r1cs.A.mul_vector(&z);
        let z_B = r1cs.B.mul_vector(&z);
        let z_C = r1cs.C.mul_vector(&z);

        let committer = Pedersen::<C>::new(n);

        let mut rng = rand::thread_rng();

        // Sample blinders
        let w_A = C::ScalarExt::random(&mut rng);
        let w_B = C::ScalarExt::random(&mut rng);
        let w_C = C::ScalarExt::random(&mut rng);

        // Commit to the vectors z_A, z_B, z_C
        let c_A = committer.commit(&z_A, w_A);
        let c_B = committer.commit(&z_B, w_B);
        let c_C = committer.commit(&z_C, w_C);

        // Append c_A, c_B, c_C to the transcript
        transcript.append_point(c_A);
        transcript.append_point(c_B);
        transcript.append_point(c_C);

        // Sample blinders
        let w_A_prime = C::ScalarExt::random(&mut rng);
        let w_B_prime = C::ScalarExt::random(&mut rng);
        let w_C_prime = C::ScalarExt::random(&mut rng);

        // Sample the blinders to blind the witness
        let r = (0..r1cs.num_vars)
            .map(|_| C::ScalarExt::random(&mut rng))
            .collect::<Vec<C::ScalarExt>>();

        // Add zeros for the public part of Z
        let mut r_padded = vec![C::ScalarExt::zero(); n - r1cs.num_vars];
        r_padded.extend_from_slice(&r);

        // Compute the matrix-vector product of all three matrices and the blinding vector r
        let r_A = r1cs.A.mul_vector(&r_padded);
        let r_B = r1cs.B.mul_vector(&r_padded);
        let r_C = r1cs.C.mul_vector(&r_padded);

        // Commit to the vectors r_A, r_B, r_C
        let c_A_prime = committer.commit(&r_A, w_A_prime);
        let c_B_prime = committer.commit(&r_B, w_B_prime);
        let c_C_prime = committer.commit(&r_C, w_C_prime);

        // Append c_A_prime, c_B_prime, c_C_prime to the transcript
        transcript.append_point(c_A_prime);
        transcript.append_point(c_B_prime);
        transcript.append_point(c_C_prime);

        // Sample blinders
        let w_1 = C::ScalarExt::random(&mut rng);
        let w_2 = C::ScalarExt::random(&mut rng);

        // Compute the cross-term
        let cross_1 = Self::add_vecs(&Self::hadamard(&z_A, &r_B), &Self::hadamard(&z_B, &r_A));
        let cross_2 = Self::hadamard(&r_A, &r_B);

        // Commit to the cross-terms
        let c_1 = committer.commit(&cross_1, w_1);
        let c_2 = committer.commit(&cross_2, w_2);

        // Append c_1, c_2 to the transcript
        transcript.append_point(c_1);
        transcript.append_point(c_2);

        // Get the challenge 𝛄
        let gamma = transcript.challenge_fe("gamma".to_string());

        // Compute s = witness + r * 𝛄
        let s = Self::add_vecs(&witness, &Self::scale_vec(&r, &gamma));

        let sigma_A = w_A + gamma * w_A_prime;
        let sigma_B = w_B + gamma * w_B_prime;
        let sigma_C = w_C + gamma * w_C_prime;

        let sigma_o = w_C + gamma * w_1 + gamma * gamma * w_2;

        NARKProof {
            c_A,
            c_B,
            c_C,
            c_A_prime,
            c_B_prime,
            c_C_prime,
            c_1,
            c_2,
            s,
            sigma_A,
            sigma_B,
            sigma_C,
            sigma_o,
        }
    }

    pub fn verify(
        r1cs: &R1CS<C::ScalarExt>,
        proof: &NARKProof<C>,
        pub_input: &[C::ScalarExt],
        transcript: &mut PoseidonTranscript<C>,
    ) {
        // Construct the vector z_prime = (pub_input, s)
        let z_prime = R1CS::construct_z(&proof.s, pub_input);

        // Compute the matrix-vector product of all three matrices and the vector z_prime
        let s_A = r1cs.A.mul_vector(&z_prime);
        let s_B = r1cs.B.mul_vector(&z_prime);
        let s_C = r1cs.C.mul_vector(&z_prime);

        // Initialize the Pedersen committer
        let pedersen = Pedersen::<C>::new(r1cs.z_len());

        // Append the commitments to the transcript
        transcript.append_point(proof.c_A);
        transcript.append_point(proof.c_B);
        transcript.append_point(proof.c_C);
        transcript.append_point(proof.c_A_prime);
        transcript.append_point(proof.c_B_prime);
        transcript.append_point(proof.c_C_prime);
        transcript.append_point(proof.c_1);
        transcript.append_point(proof.c_2);

        // Get the challenge 𝛄
        let gamma = transcript.challenge_fe("gamma".to_string());

        let lhs = proof.c_A + proof.c_A_prime * gamma;
        let rhs = pedersen.commit(&s_A, proof.sigma_A);
        assert_eq!(lhs, rhs, "c_A + c_A_prime * gamma != sigma_A");

        let lhs = proof.c_B + proof.c_B_prime * gamma;
        let rhs = pedersen.commit(&s_B, proof.sigma_B);
        assert_eq!(lhs, rhs, "c_B + c_B_prime * gamma != sigma_B");

        let lhs = proof.c_C + proof.c_C_prime * gamma;
        let rhs = pedersen.commit(&s_C, proof.sigma_C);
        assert_eq!(lhs, rhs, "c_C + c_C_prime * gamma != sigma_C");

        let lhs = proof.c_C + proof.c_1 * gamma + proof.c_2 * gamma * gamma;
        let rhs = pedersen.commit(&Self::hadamard(&s_A, &s_B), proof.sigma_o);
        assert_eq!(lhs, rhs, "c_C + c_1 * gamma + c_2 * gamma^2 != sigma_o");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use halo2curves::secp256k1::{Fq, Secp256k1Affine};
    use poseidon::sponge::IOPattern;

    #[test]
    fn test_nark() {
        let num_vars = 2usize.pow(4);
        let num_input = 5;
        let (r1cs, witness, public_input) = R1CS::<Fq>::produce_synthetic_r1cs(num_vars, num_input);

        let mut prover_transcript =
            PoseidonTranscript::<Secp256k1Affine>::new(b"test", IOPattern(vec![]));

        let proof =
            NARK::<Secp256k1Affine>::prove(&r1cs, &witness, &public_input, &mut prover_transcript);

        let mut verifier_transcript =
            PoseidonTranscript::<Secp256k1Affine>::new(b"test", IOPattern(vec![]));

        NARK::<Secp256k1Affine>::verify(&r1cs, &proof, &public_input, &mut verifier_transcript);
    }
}
