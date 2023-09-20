use core::panic;
use std::collections::BTreeMap;

use halo2curves::group::ff::PrimeField;
use halo2curves::group::Curve;
use poseidon::sponge::{IOPattern, PoseidonSponge};

use crate::traits::CurveLike;

#[derive(Clone)]
pub struct PoseidonTranscript<C: CurveLike> {
    sponge: PoseidonSponge<C::ScalarExt, 3>,
    challenges: BTreeMap<String, C::ScalarExt>, // We store challenges for later reference
}

impl<C: CurveLike> PoseidonTranscript<C> {
    pub fn new(label: &'static [u8], io_pattern: IOPattern) -> Self {
        Self {
            sponge: PoseidonSponge::new(label, io_pattern),
            challenges: BTreeMap::new(),
        }
    }
}

impl<C: CurveLike> PoseidonTranscript<C> {
    pub fn append_fe(&mut self, fe: C::ScalarExt) {
        self.sponge.absorb(&[fe]);
    }

    pub fn append_point(&mut self, p: C::CurveExt) {
        let (x, y) = p.to_affine().coordinates();

        let x_bytes = x.to_repr();
        let y_bytes = y.to_repr();

        self.append_bytes(x_bytes);
        self.append_bytes(y_bytes);
    }

    pub fn append_bytes(&mut self, bytes: [u8; 32]) {
        let mut bytes_low = bytes[0..16].to_vec();
        let mut bytes_high = bytes[16..32].to_vec();

        bytes_low.resize(32, 0);
        bytes_high.resize(32, 0);

        let fe_low = C::ScalarExt::from_repr(bytes_low.try_into().unwrap()).unwrap();
        let fe_high = C::ScalarExt::from_repr(bytes_high.try_into().unwrap()).unwrap();
        self.sponge.absorb(&[fe_low, fe_high]);
    }

    pub fn challenge_fe(&mut self, label: String) -> C::ScalarExt {
        let c = self.sponge.squeeze(1)[0];
        if label != "".to_string() {
            if self.challenges.contains_key(&label) {
                panic!("Challenge label {} already exists", label);
            }
            self.challenges.insert(label, c);
        }

        c
    }
}
