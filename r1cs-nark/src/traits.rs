use crate::halo2curves::{group::prime::PrimeCurveAffine, CurveAffine, CurveExt, FieldExt};
use poseidon::{halo2curves::secp256k1::Secp256k1Affine, traits::PoseidonField};
use std::ops::{Add, Sub};

pub trait CurveLike:
    PrimeCurveAffine<Scalar = <Self as CurveLike>::ScalarExt, Curve = <Self as CurveLike>::CurveExt>
    + Default
    + Add<Output = <Self as PrimeCurveAffine>::Curve>
    + Sub<Output = <Self as PrimeCurveAffine>::Curve>
    + From<<Self as PrimeCurveAffine>::Curve>
{
    type ScalarExt: PoseidonField + FieldExt<Repr = [u8; 32]>;
    type CurveExt: CurveExt<AffineExt = Self, ScalarExt = <Self as CurveLike>::ScalarExt>;
    type Base: FieldExt<Repr = [u8; 32]>;

    fn coordinates(&self) -> (Self::Base, Self::Base);
}

impl CurveLike for Secp256k1Affine {
    type ScalarExt = <Secp256k1Affine as CurveAffine>::ScalarExt;
    type CurveExt = <Secp256k1Affine as CurveAffine>::CurveExt;
    type Base = <Secp256k1Affine as CurveAffine>::Base;

    fn coordinates(&self) -> (Self::Base, Self::Base) {
        let coords = <Self as CurveAffine>::coordinates(self).unwrap();
        (*coords.x(), *coords.y())
    }
}
