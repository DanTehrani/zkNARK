use crate::traits::CurveLike;
use halo2curves::group::Group;

pub struct Pedersen<C: CurveLike> {
    G: Vec<C::CurveExt>,
    H: C::CurveExt,
}

impl<C: CurveLike> Pedersen<C> {
    pub fn new(n: usize) -> Self {
        // Compute the basis G and H
        // TODO: This is absolutely not secure!
        let g = C::generator();
        let G = (0..n)
            .map(|i| g * C::ScalarExt::from(i as u64))
            .collect::<Vec<C::CurveExt>>();
        let H = g * C::ScalarExt::from(33);

        Self { G, H }
    }

    pub fn commit(&self, x: &[C::ScalarExt], blinder: C::ScalarExt) -> C::CurveExt {
        assert_eq!(self.G.len(), x.len());

        let pairs = self
            .G
            .iter()
            .enumerate()
            .map(|(i, g)| (x[i], (*g).into()))
            .collect::<Vec<(C::ScalarExt, C::Curve)>>();

        // TODO: Use multiexp
        let mut com = pairs
            .iter()
            .map(|(x, g)| *g * *x)
            .fold(C::CurveExt::identity(), |acc, x| acc + x);

        // com += multiexp(&pairs);
        com += self.H * blinder;

        com.into()
    }
}
