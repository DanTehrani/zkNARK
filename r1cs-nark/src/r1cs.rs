use crate::halo2curves::FieldExt;

#[derive(Clone, Debug)]
pub struct SparseMatrixEntry<F: FieldExt> {
    pub row: usize,
    pub col: usize,
    pub val: F,
}

#[derive(Clone)]
pub struct Matrix<F: FieldExt> {
    pub entries: Vec<SparseMatrixEntry<F>>,
    pub num_cols: usize,
    pub num_rows: usize,
}

impl<F: FieldExt> Matrix<F> {
    pub const fn empty() -> Self {
        Self {
            entries: vec![],
            num_cols: 0,
            num_rows: 0,
        }
    }

    pub fn new(entries: Vec<SparseMatrixEntry<F>>, num_cols: usize, num_rows: usize) -> Self {
        Self {
            entries,
            num_cols,
            num_rows,
        }
    }

    pub fn mul_vector(&self, vec: &[F]) -> Vec<F> {
        debug_assert_eq!(vec.len(), self.num_cols);
        let mut result = vec![F::zero(); self.num_rows];
        let entries = &self.entries;
        for i in 0..entries.len() {
            let row = entries[i].row;
            let col = entries[i].col;
            let val = entries[i].val;
            result[row] += val * vec[col];
        }
        result
    }
}

#[derive(Clone)]
pub struct R1CS<F: FieldExt> {
    pub A: Matrix<F>,
    pub B: Matrix<F>,
    pub C: Matrix<F>,
    pub num_vars: usize,
    pub num_input: usize,
}

impl<F: FieldExt> R1CS<F> {
    pub const fn empty() -> Self {
        Self {
            A: Matrix::empty(),
            B: Matrix::empty(),
            C: Matrix::empty(),
            num_vars: 0,
            num_input: 0,
        }
    }

    pub fn hadamard_prod(a: &[F], b: &[F]) -> Vec<F> {
        assert_eq!(a.len(), b.len());
        let mut result = vec![F::zero(); a.len()];
        for i in 0..a.len() {
            result[i] = a[i] * b[i];
        }
        result
    }

    pub fn num_cons(&self) -> usize {
        self.A.entries.len()
    }

    pub fn z_len(&self) -> usize {
        self.num_vars.next_power_of_two() * 2
    }

    pub fn construct_z(witness: &[F], public_input: &[F]) -> Vec<F> {
        assert!(witness.len() >= public_input.len());
        // Z = (1, io, witness)
        let n = witness.len().next_power_of_two();
        let mut z = vec![];
        z.push(F::one());
        z.extend_from_slice(public_input);
        z.resize(n, F::zero());

        z.extend_from_slice(witness);
        z.resize(n * 2, F::zero());

        z
    }

    pub fn produce_synthetic_r1cs(num_vars: usize, num_input: usize) -> (Self, Vec<F>, Vec<F>) {
        let mut public_input = Vec::with_capacity(num_input);
        let mut witness = Vec::with_capacity(num_vars);

        for i in 0..num_input {
            public_input.push(F::from((i + 1) as u64));
        }

        for i in 0..num_vars {
            witness.push(F::from((i + 1) as u64));
        }

        let z = Self::construct_z(&witness, &public_input);

        let mut A_entries: Vec<SparseMatrixEntry<F>> = vec![];
        let mut B_entries: Vec<SparseMatrixEntry<F>> = vec![];
        let mut C_entries: Vec<SparseMatrixEntry<F>> = vec![];

        // Constrain the variables
        let witness_start_index = num_vars.next_power_of_two();
        for i in witness_start_index..(witness_start_index + num_vars) {
            let A_col = i;
            let B_col = (i + 1) % (witness_start_index + num_vars);
            let C_col = (i + 2) % (witness_start_index + num_vars);

            // For the i'th constraint,
            // add the value 1 at the (i % num_vars)th column of A, B.
            // Compute the corresponding C_column value so that A_i * B_i = C_i
            // we apply multiplication since the Hadamard product is computed for Az ・ Bz,

            // We only _enable_ a single variable in each constraint.
            let AB = if z[C_col] == F::zero() {
                F::zero()
            } else {
                F::one()
            };

            A_entries.push(SparseMatrixEntry {
                row: i,
                col: A_col,
                val: AB,
            });
            B_entries.push(SparseMatrixEntry {
                row: i,
                col: B_col,
                val: AB,
            });
            C_entries.push(SparseMatrixEntry {
                row: i,
                col: C_col,
                val: if z[C_col] == F::zero() {
                    F::zero()
                } else {
                    (z[A_col] * z[B_col]) * z[C_col].invert().unwrap()
                },
            });
        }

        // Constrain the public inputs
        let input_index_start = 1;
        for i in input_index_start..(input_index_start + num_input) {
            let A_col = i;
            let B_col = (i + 1) % input_index_start + num_input;
            let C_col = (i + 2) % input_index_start + num_input;

            // For the i'th constraint,
            // add the value 1 at the (i % num_vars)th column of A, B.
            // Compute the corresponding C_column value so that A_i * B_i = C_i
            // we apply multiplication since the Hadamard product is computed for Az ・ Bz,

            // We only _enable_ a single variable in each constraint.
            let AB = if z[C_col] == F::zero() {
                F::zero()
            } else {
                F::one()
            };

            A_entries.push(SparseMatrixEntry {
                row: i,
                col: A_col,
                val: AB,
            });
            B_entries.push(SparseMatrixEntry {
                row: i,
                col: B_col,
                val: AB,
            });
            C_entries.push(SparseMatrixEntry {
                row: i,
                col: C_col,
                val: if z[C_col] == F::zero() {
                    F::zero()
                } else {
                    (z[A_col] * z[B_col]) * z[C_col].invert().unwrap()
                },
            });
        }

        let num_cols = z.len();
        let num_rows = z.len();

        let A = Matrix::new(A_entries, num_cols, num_rows);
        let B = Matrix::new(B_entries, num_cols, num_rows);
        let C = Matrix::new(C_entries, num_cols, num_rows);

        (
            Self {
                A,
                B,
                C,
                num_vars,
                num_input,
            },
            witness,
            public_input,
        )
    }

    pub fn is_sat(&self, witness: &[F], public_input: &[F]) -> bool {
        let z = Self::construct_z(witness, public_input);
        let Az = self.A.mul_vector(&z);
        let Bz = self.B.mul_vector(&z);
        let Cz = self.C.mul_vector(&z);

        Self::hadamard_prod(&Az, &Bz) == Cz
    }
}
