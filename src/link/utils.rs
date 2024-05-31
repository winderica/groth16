//! Utils for matrix and vector operations

use ark_ec::{pairing::Pairing, AffineRepr, CurveGroup, ScalarMul, VariableBaseMSM};
use ark_ff::Zero;
use ark_std::{
    marker::PhantomData,
    ops::{AddAssign, Mul},
    vec,
};
use rayon::prelude::*;

// a column is a vector of CoeffPos-s
type Col<T> = Vec<(T, usize)>;

// TODO: One could consider a cache-friendlier implementation for the 2-row case

/// Column-Major Sparse Matrix
#[derive(Clone, Debug)]
pub struct SparseMatrix<T> {
    cols: Vec<Col<T>>, // a vector of columns
    pub nr: usize,     // no. of rows
    pub nc: usize,     // no. of columns
}

impl<T: Copy> SparseMatrix<T> {
    // NB: Given column by column
    pub fn new(nr: usize, nc: usize) -> SparseMatrix<T> {
        SparseMatrix {
            cols: vec![vec![]; nc],
            nr,
            nc,
        }
    }

    /// Insert value `v` in the column index `c` at row index `r`
    pub fn insert_val(&mut self, r: usize, c: usize, v: T) {
        self.cols[c].push((v, r));
    }

    /// insert a continuous sequence of values at row r starting from c_offset
    pub fn insert_row_slice(&mut self, r: usize, c_offset: usize, vs: Vec<T>) {
        // NB: could be improved in efficiency by first extending the vector
        for (i, x) in vs.into_iter().enumerate() {
            self.insert_val(r, c_offset + i, x);
        }
    }
}

pub struct SparseLinAlgebra<PE: Pairing> {
    pairing_engine_type: PhantomData<PE>,
}

impl<PE: Pairing> SparseLinAlgebra<PE> {
    /// Inner products of all columns of a sparse matrix and another (sparse)
    /// vector to compute the matrix multiplication `m^T \dot v` where `m^T`
    /// is the transpose of `m`. v has dimensions `v.len() x 1` and m has
    /// dimensions `nr x nc`. Returns a matrix of dimension `nr x 1`
    pub fn sparse_vector_matrix_mult(
        v: &[PE::ScalarField],
        m: &SparseMatrix<PE::G1Affine>,
    ) -> Vec<PE::G1Affine> {
        PE::G1::normalize_batch(
            &cfg_iter!(m.cols)
                .map(|w| w.iter().map(|(g, i)| g.mul(v[*i])).sum::<PE::G1>())
                .collect::<Vec<_>>(),
        )
    }
}

/// MSM between a scalar vector and a G1 vector
pub fn inner_product<PE: Pairing>(a: &[PE::ScalarField], b: &[PE::G1Affine]) -> PE::G1Affine {
    PE::G1::msm_unchecked(b, &a).into_affine()
}

/// Scale given vector `v` by scalar `a`
pub fn scale_vector<PE: Pairing>(
    a: PE::ScalarField,
    v: &[PE::ScalarField],
) -> Vec<PE::ScalarField> {
    cfg_iter!(v).map(|x| a * x).collect()
}

/// Given a group element `g` and vector `multiples` of scalars, returns a
/// vector with elements `v_i * g`
pub fn multiples_of_g<G: AffineRepr>(g: &G, multiples: &[G::ScalarField]) -> Vec<G> {
    g.into_group().batch_mul(multiples)
}
