package com.oktsrl.math.impls.blas;

import jdistlib.rng.MersenneTwister;
import jdistlib.rng.RandomEngine;

import org.apache.commons.math3.stat.correlation.Covariance;
import org.jblas.Decompose;
import org.jblas.DoubleMatrix;

import com.oktsrl.math.Indices;
import com.oktsrl.math.MatrixOKT;
import com.oktsrl.math.Multiplication;
import com.oktsrl.math.Summation;
import com.oktsrl.math.impls.AbstractMatrixFactoryOKT;

public final class BLASMatrixFactoryOKT extends AbstractMatrixFactoryOKT {

	class BlasMultiplication implements Multiplication {
		DoubleMatrix acc;

		@Override
		public MatrixOKT getResult() {
			if (acc == null)
				return null;
			return new BlasMatrixOKT(acc, BLASMatrixFactoryOKT.this);
		}

		@Override
		public void mul(MatrixOKT... mats) {
			int index = 0;
			if (acc == null) {
				acc = ((BlasMatrixOKT) mats[0]).matrix.dup();
				index = 1;
			}
			for (int j = index; j < mats.length; j++)
				acc.muli(((BlasMatrixOKT) mats[j]).matrix);
		}

		@Override
		public void reset() {
			acc = null;
		}
	}

	class BlasSummation implements Summation {
		DoubleMatrix acc;

		@Override
		public void add(MatrixOKT... mats) {
			int index = 0;
			if (acc == null) {
				acc = ((BlasMatrixOKT) mats[0]).matrix.dup();
				index = 1;
			}
			for (int j = index; j < mats.length; j++)
				acc.addi(((BlasMatrixOKT) mats[j]).matrix);
		}

		@Override
		public MatrixOKT getResult() {
			if (acc == null)
				return null;
			return new BlasMatrixOKT(acc, BLASMatrixFactoryOKT.this);
		}

		@Override
		public void reset() {
			acc = null;
		}
	}

	private static final long serialVersionUID = 5864195360672832126L;

	@SuppressWarnings("unused")
	private final RandomEngine randomEngine = new MersenneTwister();

	@Override
	public MatrixOKT cholesky(MatrixOKT mat) {
		return new BlasMatrixOKT(
				Decompose.cholesky(((BlasMatrixOKT) mat).matrix), this);
	}

	@Override
	public MatrixOKT covariance(MatrixOKT mat) {
		final DoubleMatrix cov = new DoubleMatrix(new Covariance(
				((BlasMatrixOKT) mat).matrix.toArray2()).getCovarianceMatrix()
				.getData());
		return new BlasMatrixOKT(cov, this);
	}

	@Override
	public MatrixOKT create(double[][] values) {
		return new BlasMatrixOKT(new DoubleMatrix(values), this);
	}

	@Override
	public MatrixOKT create(int size) {
		return new BlasMatrixOKT(new DoubleMatrix(size, size), this);
	}

	@Override
	public MatrixOKT create(int rows, int columns) {
		return new BlasMatrixOKT(new DoubleMatrix(rows, columns), this);
	}

	@Override
	public Indices createIndices(int... indices) {
		return new BlasIndices(indices);
	}

	@Override
	public MatrixOKT eye(int size) {
		return new BlasMatrixOKT(DoubleMatrix.eye(size), this);
	}

	@Override
	public Multiplication getMultiplication() {
		return new BlasMultiplication();
	}

	@Override
	public Summation getSummation() {
		return new BlasSummation();
	}

	@Override
	public MatrixOKT mean(MatrixOKT mat) {
		return new BlasMatrixOKT(((BlasMatrixOKT) mat).matrix.columnMeans(),
				this);
	}

	@Override
	public double mean(MatrixOKT mat, double target, int policy) {
		final DoubleMatrix matrix = ((BlasMatrixOKT) mat).matrix;
		final Indices[] indices = mat.find(target, policy);
		final int n = indices[0].count();
		double mean = 0;
		for (int j = 0; j < n; j++)
			mean += matrix.get(indices[0].get(j), indices[1].get(j));
		return mean / n;
	}

	@Override
	public MatrixOKT ones(int size) {
		return new BlasMatrixOKT(DoubleMatrix.ones(size, size), this);
	}

	@Override
	public MatrixOKT ones(int rows, int columns) {
		return new BlasMatrixOKT(DoubleMatrix.ones(rows, columns), this);
	}

	@Override
	public MatrixOKT randn(int size) {
		return new BlasMatrixOKT(DoubleMatrix.randn(size, size), this);
	}

	@Override
	public MatrixOKT randn(int rows, int columns) {
		return new BlasMatrixOKT(DoubleMatrix.randn(rows, columns), this);
	}

	@Override
	public MatrixOKT sparse(double[][] values) {
		return new BlasMatrixOKT(new DoubleMatrix(values), this);
	}

	@Override
	public MatrixOKT sparse(int size) {
		return zeros(size, size);
	}

	@Override
	public MatrixOKT sparse(int rows, int columns) {
		return zeros(rows, columns);
	}

	@Override
	public MatrixOKT zeros(int size) {
		return new BlasMatrixOKT(DoubleMatrix.zeros(size, size), this);
	}

	@Override
	public MatrixOKT zeros(int rows, int columns) {
		return new BlasMatrixOKT(DoubleMatrix.zeros(rows, columns), this);
	}
}