package com.oktsrl.math.impls.ujmp;

import java.util.Iterator;

import jdistlib.rng.MersenneTwister;
import jdistlib.rng.RandomEngine;

import org.ujmp.core.Matrix;
import org.ujmp.core.calculation.Calculation;
import org.ujmp.core.matrix.SparseMatrix;

import com.oktsrl.math.Indices;
import com.oktsrl.math.MatrixOKT;
import com.oktsrl.math.Multiplication;
import com.oktsrl.math.Summation;
import com.oktsrl.math.impls.AbstractMatrixFactoryOKT;

public class UJMPMatrixFactoryOKT extends AbstractMatrixFactoryOKT {

	class UJMPMultiplication implements Multiplication {
		Matrix acc;

		@Override
		public MatrixOKT getResult() {
			if (acc == null)
				return null;
			return new UJMPMatrixOKT(acc, UJMPMatrixFactoryOKT.this);
		}

		@Override
		public void mul(MatrixOKT... mats) {
			int index = 0;
			if (acc == null) {
				acc = getMatrixClone(((UJMPMatrixOKT) mats[0]).matrix);
				index = 1;
			}
			for (int j = index; j < mats.length; j++)
				acc.mtimes(Calculation.Ret.ORIG, false,
						((UJMPMatrixOKT) mats[j]).matrix);
		}

		@Override
		public void reset() {
			acc = null;
		}
	}

	class UJMPSummation implements Summation {
		Matrix acc;

		@Override
		public void add(MatrixOKT... mats) {
			int index = 0;
			if (acc == null) {
				acc = getMatrixClone(((UJMPMatrixOKT) mats[0]).matrix);
				index = 1;
			}
			for (int j = index; j < mats.length; j++)
				acc.plus(Calculation.Ret.ORIG, false,
						((UJMPMatrixOKT) mats[j]).matrix);
		}

		@Override
		public MatrixOKT getResult() {
			if (acc == null)
				return null;
			return new UJMPMatrixOKT(acc, UJMPMatrixFactoryOKT.this);
		}

		@Override
		public void reset() {
			acc = null;
		}
	}

	private static final long serialVersionUID = 5828217981816793997L;

	protected static Matrix getMatrixClone(Matrix matrix) {
		if (matrix.isSparse()) {
			final Matrix ret = SparseMatrix.factory.zeros(matrix.getRowCount(),
					matrix.getColumnCount());
			final Iterator<long[]> iterator = matrix.availableCoordinates()
					.iterator();
			long[] coordinate;
			while (iterator.hasNext()) {
				coordinate = iterator.next();
				ret.setAsDouble(matrix.getAsDouble(coordinate), coordinate);
			}
			return ret;
		}
		return matrix.clone();
	}

	@SuppressWarnings("unused")
	private final RandomEngine randomEngine = new MersenneTwister();

	@Override
	public MatrixOKT cholesky(MatrixOKT mat) {
		return new UJMPMatrixOKT(((UJMPMatrixOKT) mat).matrix.chol()
				.transpose(), this);
	}

	@Override
	public MatrixOKT covariance(MatrixOKT mat) {
		return new UJMPMatrixOKT(((UJMPMatrixOKT) mat).matrix.cov(
				Calculation.Ret.NEW, false), this);
	}

	@Override
	public MatrixOKT create(double[][] values) {
		final Matrix mat = SparseMatrix.factory.zeros(values.length,
				values[0].length);
		for (int r = 0; r < values.length; r++)
			for (int c = 0; c < values[0].length; c++)
				mat.setAsDouble(values[r][c], r, c);

		return new UJMPMatrixOKT(mat, this);
	}

	@Override
	public MatrixOKT create(int size) {
		return zeros(size, size);
	}

	@Override
	public MatrixOKT create(int rows, int columns) {
		return zeros(rows, columns);
	}

	@Override
	public Indices createIndices(int... indices) {
		return new UJMPIndices(indices);
	}

	@Override
	public MatrixOKT eye(int size) {
		return new UJMPMatrixOKT(SparseMatrix.factory.eye(size, size), this);
	}

	@Override
	public Multiplication getMultiplication() {
		return new UJMPMultiplication();
	}

	@Override
	public Summation getSummation() {
		return new UJMPSummation();
	}

	@Override
	public MatrixOKT mean(MatrixOKT mat) {
		return new UJMPMatrixOKT(((UJMPMatrixOKT) mat).matrix.mean(
				Calculation.Ret.NEW, 0, false), this);
	}

	@Override
	public double mean(MatrixOKT mat, double target, int policy) {
		final Matrix matrix = ((UJMPMatrixOKT) mat).matrix;
		final Indices[] indices = mat.find(target, policy);
		final int n = indices[0].count();
		double mean = 0;
		for (int j = 0; j < n; j++)
			mean += matrix.getAsDouble(indices[0].get(j), indices[1].get(j));
		return mean / n;
	}

	@Override
	public MatrixOKT ones(int size) {
		return new UJMPMatrixOKT(SparseMatrix.factory.ones(size, size), this);
	}

	@Override
	public MatrixOKT ones(int rows, int columns) {
		return new UJMPMatrixOKT(SparseMatrix.factory.ones(rows, columns), this);
	}

	@Override
	public MatrixOKT randn(int size) {
		return new UJMPMatrixOKT(SparseMatrix.factory.randn(size, size), this);
	}

	@Override
	public MatrixOKT randn(int rows, int columns) {
		return new UJMPMatrixOKT(SparseMatrix.factory.randn(rows, columns),
				this);
	}

	@Override
	public MatrixOKT sparse(double[][] values) {
		final Matrix mat = SparseMatrix.factory.zeros(values.length,
				values[0].length);
		for (int r = 0; r < values.length; r++)
			for (int c = 0; c < values[0].length; c++)
				mat.setAsDouble(values[r][c], r, c);

		return new UJMPMatrixOKT(mat, this);
	}

	@Override
	public MatrixOKT sparse(int size) {
		return new UJMPMatrixOKT(SparseMatrix.factory.zeros(size, size), this);
	}

	@Override
	public MatrixOKT sparse(int rows, int columns) {
		return new UJMPMatrixOKT(SparseMatrix.factory.zeros(rows, columns),
				this);
	}

	@Override
	public MatrixOKT zeros(int size) {
		return new UJMPMatrixOKT(SparseMatrix.factory.zeros(size, size), this);
	}

	@Override
	public MatrixOKT zeros(int rows, int columns) {
		return new UJMPMatrixOKT(SparseMatrix.factory.zeros(rows, columns),
				this);
	}
}