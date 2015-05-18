package com.oktsrl.math.impls;

import org.jblas.util.Permutations;
import org.jblas.util.Random;

import com.oktsrl.math.Indices;
import com.oktsrl.math.MatrixFactoryOKT;
import com.oktsrl.math.MatrixOKT;

public abstract class AbstractMatrixFactoryOKT implements MatrixFactoryOKT {

	private static final long serialVersionUID = -1333824378110987687L;

	public Indices indices(int start, int step, int end) {
		final int n = (end - start) / step + 1;
		final int[] i = new int[n];
		for (int j = 0, ii = start; j < n; j++, ii += step)
			i[j] = ii;
		return createIndices(i);
	}

	@Override
	public int nextRandomInt(int max) {
		return Random.nextInt(max);
	}

	@Override
	public int[] randperm(int size) {
		return Permutations.randomPermutation(size);
	}

	@Override
	public double sigmoid(double x) {
		return 1d / (1d + Math.exp(-x));
	}

	@Override
	public MatrixOKT wishart(MatrixOKT Sigma, int df) {
		return wishart(Sigma, df, null);
	}

	@Override
	public MatrixOKT wishart(MatrixOKT Sigma, int df, MatrixOKT D) {
		if (D == null)
			D = cholesky(Sigma);
		if (D == null)
			return null;
		final int p = D.rowsCount();
		final MatrixOKT Z = randn(df, p).mul(D);
		return Z.transpose().mul(Z);
	}
}