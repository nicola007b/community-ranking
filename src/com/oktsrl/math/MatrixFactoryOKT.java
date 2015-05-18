package com.oktsrl.math;

import java.io.Serializable;

public interface MatrixFactoryOKT extends Serializable {

	public MatrixOKT cholesky(MatrixOKT mat);

	public MatrixOKT covariance(MatrixOKT mat);

	public MatrixOKT create(double[][] values);

	public MatrixOKT create(int size);

	public MatrixOKT create(int rows, int columns);

	public Indices createIndices(int... indices);

	public MatrixOKT eye(int size);

	public Multiplication getMultiplication();

	public Summation getSummation();

	public MatrixOKT mean(MatrixOKT mat);

	public double mean(MatrixOKT mat, double target, int policy);

	public int nextRandomInt(int max);

	public MatrixOKT ones(int size);

	public MatrixOKT ones(int rows, int columns);

	public MatrixOKT randn(int size);

	public MatrixOKT randn(int rows, int columns);

	public int[] randperm(int size);

	public double sigmoid(double x);

	public MatrixOKT sparse(double[][] values);

	public MatrixOKT sparse(int size);

	public MatrixOKT sparse(int rows, int columns);

	public MatrixOKT wishart(MatrixOKT Sigma, int df);

	public MatrixOKT wishart(MatrixOKT Sigma, int df, MatrixOKT D);

	public MatrixOKT zeros(int size);

	// public MatrixOKT chol(MatrixOKT a);
	public MatrixOKT zeros(int rows, int columns);
}