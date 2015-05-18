package com.oktsrl.math;

import java.io.Serializable;

public interface MatrixOKT extends Serializable {
	public static final int GREATER = 0x01;
	public static final int LESS = 0x02;
	public static final int EQUAL = 0x04;
	public static final int NOT = 0x08;

	public MatrixOKT abs();

	public MatrixOKT columns(boolean linked, Indices indices);

	public MatrixOKT columns(boolean linked, int... columns);

	public MatrixOKT columns(Indices indices);

	public MatrixOKT columns(int... columns);

	public int columnsCount();

	public double det();

	public MatrixOKT div(double value);

	public MatrixOKT div(MatrixOKT mat);

	public MatrixOKT divMe(double value);

	public MatrixOKT divMe(MatrixOKT mat);

	public double dot(MatrixOKT mat);

	public double dotColumnColumn(int row, int column);

	public double dotColumnColumn(int row, MatrixOKT mat, int column);

	public double dotColumnRow(int row, int column);

	public double dotColumnRow(int row, MatrixOKT mat, int column);

	public double dotRowColumn(int row, int column);

	public double dotRowColumn(int row, MatrixOKT mat, int column);

	public double dotRowRow(int row, int column);

	public double dotRowRow(int row, MatrixOKT mat, int column);

	public Indices[] find(double target, int policy);

	public Indices findColumnIndices(double target, int row, int policy);

	public Indices findIndices(double target, int policy);

	public Indices findRowIndices(double target, int column, int policy);

	public double get(int row, int column);

	public MatrixOKT getCopy();

	public MatrixFactoryOKT getFactory();

	public MatrixOKT inverse();

	public boolean isEmpty();

	public double max();

	public double min();

	public MatrixOKT mul(double value);

	public MatrixOKT mul(MatrixOKT mat);

	public MatrixOKT mulMe(double value);

	public MatrixOKT mulMe(MatrixOKT mat);

	public int nnz();

	public void print(String name);

	public void put(Indices indices, double value);

	public MatrixOKT put(Indices rowIndices, Indices columnIndices,
			MatrixOKT mat);

	public MatrixOKT putColumn(int column, MatrixOKT mat);

	public MatrixOKT putRow(int row, MatrixOKT mat);

	public MatrixOKT replace(double target, double value);

	public MatrixOKT rows(boolean linked, Indices indices);

	public MatrixOKT rows(boolean linked, int... rows);

	public MatrixOKT rows(Indices indices);

	public MatrixOKT rows(int... rows);

	public int rowsCount();

	public void set(int row, int column, double value);

	public MatrixOKT sub(double value);

	public MatrixOKT sub(MatrixOKT mat);

	public MatrixOKT submatrix(Indices rowIndices, Indices columnIndices);

	public MatrixOKT subMe(double value);

	public MatrixOKT subMe(MatrixOKT mat);

	public MatrixOKT sum(double value);

	public MatrixOKT sum(MatrixOKT mat);

	public MatrixOKT sumMe(double value);

	public MatrixOKT sumMe(MatrixOKT mat);

	public MatrixOKT transpose();

	public MatrixOKT transpose(boolean link);
}