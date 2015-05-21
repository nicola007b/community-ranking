package com.oktsrl.math.impls.blas;

import org.jblas.Decompose;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.jblas.Solve;

import com.oktsrl.math.Indices;
import com.oktsrl.math.MatrixFactoryOKT;
import com.oktsrl.math.MatrixOKT;

public final class BlasMatrixOKT implements MatrixOKT {

	private static final long serialVersionUID = 393145071641567401L;

	public DoubleMatrix matrix;
	protected MatrixFactoryOKT factory;

	public BlasMatrixOKT(DoubleMatrix matrix, MatrixFactoryOKT factory) {
		this.matrix = matrix;
		this.factory = factory;
	}

	public BlasMatrixOKT(MatrixFactoryOKT factory) {
		this.factory = factory;

	}

	@Override
	public MatrixOKT abs() {
		return new BlasMatrixOKT(MatrixFunctions.abs(matrix), factory);
	}

	@Override
	public MatrixOKT columns(boolean linked, Indices indices) {
		return columns(indices);
	}

	@Override
	public MatrixOKT columns(boolean linked, int... columns) {
		return columns(columns);
	}

	@Override
	public MatrixOKT columns(Indices indices) {
		return new BlasMatrixOKT(matrix.getColumns(indices.toArray()), factory);
	}

	@Override
	public MatrixOKT columns(int... columns) {
		return columns(factory.createIndices(columns));
	}

	@Override
	public int columnsCount() {
		return matrix.columns;
	}

	@Override
	public double det() {
		final Decompose.LUDecomposition<DoubleMatrix> lud = Decompose
				.lu(matrix);
		final int[] indices = lud.p.eq(1).findIndices();
		int positive = 0;
		int negative = 0;
		for (final int index : indices)
			if (index % 2 == 1)
				positive++;
			else
				negative++;
		return (positive > negative ? +1 : -1) * lud.u.diag().prod();
	}

	@Override
	public MatrixOKT div(double value) {
		return new BlasMatrixOKT(matrix.div(value), factory);
	}

	@Override
	public MatrixOKT div(MatrixOKT mat) {
		return new BlasMatrixOKT(matrix.mmul(Solve
				.pinv(((BlasMatrixOKT) mat).matrix)), factory);
	}

	@Override
	public MatrixOKT divMe(double value) {
		matrix.divi(value);
		return this;
	}

	@Override
	public MatrixOKT divMe(MatrixOKT mat) {
		matrix.mmuli(Solve.pinv(((BlasMatrixOKT) mat).matrix));
		return this;
	}

	@Override
	public double dot(MatrixOKT mat) {
		return matrix.dot(((BlasMatrixOKT) mat).matrix);
	}

	@Override
	public double dotColumnColumn(int column1, int column2) {
		return matrix.getColumn(column1).dot(matrix.getRow(column2));
	}

	@Override
	public double dotColumnColumn(int column1, MatrixOKT mat, int column2) {
		return matrix.getColumn(column1).dot(
				((BlasMatrixOKT) mat).matrix.getRow(column2));
	}

	@Override
	public double dotColumnRow(int row, int column) {
		return dotColumnRow(row, this, column);
	}

	@Override
	public double dotColumnRow(int column, MatrixOKT mat, int row) {
		DoubleMatrix m1, m2;
		if (matrix.columns == 1 || matrix.rows == 1)
			m1 = matrix;
		else
			m1 = ((BlasMatrixOKT) columns(column)).matrix;
		if (((BlasMatrixOKT) mat).matrix.columns == 1
				|| ((BlasMatrixOKT) mat).matrix.rows == 1)
			m2 = ((BlasMatrixOKT) mat).matrix;
		else
			m2 = ((BlasMatrixOKT) mat.rows(row)).matrix;
		return m1.dot(m2);
	}

	@Override
	public double dotRowColumn(int row, int column) {
		return dotRowColumn(row, this, column);
	}

	@Override
	public double dotRowColumn(int row, MatrixOKT mat, int column) {
		DoubleMatrix m1, m2;
		if (matrix.columns == 1 || matrix.rows == 1)
			m1 = matrix;
		else
			m1 = ((BlasMatrixOKT) rows(row)).matrix;
		if (((BlasMatrixOKT) mat).matrix.columns == 1
				|| ((BlasMatrixOKT) mat).matrix.rows == 1)
			m2 = ((BlasMatrixOKT) mat).matrix;
		else
			m2 = ((BlasMatrixOKT) mat.columns(column)).matrix;
		return m1.dot(m2);
	}

	@Override
	public double dotRowRow(int row1, int row2) {
		return matrix.getRow(row1).dot(matrix.getRow(row2));
	}

	@Override
	public double dotRowRow(int row1, MatrixOKT mat, int row2) {
		return matrix.getRow(row1).dot(
				((BlasMatrixOKT) mat).matrix.getRow(row2));
	}

	@Override
	@SuppressWarnings("unused")
	public Indices[] find(double target, int policy) {
		// gli indici sono disposti per colonna
		final Indices indices = findIndices(target, policy);
		final int nr = matrix.rows;
		final int nc = matrix.columns;
		final int[] ir = new int[indices.count()];
		final int[] ic = new int[indices.count()];
		int pos = 0;
		for (final int index : indices.toArray()) {
			ir[pos] = index % nr;
			ic[pos++] = index / nr;
		}
		return new Indices[] { new BlasIndices(ir), new BlasIndices(ic) };
	}

	@Override
	public Indices findColumnIndices(double target, int row, int policy) {
		DoubleMatrix tmp;
		if ((policy & GREATER) == GREATER) {
			if ((policy & EQUAL) == EQUAL)
				tmp = matrix.getRow(row).ge(target);
			else
				tmp = matrix.getRow(row).gt(target);
		} else if ((policy & LESS) == LESS) {
			if ((policy & EQUAL) == EQUAL)
				tmp = matrix.getRow(row).le(target);
			else
				tmp = matrix.getRow(row).lt(target);
		} else if ((policy & EQUAL) == EQUAL)
			tmp = matrix.getRow(row).eq(target);
		else if ((policy & NOT) == NOT)
			tmp = matrix.getRow(row).ne(target);
		else
			return new BlasIndices(new int[0]);

		return new BlasIndices(tmp.findIndices());
	}

	@Override
	public Indices findIndices(double target, int policy) {
		int[] indices;
		if ((policy & GREATER) == GREATER) {
			if ((policy & EQUAL) == EQUAL)
				indices = matrix.ge(target).findIndices();
			else
				indices = matrix.gt(target).findIndices();
		} else if ((policy & LESS) == LESS) {
			if ((policy & EQUAL) == EQUAL)
				indices = matrix.le(target).findIndices();
			else
				indices = matrix.lt(target).findIndices();
		} else if ((policy & EQUAL) == EQUAL)
			indices = matrix.eq(target).findIndices();
		else if ((policy & NOT) == NOT)
			indices = matrix.ne(target).findIndices();
		else
			return new BlasIndices(new int[0]);

		return new BlasIndices(indices);
	}

	@Override
	public Indices findRowIndices(double target, int column, int policy) {
		DoubleMatrix tmp;
		if ((policy & GREATER) == GREATER) {
			if ((policy & EQUAL) == EQUAL)
				tmp = matrix.getColumn(column).ge(target);
			else
				tmp = matrix.getColumn(column).gt(target);
		} else if ((policy & LESS) == LESS) {
			if ((policy & EQUAL) == EQUAL)
				tmp = matrix.getColumn(column).le(target);
			else
				tmp = matrix.getColumn(column).lt(target);
		} else if ((policy & EQUAL) == EQUAL)
			tmp = matrix.getColumn(column).eq(target);
		else if ((policy & NOT) == NOT)
			tmp = matrix.getColumn(column).ne(target);
		else
			return new BlasIndices(new int[0]);
		return new BlasIndices(tmp.findIndices());
	}

	@Override
	public double get(int row, int column) {
		return matrix.get(row, column);
	}

	@Override
	public MatrixOKT getCopy() {
		// DoubleMatrix mat= new DoubleMatrix(matrix.rows, matrix.columns);
		return new BlasMatrixOKT(matrix.dup(), factory);
	}

	@Override
	public MatrixFactoryOKT getFactory() {
		return factory;
	}

	@Override
	public MatrixOKT inverse() {
		return new BlasMatrixOKT(Solve.pinv(matrix), factory);
	}

	@Override
	public boolean isEmpty() {
		return matrix.isEmpty();
	}

	@Override
	public double max() {
		return matrix.max();
	}

	@Override
	public double min() {
		return matrix.min();
	}

	@Override
	public MatrixOKT mul(double value) {
		return new BlasMatrixOKT(matrix.mul(value), factory);
	}

	@Override
	public MatrixOKT mul(MatrixOKT mat) {
		return new BlasMatrixOKT(matrix.mmul(((BlasMatrixOKT) mat).matrix),
				factory);
	}

	@Override
	public MatrixOKT mulMe(double value) {
		matrix.muli(value);
		return this;
	}

	@Override
	public MatrixOKT mulMe(MatrixOKT mat) {
		matrix.mmuli(((BlasMatrixOKT) mat).matrix);
		return this;
	}

	@Override
	public int nnz() {
		return matrix.ne(0).findIndices().length;
	}

	@Override
	public void print(String name) {
		System.out.println("Matrix " + name + ":");
		matrix.print();
	}

	@Override
	public void put(Indices indices, double value) {
		matrix.put(indices.toArray(), value);
	}

	@Override
	public MatrixOKT put(Indices rowIndices, Indices columnIndices,
			MatrixOKT mat) {
		matrix.put(rowIndices.toArray(), columnIndices.toArray(),
				((BlasMatrixOKT) mat).matrix);
		return this;
	}

	@Override
	public MatrixOKT putColumn(int column, MatrixOKT mat) {
		matrix.putColumn(column, ((BlasMatrixOKT) mat).matrix);
		return this;
	}

	@Override
	public MatrixOKT putRow(int row, MatrixOKT mat) {
		final DoubleMatrix m = ((BlasMatrixOKT) mat).matrix;
		// if( m.isRowVector() )
		// matrix.putRow(row, m.transpose());
		// else
		matrix.putRow(row, m);
		return this;
	}

	@Override
	public MatrixOKT replace(double target, double value) {
		matrix.put(matrix.eq(value).findIndices(), value);
		return this;
	}

	@Override
	public MatrixOKT rows(boolean linked, Indices indices) {
		return rows(indices);
	}

	@Override
	public MatrixOKT rows(boolean linked, int... rows) {
		return rows(rows);
	}

	@Override
	public MatrixOKT rows(Indices indices) {
		return new BlasMatrixOKT(matrix.getRows(indices.toArray()), factory);
	}

	@Override
	public MatrixOKT rows(int... rows) {
		return rows(factory.createIndices(rows));
	}

	@Override
	public int rowsCount() {
		return matrix.rows;
	}

	@Override
	public void set(int row, int column, double value) {
		matrix.put(row, column, value);
	}

	@Override
	public MatrixOKT sub(double value) {
		return new BlasMatrixOKT(matrix.sub(value), factory);
	}

	@Override
	public MatrixOKT sub(MatrixOKT mat) {
		return new BlasMatrixOKT(matrix.sub(((BlasMatrixOKT) mat).matrix),
				factory);
	}

	@Override
	public MatrixOKT submatrix(Indices rowIndices, Indices columnIndices) {
		return rows(rowIndices).columns(columnIndices);
	}

	@Override
	public MatrixOKT subMe(double value) {
		matrix.subi(value);
		return this;
	}

	@Override
	public MatrixOKT subMe(MatrixOKT mat) {
		matrix.subi(((BlasMatrixOKT) mat).matrix);
		return this;
	}

	@Override
	public MatrixOKT sum(double value) {
		return new BlasMatrixOKT(matrix.add(value), factory);
	}

	@Override
	public MatrixOKT sum(MatrixOKT mat) {
		return new BlasMatrixOKT(matrix.add(((BlasMatrixOKT) mat).matrix),
				factory);
	}

	@Override
	public MatrixOKT sumMe(double value) {
		matrix.addi(value);
		return this;
	}

	@Override
	public MatrixOKT sumMe(MatrixOKT mat) {
		matrix.addi(((BlasMatrixOKT) mat).matrix);
		return this;
	}

	@Override
	public MatrixOKT transpose() {
		return new BlasMatrixOKT(matrix.transpose(), factory);
	}

	@Override
	public MatrixOKT transpose(boolean link) {
		return transpose();
	}
}
