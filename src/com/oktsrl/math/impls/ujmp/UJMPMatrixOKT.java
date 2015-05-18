package com.oktsrl.math.impls.ujmp;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;

import org.ujmp.core.Matrix;
import org.ujmp.core.calculation.Calculation;

import com.oktsrl.math.Indices;
import com.oktsrl.math.MatrixFactoryOKT;
import com.oktsrl.math.MatrixOKT;
import com.oktsrl.math.impls.AbstractMatrixOKT;

public final class UJMPMatrixOKT extends AbstractMatrixOKT {

	interface Policy {
		boolean isTruth(double value);
	}

	private static final long serialVersionUID = 2006854950752167789L;

	public Matrix matrix;

	public UJMPMatrixOKT(Matrix matrix, MatrixFactoryOKT factory) {
		// TEST OK
		super(factory);
		this.matrix = matrix;
	}

	@Override
	public MatrixOKT abs() {
		// TEST OK
		return new UJMPMatrixOKT(matrix.abs(Calculation.Ret.NEW), factory);
	}

	@Override
	public MatrixOKT columns(boolean linked, Indices indices) {
		return new UJMPMatrixOKT(matrix.selectColumns(Calculation.Ret.LINK,
				indices.toLongArray()), factory);
	}

	@Override
	public MatrixOKT columns(boolean linked, int... columns) {
		final long[] tmp = new long[columns.length];
		for (int j = 0; j < tmp.length; j++)
			tmp[j] = columns[j];
		return new UJMPMatrixOKT(
				matrix.selectColumns(Calculation.Ret.LINK, tmp), factory);
	}

	@Override
	public MatrixOKT columns(Indices indices) {
		// TEST OK
		return new UJMPMatrixOKT(matrix.selectColumns(Calculation.Ret.NEW,
				indices.toLongArray()), factory);
	}

	@Override
	public MatrixOKT columns(int... columns) {
		// TEST OK
		final long[] tmp = new long[columns.length];
		for (int j = 0; j < tmp.length; j++)
			tmp[j] = columns[j];
		return new UJMPMatrixOKT(
				matrix.selectColumns(Calculation.Ret.NEW, tmp), factory);
	}

	@Override
	public int columnsCount() {
		// TEST OK
		return (int) matrix.getColumnCount();
	}

	private Indices convertToIndices(Matrix boolMatrix) {
		// TEST OK
		final Iterator<long[]> iterator = boolMatrix.nonZeroCoordinates()
				.iterator();
		final ArrayList<Number> indices = new ArrayList<Number>();
		long[] xy;
		while (iterator.hasNext()) {
			xy = iterator.next();
			if (boolMatrix.getAsBoolean(xy[0], xy[1]))
				indices.add(Math.max(xy[0], xy[1]));
		}
		return new UJMPIndices(indices);
	}

	@Override
	public double det() {
		// TEST OK
		return matrix.det();
	}

	@Override
	public MatrixOKT div(double value) {
		// TEST OK
		return new UJMPMatrixOKT(matrix.divide(value), factory);
	}

	@Override
	public MatrixOKT div(MatrixOKT mat) {
		// TEST OK
		return new UJMPMatrixOKT(matrix.mtimes(((UJMPMatrixOKT) mat).matrix
				.inv()), factory);
	}

	@Override
	public MatrixOKT divMe(double value) {
		// TEST OK
		matrix = matrix.divide(Calculation.Ret.ORIG, false, value);
		return this;
	}

	@Override
	public MatrixOKT divMe(MatrixOKT mat) {
		// TEST OK
		matrix = matrix.mtimes(Calculation.Ret.ORIG, false,
				((UJMPMatrixOKT) mat).matrix.inv());
		return this;
	}

	@Override
	public double dot(MatrixOKT mat) {
		// TEST OK
		final Matrix m = ((UJMPMatrixOKT) mat).matrix;
		if (matrix.getRowCount() == 1) {
			if (m.getRowCount() == 1)
				return matrix.mtimes(m.transpose(Calculation.Ret.LINK))
						.getAsDouble(0, 0);
			else if (m.getColumnCount() == 1)
				return matrix.mtimes(m).getAsDouble(0, 0);
		} else if (matrix.getColumnCount() == 1)
			if (m.getRowCount() == 1)
				return m.mtimes(matrix).getAsDouble(0, 0);
			else if (m.getColumnCount() == 1)
				return matrix.transpose(Calculation.Ret.LINK).mtimes(m)
						.getAsDouble(0, 0);
		return Double.NaN;
	}

	@Override
	public double dotColumnColumn(int col1, int col2) {
		// TEST OK
		final Matrix mcol1 = matrix.selectColumns(Calculation.Ret.LINK, col1);
		final Matrix mcol2 = matrix.selectColumns(Calculation.Ret.LINK, col2);
		return mcol1.mtimes(mcol2).getAsDouble(0, 0);
	}

	@Override
	public double dotColumnColumn(int col1, MatrixOKT mat, int col2) {
		// TEST OK
		final Matrix mcol1 = matrix.selectColumns(Calculation.Ret.LINK, col1);
		final Matrix mcol2 = ((UJMPMatrixOKT) mat).matrix.selectColumns(
				Calculation.Ret.LINK, col2);
		return mcol1.mtimes(mcol2).getAsDouble(0, 0);
	}

	@Override
	public double dotColumnRow(int column, int row) {
		// TEST OK
		final Matrix mcolumn = matrix.selectColumns(Calculation.Ret.LINK,
				column);
		final Matrix mrow = matrix.selectRows(Calculation.Ret.LINK, row);
		return mrow.mtimes(mcolumn).getAsDouble(0, 0);
	}

	@Override
	public double dotColumnRow(int column, MatrixOKT mat, int row) {
		// TEST OK
		final Matrix mcolumn = matrix.selectColumns(Calculation.Ret.LINK,
				column);
		final Matrix mrow = ((UJMPMatrixOKT) mat).matrix.selectRows(
				Calculation.Ret.LINK, row);
		return mrow.mtimes(mcolumn).getAsDouble(0, 0);
	}

	@Override
	public double dotRowColumn(int row, int column) {
		// TEST OK
		final Matrix mrow = matrix.selectRows(Calculation.Ret.LINK, row);
		final Matrix mcolumn = matrix.selectColumns(Calculation.Ret.LINK,
				column);
		return mrow.mtimes(mcolumn).getAsDouble(0, 0);
	}

	@Override
	public double dotRowColumn(int row, MatrixOKT mat, int column) {
		// TEST OK
		final Matrix mrow = matrix.selectRows(Calculation.Ret.LINK, row);
		final Matrix mcolumn = ((UJMPMatrixOKT) mat).matrix.selectColumns(
				Calculation.Ret.LINK, column);
		return mrow.mtimes(mcolumn).getAsDouble(0, 0);
	}

	@Override
	public double dotRowRow(int row1, int row2) {
		// TEST OK
		final Matrix mrow1 = matrix.selectRows(Calculation.Ret.LINK, row1);
		final Matrix mrow2 = matrix.selectRows(Calculation.Ret.LINK, row2);
		return mrow1.mtimes(mrow2).getAsDouble(0, 0);
	}

	@Override
	public double dotRowRow(int row1, MatrixOKT mat, int row2) {
		// TEST OK
		final Matrix mrow1 = matrix.selectRows(Calculation.Ret.LINK, row1);
		final Matrix mrow2 = ((UJMPMatrixOKT) mat).matrix.selectRows(
				Calculation.Ret.LINK, row2);
		return mrow1.mtimes(mrow2).getAsDouble(0, 0);
	}

	@Override
	@SuppressWarnings("unused")
	public Indices[] find(double target, int policy) {
		// TEST OK
		// trova tuti gli indndici assoluti
		// Matrix boolMat= solvePolicy(matrix, target, policy);
		/*
		 * Iterator<long[]> iterator= boolMat.nonZeroCoordinates().iterator();
		 * ArrayList<Number> rIndices= new ArrayList<Number>();
		 * ArrayList<Number> cIndices= new ArrayList<Number>();
		 */
		final int blockResize = 1000;
		// Matrix boolMat= solvePolicy(matrix, target, policy);
		final Iterable<long[]> iterable = matrix.nonZeroCoordinates();
		final Policy pol = getPolicy(target, policy);

		final long[] xx = new long[((Collection<?>) iterable).size()];
		final long[] yy = new long[((Collection<?>) iterable).size()];
		final Iterator<long[]> iterator = iterable.iterator();
		int count = 0;
		long[] xy;
		while (iterator.hasNext()) {
			xy = iterator.next();
			if (pol.isTruth(matrix.getAsDouble(xy[0], xy[1]))) {
				xx[count] = xy[0];
				yy[count++] = xy[1];
			}
		}
		/*
		 * for(int row= 0; row<rowsCount(); row++) { Matrix rowMat=
		 * solvePolicy(matrix.selectRows(Calculation.Ret.LINK, row), target,
		 * policy); Iterator<long[]> iterator=
		 * rowMat.nonZeroCoordinates().iterator(); while (iterator.hasNext()) {
		 * if( count>=xx.length ) { xx= Arrays.copyOf(xx,
		 * xx.length+blockResize); yy= Arrays.copyOf(yy, yy.length+blockResize);
		 * } xy= iterator.next(); xx[count]= xy[0]; yy[count++]= xy[1]; } }
		 */

		/*
		 * long[] xy; while (iterator.hasNext()) { xy= iterator.next(); if(
		 * boolMat.getAsBoolean(xy[0], xy[1]) ) { rIndices.add(xy[0]);
		 * cIndices.add(xy[1]); } }
		 */
		return new Indices[] { new UJMPIndices(Arrays.copyOf(xx, count)),
				new UJMPIndices(Arrays.copyOf(yy, count)) };
	}

	@Override
	public Indices findColumnIndices(double target, int row, int policy) {
		// TEST OK
		final Matrix dimension = solvePolicy(
				matrix.selectRows(Calculation.Ret.LINK, row), target, policy);
		return dimension == null ? new UJMPIndices(new int[0])
				: convertToIndices(dimension);
	}

	@Override
	public Indices findIndices(double target, int policy) {
		// TEST OK
		// trova tuti gli indndici assoluti
		final Matrix boolMat = solvePolicy(matrix, target, policy);
		final Iterator<long[]> iterator = boolMat.nonZeroCoordinates()
				.iterator();
		final ArrayList<Number> indices = new ArrayList<Number>();
		long[] xy;
		final int nc = (int) boolMat.getRowCount();
		while (iterator.hasNext()) {
			xy = iterator.next();
			if (boolMat.getAsBoolean(xy[0], xy[1]))
				indices.add(xy[0] + xy[1] * nc);
		}
		return new UJMPIndices(indices);
	}

	@Override
	public Indices findRowIndices(double target, int column, int policy) {
		// TEST OK
		final Matrix dimension = solvePolicy(
				matrix.selectColumns(Calculation.Ret.LINK, column), target,
				policy);
		return dimension == null ? new UJMPIndices(new int[0])
				: convertToIndices(dimension);
	}

	@Override
	public double get(int row, int column) {
		// TEST OK
		return matrix.getAsDouble(row, column);
	}

	@Override
	public MatrixOKT getCopy() {
		// TEST OK
		return new UJMPMatrixOKT(UJMPMatrixFactoryOKT.getMatrixClone(matrix),
				factory);
	}

	private Policy getPolicy(final double target, int policy) {
		// TEST OK
		if ((policy & GREATER) == GREATER) {
			if ((policy & EQUAL) == EQUAL)
				return value -> value >= target;
			else
				return value -> value > target;
		} else if ((policy & LESS) == LESS) {
			if ((policy & EQUAL) == EQUAL)
				return value -> value <= target;
			else
				return value -> value < target;
		} else if ((policy & EQUAL) == EQUAL)
			if ((policy & NOT) == NOT)
				return value -> value != target;
			else
				return value -> value == target;
		return null;
	}

	@Override
	public MatrixOKT inverse() {
		// TEST OK
		return new UJMPMatrixOKT(matrix.inv(), factory);
	}

	@Override
	public boolean isEmpty() {
		// TEST OK
		return matrix.isEmpty();
	}

	@Override
	public double max() {
		return matrix.getMaxValue();
	}

	@Override
	public double min() {
		return matrix.getMinValue();
	}

	@Override
	public MatrixOKT mul(double value) {
		// TEST OK
		return new UJMPMatrixOKT(matrix.times(value), factory);
	}

	@Override
	public MatrixOKT mul(MatrixOKT mat) {
		// TEST OK
		return new UJMPMatrixOKT(matrix.mtimes(((UJMPMatrixOKT) mat).matrix),
				factory);
	}

	@Override
	public MatrixOKT mulMe(double value) {
		// TEST OK
		matrix = matrix.mtimes(Calculation.Ret.ORIG, false, value);
		return this;
	}

	@Override
	public MatrixOKT mulMe(MatrixOKT mat) {
		// TEST OK
		matrix = matrix.mtimes(Calculation.Ret.ORIG, false,
				((UJMPMatrixOKT) mat).matrix);
		return this;
	}

	@Override
	public int nnz() {
		// TEST OK
		return ((Collection<?>) matrix.nonZeroCoordinates()).size();
	}

	@Override
	public void print(String name) {
		// TEST OK
		System.out.println("Matrix " + name + ":");
		for (int i = 0; i < matrix.getRowCount(); i++) {
			for (int j = 0; j < matrix.getColumnCount(); j++)
				System.out.print(matrix.getAsDouble(i, j) + "\t");
			System.out.println();
		}
	}

	@Override
	public void put(Indices indices, double value) {
		// TEST OK
		final int n = rowsCount();
		for (final int j : indices.toArray())
			matrix.setAsDouble(value, j % n, j / n);
	}

	@Override
	public MatrixOKT put(Indices rowIndices, Indices columnIndices,
			MatrixOKT mat) {
		// TEST OK
		final Matrix m = ((UJMPMatrixOKT) mat).matrix;
		Matrix ref = matrix.selectRows(Calculation.Ret.LINK,
				rowIndices.toLongArray());
		ref = ref.selectColumns(Calculation.Ret.LINK,
				columnIndices.toLongArray());
		for (int i = 0; i < rowIndices.count(); i++)
			for (int j = 0; j < columnIndices.count(); j++)
				ref.setAsDouble(m.getAsDouble(i, j), i, j);
		return this;
	}

	@Override
	public MatrixOKT putColumn(int column, MatrixOKT mat) {
		// TEST OK
		final Matrix m = ((UJMPMatrixOKT) mat).matrix;
		final Matrix ref = matrix.selectColumns(Calculation.Ret.LINK, column);
		if (m.isRowVector())
			for (int i = 0; i < m.getRowCount(); i++)
				ref.setAsDouble(m.getAsDouble(0, i), i, 0);
		else
			for (int j = 0; j < m.getColumnCount(); j++)
				ref.setAsDouble(m.getAsDouble(0, j), j, 0);
		return this;
	}

	@Override
	public MatrixOKT putRow(int row, MatrixOKT mat) {
		// TEST OK
		final Matrix m = ((UJMPMatrixOKT) mat).matrix;
		final Matrix ref = matrix.selectRows(Calculation.Ret.LINK, row);
		int i, j;
		try {
			if (m.isRowVector())
				for (i = 0; i < m.getRowCount(); i++)
					ref.setAsDouble(m.getAsDouble(i, 0), 0, i);
			else
				for (j = 0; j < m.getColumnCount(); j++)
					ref.setAsDouble(m.getAsDouble(0, j), 0, j);
			return this;
		} catch (final Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	@Override
	public MatrixOKT replace(double target, double value) {
		// TEST OK
		matrix.replace(Calculation.Ret.ORIG, target, value);
		return this;
	}

	@Override
	public MatrixOKT rows(boolean linked, Indices indices) {
		return new UJMPMatrixOKT(matrix.selectRows(Calculation.Ret.LINK,
				indices.toLongArray()), factory);
	}

	@Override
	public MatrixOKT rows(boolean linked, int... rows) {
		final long[] tmp = new long[rows.length];
		for (int j = 0; j < tmp.length; j++)
			tmp[j] = rows[j];
		return new UJMPMatrixOKT(matrix.selectRows(Calculation.Ret.LINK, tmp),
				factory);
	}

	@Override
	public MatrixOKT rows(Indices indices) {
		// TEST OK
		return new UJMPMatrixOKT(matrix.selectRows(Calculation.Ret.NEW,
				indices.toLongArray()), factory);
	}

	@Override
	public MatrixOKT rows(int... rows) {
		// TEST OK
		final long[] tmp = new long[rows.length];
		for (int j = 0; j < tmp.length; j++)
			tmp[j] = rows[j];
		final Matrix m = matrix.selectRows(Calculation.Ret.LINK, tmp);
		return new UJMPMatrixOKT(matrix.selectRows(Calculation.Ret.NEW, tmp),
				factory);
	}

	@Override
	public int rowsCount() {
		// TEST OK
		return (int) matrix.getRowCount();
	}

	@Override
	public void set(int row, int column, double value) {
		// TEST OK
		matrix.setAsDouble(value, row, column);
	}

	private Matrix solvePolicy(Matrix dimension, double target, int policy) {
		// TEST OK
		if ((policy & GREATER) == GREATER) {
			if ((policy & EQUAL) == EQUAL)
				return dimension.ge(Calculation.Ret.LINK, target);
			else
				return dimension.gt(Calculation.Ret.LINK, target);
		} else if ((policy & LESS) == LESS) {
			if ((policy & EQUAL) == EQUAL)
				return dimension.le(Calculation.Ret.LINK, target);
			else
				return dimension.lt(Calculation.Ret.LINK, target);
		} else if ((policy & EQUAL) == EQUAL)
			if ((policy & NOT) == NOT)
				return dimension.ne(Calculation.Ret.LINK, target);
			else
				return dimension.eq(Calculation.Ret.LINK, target);
		return null;
	}

	@Override
	public MatrixOKT sub(double value) {
		// TEST OK
		return new UJMPMatrixOKT(matrix.minus(value), factory);
	}

	@Override
	public MatrixOKT sub(MatrixOKT mat) {
		// TEST OK
		return new UJMPMatrixOKT(matrix.minus(((UJMPMatrixOKT) mat).matrix),
				factory);
	}

	@Override
	public MatrixOKT submatrix(Indices rowIndices, Indices columnIndices) {
		// TEST OK
		if (rowIndices.count() == 1) {
			final Matrix m = matrix.selectRows(Calculation.Ret.LINK,
					rowIndices.get(0)).selectColumns(Calculation.Ret.NEW,
					columnIndices.toLongArray());
			return new UJMPMatrixOKT(m, factory);
		} else if (columnIndices.count() == 1) {
			final Matrix m = matrix.selectColumns(Calculation.Ret.LINK,
					columnIndices.get(0)).selectRows(Calculation.Ret.NEW,
					rowIndices.toLongArray());
			return new UJMPMatrixOKT(m, factory);
		}
		return rows(rowIndices).columns(columnIndices);
	}

	@Override
	public MatrixOKT subMe(double value) {
		// TEST OK
		matrix = matrix.minus(Calculation.Ret.ORIG, false, value);
		return this;
	}

	@Override
	public MatrixOKT subMe(MatrixOKT mat) {
		// TEST OK
		matrix.minus(Calculation.Ret.ORIG, false, ((UJMPMatrixOKT) mat).matrix);
		return this;
	}

	@Override
	public MatrixOKT sum(double value) {
		// TEST OK
		return new UJMPMatrixOKT(matrix.plus(value), factory);
	}

	@Override
	public MatrixOKT sum(MatrixOKT mat) {
		// TEST OK
		return new UJMPMatrixOKT(matrix.plus(((UJMPMatrixOKT) mat).matrix),
				factory);
	}

	@Override
	public MatrixOKT sumMe(double value) {
		// TEST OK
		return new UJMPMatrixOKT(
				matrix.plus(Calculation.Ret.ORIG, false, value), factory);
	}

	@Override
	public MatrixOKT sumMe(MatrixOKT mat) {
		// TEST OK
		matrix = matrix.plus(Calculation.Ret.ORIG, false,
				((UJMPMatrixOKT) mat).matrix);
		return this;
	}

	@Override
	public MatrixOKT transpose() {
		// TEST OK
		return new UJMPMatrixOKT(matrix.transpose(), factory);
	}

	@Override
	public MatrixOKT transpose(boolean link) {
		// TEST OK
		return new UJMPMatrixOKT(matrix.transpose(Calculation.Ret.LINK),
				factory);
	}
}