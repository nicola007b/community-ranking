package com.oktsrl.utils;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

public class BidimensionalIndex implements Serializable {

	private static final long serialVersionUID = -598627181947933957L;

	private final IntArrayList rowIds;
	private final IntArrayList columnIds;
	private final Map<Integer, Integer> rowIndexes;
	private final Map<Integer, Integer> columnIndexes;

	public BidimensionalIndex() {
		rowIds = new IntArrayList(1024);
		columnIds = new IntArrayList(1024);
		rowIndexes = new HashMap<Integer, Integer>(1024);
		columnIndexes = new HashMap<Integer, Integer>(1024);
	}

	public int addColumn(final int columnId) {
		final Integer i = columnIndexes.get(columnId);

		if (i == null) {
			final int size = columnIds.size();
			columnIndexes.put(columnId, size);
			columnIds.add(columnId);
			return size;
		}

		return i;
	}

	public int addRow(final int rowId) {
		final Integer i = rowIndexes.get(rowId);

		if (i == null) {
			final int size = rowIds.size();
			rowIndexes.put(rowId, size);
			rowIds.add(rowId);
			return size;
		}

		return i;
	}

	public int columnId(final int columnIndex) {
		return columnIds.get(columnIndex);
	}

	public int columnIndex(final int columnId) {
		return columnIndexes.get(columnId);
	}

	@Override
	protected void finalize() throws Throwable {
		if (rowIndexes != null)
			rowIndexes.clear();

		if (columnIndexes != null)
			columnIndexes.clear();
	}

	public boolean hasColumnId(final int columnId) {
		return columnIndexes.containsKey(columnId);
	}

	public boolean hasColumnIndex(final int columnIndex) {
		return columnIndex > 0 && columnIndex < rowIds.size();
	}

	public boolean hasRowId(final int rowId) {
		return rowIndexes.containsKey(rowId);
	}

	public boolean hasRowIndex(final int rowIndex) {
		return rowIndex > 0 && rowIndex < rowIds.size();
	}

	public int nColumns() {
		return columnIds.size();
	}

	public int nRows() {
		return rowIds.size();
	}

	public int rowId(final int rowIndex) {
		return rowIds.get(rowIndex);
	}

	public int rowIndex(final int rowId) {
		return rowIndexes.get(rowId);
	}
}