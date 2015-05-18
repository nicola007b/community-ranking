package com.oktsrl.math.impls.ujmp;

import java.util.ArrayList;

import com.oktsrl.math.Indices;

public final class UJMPIndices implements Indices {
	protected long[] indices;

	public UJMPIndices(ArrayList<Number> indices) {
		this.indices = new long[indices.size()];
		int j = 0;
		for (final Number index : indices)
			this.indices[j++] = index.longValue();
	}

	public UJMPIndices(int[] indices) {
		this.indices = new long[indices.length];
		for (int j = 0; j < indices.length; j++)
			this.indices[j] = indices[j];
	}

	public UJMPIndices(long[] indices) {
		this.indices = indices;
	}

	@Override
	public int count() {
		return indices.length;
	}

	@Override
	public int get(int pos) {
		return (int) indices[pos];
	}

	@Override
	public Object getSource() {
		return indices;
	}

	@Override
	public boolean isEmpty() {
		return indices == null || indices.length == 0;
	}

	@Override
	public int[] toArray() {
		final int[] res = new int[indices.length];
		for (int j = 0; j < res.length; j++)
			res[j] = (int) indices[j];
		return res;
	}

	@Override
	public long[] toLongArray() {
		return indices;
	}

	@Override
	public String toString() {
		final StringBuilder builder = new StringBuilder("[");
		boolean isNotFirst = false;
		for (final Number index : indices) {
			if (isNotFirst)
				builder.append(",");
			else
				isNotFirst = true;
			builder.append(index);
		}
		builder.append("]");
		return builder.toString();
	}
}