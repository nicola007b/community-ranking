package com.oktsrl.math.impls;

import com.oktsrl.math.Indices;

public abstract class AbstractIndices implements Indices {
	protected int[] indices;

	protected AbstractIndices(int[] indices) {
		this.indices = indices;
	}

	@Override
	public int count() {
		return indices.length;
	}

	@Override
	public int get(int pos) {
		return indices[pos];
	}

	@Override
	public Object getSource() {
		return indices;
	}

	@Override
	public boolean isEmpty() {
		return indices.length == 0;
	}

	@Override
	public int[] toArray() {
		return indices;
	}

	@Override
	public long[] toLongArray() {
		return new long[0];
	}

	@Override
	public String toString() {
		final StringBuilder builder = new StringBuilder("[");
		boolean isNotFirst = false;
		for (final int index : indices) {
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