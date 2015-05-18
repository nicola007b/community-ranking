package com.oktsrl.math;

public interface Indices {

	public int count();

	public int get(int pos);

	public Object getSource();

	public boolean isEmpty();

	public int[] toArray();

	public long[] toLongArray();
}