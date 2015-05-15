package com.oktsrl.utils;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;

public class IntArrayList implements Serializable {

	private static final long serialVersionUID = 178637516575007260L;

	/**
	 * The array buffer into which the elements of the ArrayList are stored. The
	 * capacity of the ArrayList is the length of this array buffer.
	 */
	private int[] elementData;

	/**
	 * The size of the ArrayList (the number of elements it contains).
	 * 
	 * @serial
	 */
	private int size;

	/**
	 * Constructs an empty list with an initial capacity of ten.
	 */
	public IntArrayList() {
		this(10);
	}

	public IntArrayList(final Collection<Integer> c) {
		elementData = new int[c.size()];
		size = elementData.length;

		int j = 0;

		for (final int i : c)
			elementData[j++] = i;
	}

	/**
	 * Constructs an empty list with the specified initial capacity.
	 * 
	 * @param initialCapacity
	 *            the initial capacity of the list
	 * @exception IllegalArgumentException
	 *                if the specified initial capacity is negative
	 */
	public IntArrayList(final int initialCapacity) {
		if (initialCapacity < 0)
			throw new IllegalArgumentException("Illegal Capacity: "
					+ initialCapacity);

		elementData = new int[initialCapacity];
	}

	/**
	 * Constructs a list containing the elements of the specified collection, in
	 * the order they are returned by the collection's iterator.
	 * 
	 * @param c
	 *            the collection whose elements are to be placed into this list
	 * @throws NullPointerException
	 *             if the specified collection is null
	 */
	public IntArrayList(final IntArrayList c) {
		elementData = new int[c.size];
		size = elementData.length;

		System.arraycopy(c.elementData, 0, elementData, 0, size);
	}

	/**
	 * Appends the specified element to the end of this list.
	 * 
	 * @param e
	 *            element to be appended to this list
	 * @return <tt>true</tt> (as specified by {@link Collection#add})
	 */
	public boolean add(final int i) {
		ensureCapacity(size + 1); // Increments modCount!!
		elementData[size++] = i;
		return true;
	}

	/**
	 * Inserts the specified element at the specified position in this list.
	 * Shifts the element currently at that position (if any) and any subsequent
	 * elements to the right (adds one to their indices).
	 * 
	 * @param index
	 *            index at which the specified element is to be inserted
	 * @param element
	 *            element to be inserted
	 * @throws IndexOutOfBoundsException
	 *             {@inheritDoc}
	 */
	public void add(final int index, final int element) {
		if (index > size || index < 0)
			throw new IndexOutOfBoundsException("Index: " + index + ", Size: "
					+ size);

		ensureCapacity(size + 1); // Increments modCount!!
		System.arraycopy(elementData, index, elementData, index + 1, size
				- index);
		elementData[index] = element;
		size++;
	}

	public boolean addAll(final int index, final Collection<Integer> c) {
		if (index > size || index < 0)
			throw new IndexOutOfBoundsException("Index: " + index + ", Size: "
					+ size);

		final int numNew = c.size();
		ensureCapacity(size + numNew); // Increments modCount

		final int numMoved = size - index;
		if (numMoved > 0)
			System.arraycopy(elementData, index, elementData, index + numNew,
					numMoved);

		int j = index;

		for (final int i : c)
			elementData[j++] = i;

		size += numNew;
		return numNew != 0;
	}

	/**
	 * Inserts all of the elements in the specified collection into this list,
	 * starting at the specified position. Shifts the element currently at that
	 * position (if any) and any subsequent elements to the right (increases
	 * their indices). The new elements will appear in the list in the order
	 * that they are returned by the specified collection's iterator.
	 * 
	 * @param index
	 *            index at which to insert the first element from the specified
	 *            collection
	 * @param c
	 *            collection containing elements to be added to this list
	 * @return <tt>true</tt> if this list changed as a result of the call
	 * @throws IndexOutOfBoundsException
	 *             {@inheritDoc}
	 * @throws NullPointerException
	 *             if the specified collection is null
	 */
	public boolean addAll(final int index, final IntArrayList c) {
		if (index > size || index < 0)
			throw new IndexOutOfBoundsException("Index: " + index + ", Size: "
					+ size);

		final int[] a = c.elementData;
		final int numNew = a.length;
		ensureCapacity(size + numNew); // Increments modCount

		final int numMoved = size - index;
		if (numMoved > 0)
			System.arraycopy(elementData, index, elementData, index + numNew,
					numMoved);

		System.arraycopy(a, 0, elementData, index, numNew);
		size += numNew;

		return numNew != 0;
	}

	/**
	 * Appends all of the elements in the specified collection to the end of
	 * this list, in the order that they are returned by the specified
	 * collection's Iterator. The behavior of this operation is undefined if the
	 * specified collection is modified while the operation is in progress.
	 * (This implies that the behavior of this call is undefined if the
	 * specified collection is this list, and this list is nonempty.)
	 * 
	 * @param c
	 *            collection containing elements to be added to this list
	 * @return <tt>true</tt> if this list changed as a result of the call
	 * @throws NullPointerException
	 *             if the specified collection is null
	 */
	public boolean addAll(final IntArrayList c) {
		final int[] a = c.elementData;
		final int numNew = a.length;
		ensureCapacity(size + numNew); // Increments modCount
		System.arraycopy(a, 0, elementData, size, numNew);
		size += numNew;
		return numNew != 0;
	}

	/**
	 * Removes all of the elements from this list. The list will be empty after
	 * this call returns.
	 */
	public void clear() {
		size = 0;
	}

	/**
	 * Returns <tt>true</tt> if this list contains the specified element. More
	 * formally, returns <tt>true</tt> if and only if this list contains at
	 * least one element <tt>e</tt> such that
	 * <tt>(o==null&nbsp;?&nbsp;e==null&nbsp;:&nbsp;o.equals(e))</tt>.
	 * 
	 * @param o
	 *            element whose presence in this list is to be tested
	 * @return <tt>true</tt> if this list contains the specified element
	 */
	public boolean contains(final int i) {
		return indexOf(i) >= 0;
	}

	/**
	 * Increases the capacity of this <tt>ArrayList</tt> instance, if necessary,
	 * to ensure that it can hold at least the number of elements specified by
	 * the minimum capacity argument.
	 * 
	 * @param minCapacity
	 *            the desired minimum capacity
	 */
	public void ensureCapacity(final int minCapacity) {
		final int oldCapacity = elementData.length;

		if (minCapacity > oldCapacity) {
			int newCapacity = oldCapacity * 3 / 2 + 1;

			if (newCapacity < minCapacity)
				newCapacity = minCapacity;

			elementData = Arrays.copyOf(elementData, newCapacity);
		}
	}

	/*
	 * Private remove method that skips bounds checking and does not return the
	 * value removed.
	 */
	private void fastRemove(final int index) {
		final int numMoved = size - index - 1;

		if (numMoved > 0)
			System.arraycopy(elementData, index + 1, elementData, index,
					numMoved);

		--size;
	}

	/**
	 * Returns the element at the specified position in this list.
	 * 
	 * @param index
	 *            index of the element to return
	 * @return the element at the specified position in this list
	 * @throws IndexOutOfBoundsException
	 *             {@inheritDoc}
	 */
	public int get(final int index) {
		RangeCheck(index);

		return elementData[index];
	}

	/**
	 * Returns the index of the first occurrence of the specified element in
	 * this list, or -1 if this list does not contain the element. More
	 * formally, returns the lowest index <tt>i</tt> such that
	 * <tt>(o==null&nbsp;?&nbsp;get(i)==null&nbsp;:&nbsp;o.equals(get(i)))</tt>,
	 * or -1 if there is no such index.
	 */
	public int indexOf(final int element) {
		for (int i = 0; i < size; i++)
			if (element == elementData[i])
				return i;

		return -1;
	}

	/**
	 * Returns <tt>true</tt> if this list contains no elements.
	 * 
	 * @return <tt>true</tt> if this list contains no elements
	 */
	public boolean isEmpty() {
		return size == 0;
	}

	/**
	 * Returns the index of the last occurrence of the specified element in this
	 * list, or -1 if this list does not contain the element. More formally,
	 * returns the highest index <tt>i</tt> such that
	 * <tt>(o==null&nbsp;?&nbsp;get(i)==null&nbsp;:&nbsp;o.equals(get(i)))</tt>,
	 * or -1 if there is no such index.
	 */
	public int lastIndexOf(final int element) {
		for (int i = size - 1; i >= 0; i--)
			if (element == elementData[i])
				return i;

		return -1;
	}

	/**
	 * Checks if the given index is in range. If not, throws an appropriate
	 * runtime exception. This method does *not* check if the index is negative:
	 * It is always used immediately prior to an array access, which throws an
	 * ArrayIndexOutOfBoundsException if index is negative.
	 */
	private void RangeCheck(final int index) {
		if (index >= size)
			throw new IndexOutOfBoundsException("Index: " + index + ", Size: "
					+ size);
	}

	/**
	 * Removes the first occurrence of the specified element from this list, if
	 * it is present. If the list does not contain the element, it is unchanged.
	 * More formally, removes the element with the lowest index <tt>i</tt> such
	 * that
	 * <tt>(o==null&nbsp;?&nbsp;get(i)==null&nbsp;:&nbsp;o.equals(get(i)))</tt>
	 * (if such an element exists). Returns <tt>true</tt> if this list contained
	 * the specified element (or equivalently, if this list changed as a result
	 * of the call).
	 * 
	 * @param o
	 *            element to be removed from this list, if present
	 * @return <tt>true</tt> if this list contained the specified element
	 */
	public boolean remove(final int element) {
		for (int index = 0; index < size; index++)
			if (element == elementData[index]) {
				fastRemove(index);
				return true;
			}

		return false;
	}

	/**
	 * Removes the element at the specified position in this list. Shifts any
	 * subsequent elements to the left (subtracts one from their indices).
	 * 
	 * @param index
	 *            the index of the element to be removed
	 * @return the element that was removed from the list
	 * @throws IndexOutOfBoundsException
	 *             {@inheritDoc}
	 */
	public int removeByIndex(final int index) {
		RangeCheck(index);

		final int oldValue = elementData[index];
		final int numMoved = size - index - 1;
		--size;

		if (numMoved > 0)
			System.arraycopy(elementData, index + 1, elementData, index,
					numMoved);

		return oldValue;
	}

	/**
	 * Removes from this list all of the elements whose index is between
	 * <tt>fromIndex</tt>, inclusive, and <tt>toIndex</tt>, exclusive. Shifts
	 * any succeeding elements to the left (reduces their index). This call
	 * shortens the list by <tt>(toIndex - fromIndex)</tt> elements. (If
	 * <tt>toIndex==fromIndex</tt>, this operation has no effect.)
	 * 
	 * @param fromIndex
	 *            index of first element to be removed
	 * @param toIndex
	 *            index after last element to be removed
	 * @throws IndexOutOfBoundsException
	 *             if fromIndex or toIndex out of range (fromIndex &lt; 0 ||
	 *             fromIndex &gt;= size() || toIndex &gt; size() || toIndex &lt;
	 *             fromIndex)
	 */
	protected void removeRange(final int fromIndex, final int toIndex) {
		System.arraycopy(elementData, toIndex, elementData, fromIndex, size
				- toIndex);
		size = size - (toIndex - fromIndex);
	}

	// FIXME da controllare
	public IntArrayList retainAll(final IntArrayList other) {
		final int[] tmp1 = new int[elementData.length];
		final int[] tmp2 = new int[other.elementData.length];

		System.arraycopy(elementData, 0, tmp1, 0, elementData.length);
		System.arraycopy(other.elementData, 0, tmp2, 0,
				other.elementData.length);

		int[] min;
		int[] max;

		if (tmp1.length < tmp2.length) {
			min = tmp1;
			max = tmp2;
		} else {
			min = tmp2;
			max = tmp1;
		}

		Arrays.sort(min);
		Arrays.sort(max);

		final LinkedList<Integer> list = new LinkedList<Integer>();

		for (final int i : min)
			if (Arrays.binarySearch(max, i) >= 0)
				list.add(i);

		final int[] data = new int[list.size()];

		final Iterator<Integer> it = list.iterator();
		int i = 0;

		while (it.hasNext())
			data[i++] = it.next();

		final IntArrayList ret = new IntArrayList();
		ret.elementData = data;
		ret.size = data.length;

		return ret;
	}

	/**
	 * Replaces the element at the specified position in this list with the
	 * specified element.
	 * 
	 * @param index
	 *            index of the element to replace
	 * @param element
	 *            element to be stored at the specified position
	 * @return the element previously at the specified position
	 * @throws IndexOutOfBoundsException
	 *             {@inheritDoc}
	 */
	public int set(final int index, final int element) {
		RangeCheck(index);

		final int oldValue = elementData[index];
		elementData[index] = element;
		return oldValue;
	}

	/**
	 * Randomly permute the list using the specified source of randomness. All
	 * permutations occur with equal likelihood assuming that the source of
	 * randomness is fair.
	 * 
	 * @param rnd
	 *            the source of randomness to use to shuffle the list.
	 */
	public void shuffle(final Random rnd) {
		for (int i = size; i > 1; --i) {
			final int k = i - 1;
			final int j = rnd.nextInt(i);

			final int tmp = elementData[k];
			elementData[k] = elementData[j];
			elementData[j] = tmp;
		}
	}

	/**
	 * Returns the number of elements in this list.
	 * 
	 * @return the number of elements in this list
	 */
	public int size() {
		return size;
	}

	/**
	 * Sorts the list of ints into ascending numerical order.
	 */
	public void sort() {
		Arrays.sort(elementData);
	}

	/**
	 * Returns an array containing all of the elements in this list in proper
	 * sequence (from first to last element).
	 * 
	 * <p>
	 * The returned array will be "safe" in that no references to it are
	 * maintained by this list. (In other words, this method must allocate a new
	 * array). The caller is thus free to modify the returned array.
	 * 
	 * <p>
	 * This method acts as bridge between array-based and collection-based APIs.
	 * 
	 * @return an array containing all of the elements in this list in proper
	 *         sequence
	 */
	public int[] toArray() {
		return Arrays.copyOf(elementData, size);
	}

	@Override
	public String toString() {
		return Arrays.toString(elementData);
	}

	/**
	 * Trims the capacity of this <tt>ArrayList</tt> instance to be the list's
	 * current size. An application can use this operation to minimize the
	 * storage of an <tt>ArrayList</tt> instance.
	 */
	public void trimToSize() {
		if (size < elementData.length)
			elementData = Arrays.copyOf(elementData, size);
	}
}