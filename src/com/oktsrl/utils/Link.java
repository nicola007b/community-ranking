package com.oktsrl.utils;

public final class Link {
	private int source;
	private int destination;

	public Link() {
		// Do nothing!
	}

	public Link(final int source, final int destination) {
		this.source = source;
		this.destination = destination;
	}

	public int getDestination() {
		return destination;
	}

	public int getSource() {
		return source;
	}

	public void setDestination(final int destination) {
		this.destination = destination;
	}

	public void setSource(final int source) {
		this.source = source;
	}

	@Override
	public String toString() {
		return "(" + source + "->" + destination + ")";
	}
}