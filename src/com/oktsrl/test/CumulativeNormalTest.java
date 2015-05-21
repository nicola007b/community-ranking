package com.oktsrl.test;

import org.apache.commons.math3.distribution.NormalDistribution;

public class CumulativeNormalTest {

	public static void main(String[] args) {
		final NormalDistribution nd = new NormalDistribution();

		for (int i = -10; i <= 10; ++i)
			System.out.println(nd.cumulativeProbability(i));
	}
}