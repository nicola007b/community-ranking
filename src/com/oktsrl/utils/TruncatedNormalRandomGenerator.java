package com.oktsrl.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.LinkedList;
import java.util.Random;
import java.util.StringTokenizer;

public class TruncatedNormalRandomGenerator {
	private static final double[] x;
	private static final double[] yu;
	private static final int[] ncell;

	private static final String xFile = ".dat/tnx.txt";
	private static final String yuFile = ".dat/tnyu.txt";
	private static final String ncellFile = ".dat/tnncell.txt";

	static {
		x = readDoubles(xFile);
		yu = readDoubles(yuFile);
		ncell = readIntegers(ncellFile);
	}

	private static double[] readDoubles(final String file) {
		try {
			final BufferedReader br = new BufferedReader(new FileReader(file));
			String line;
			final LinkedList<Double> list = new LinkedList<Double>();

			while ((line = br.readLine()) != null) {
				final StringTokenizer st = new StringTokenizer(line,
						" \t\n\r\f,", false);

				while (st.hasMoreTokens())
					list.add(Double.parseDouble(st.nextToken()));
			}

			br.close();

			final double[] ret = new double[list.size()];

			int i = 0;

			for (final Double val : list)
				ret[i++] = val;

			return ret;
		} catch (final Exception ex) {
			throw new RuntimeException(ex);
		}
	}

	private static int[] readIntegers(final String file) {
		try {
			final BufferedReader br = new BufferedReader(new FileReader(file));
			String line;
			final LinkedList<Integer> list = new LinkedList<Integer>();

			while ((line = br.readLine()) != null) {
				final StringTokenizer st = new StringTokenizer(line,
						" \t\n\r\f,", false);

				while (st.hasMoreTokens())
					list.add(Integer.parseInt(st.nextToken(), 10));
			}

			br.close();

			final int[] ret = new int[list.size()];

			int i = 0;

			for (final Integer val : list)
				ret[i++] = val;

			return ret;
		} catch (final Exception ex) {
			throw new RuntimeException(ex);
		}
	}

	private final Random r;

	public TruncatedNormalRandomGenerator(final long seed) {
		r = new Random(seed);
	}

	private double chopin(final double inf, final double sup, final double xmax) {
		// Design variables

		// if kb-ka < kmin then use a rejection algorithm
		final int kmin = 5;

		// 1/h, h being the minimal interval
		final double INVH = 1631.73284006;

		// range

		// = - floor(x(1)/h)
		final int I0 = 3271;

		// = log(2*pi)
		final double ALPHA = 1.837877066409345;

		// Index of the right tail
		final int N = 4001;

		// y_l of the leftmost rectangle
		final double yl0 = 0.053513975472;

		// y_l of the rightmost rectangle
		final double ylN = 0.000914116389555;

		// Compute kInf and kSup
		int i = I0 + (int) Math.floor(inf * INVH);
		final int kInf = ncell[i];
		int kSup;

		if (sup >= xmax)
			kSup = N;
		else {
			i = I0 + (int) Math.floor(sup * INVH);
			kSup = ncell[i];
		}

		// If |b - a| is small, use rejection algorithm with a truncated
		// exponential proposal
		if (Math.abs(kSup - kInf) < kmin) {
			boolean stop = false;
			final double twoInfSq = 2 * inf * inf;
			final double expInfSup = Math.exp(-inf * (sup - inf)) - 1;
			double z = 0;

			while (!stop) {
				z = Math.log(1 + r.nextDouble() * expInfSup);
				final double e = -Math.log(r.nextDouble());
				stop = twoInfSq * e > z * z;
			}

			return inf - z / inf;
		}

		while (true) {

			// Sample integer between ka and kb
			final int k = r.nextInt(kSup - kInf + 1) + kInf;

			if (k == N + 1) {

				// Right tail
				final double lbound = x[x.length - 1];
				double z = -Math.log(r.nextDouble());
				final double e = -Math.log(r.nextDouble());

				z = z / lbound;

				if (z * z <= 2 * e && z < sup - lbound)
					// Accept this proposition, otherwise reject
					return lbound + z;

			} else if (k <= kInf + 2 || k >= kSup && sup < xmax) {

				// Two leftmost and rightmost regions
				final double sim = x[k] + (x[k + 1] - x[k]) * r.nextDouble();

				if (sim >= inf && sim <= sup) {
					// Accept this proposition, otherwise reject
					final double simy = yu[k] * r.nextDouble();

					// Compute y_l from y_k
					double ylk;

					if (k == 1)
						ylk = yl0;
					else if (k == N)
						ylk = ylN;
					else if (k <= 1954)
						ylk = yu[k - 1];
					else
						ylk = yu[k + 1];

					if (simy < ylk
							|| sim * sim + 2 * Math.log(simy) + ALPHA < 0)
						return sim;
				}

			} else {

				// All the other boxes

				final double u = r.nextDouble();
				final double simy = yu[k] * u;
				final double d = x[k + 1] - x[k];

				// Compute y_l from y_k
				double ylk;

				if (k == 1)
					ylk = yl0;
				else if (k == N)
					ylk = ylN;
				else if (k <= 1954)
					ylk = yu[k - 1];
				else
					ylk = yu[k + 1];

				if (simy < ylk) // That's what happens most of the time
					return x[k] + u * d * yu[k] / ylk;

				final double sim = x[k] + d * r.nextDouble();

				// Otherwise, check you're below the pdf curve

				if (sim * sim + 2 * Math.log(simy) + ALPHA < 0)
					return sim;
			}
		}
	}

	public double next(final double inf, final double sup) {
		return next(inf, sup, 0, 1);
	}

	public double next(final double inf, final double sup, final double mu,
			final double sigma) {

		return rtstdnorm((inf - mu) / sigma, (sup - mu) / sigma) * sigma + mu;
	}

	/**
	 * Pseudorandom numbers from a truncated (normalized) Gaussian distribution
	 * (i.e. rtnorm(a,b,0,1))
	 * 
	 * @param inf
	 * @param sup
	 * @return
	 */
	private double rtstdnorm(final double inf, final double sup) {
		if (inf >= sup)
			throw new RuntimeException("\"sup\" must be greater than \"inf\"");

		if (Math.abs(inf) > Math.abs(sup))
			return -rtstdnorm(-sup, -inf);

		// Left and right limits
		final double xmin = -2.00443204036;
		final double xmax = 3.48672170399;

		// If a in the right tail (a > xmax), use rejection algorithm with a
		// truncated exponential proposal
		if (inf > xmax) {
			boolean stop = false;
			final double twoInfSq = 2 * inf * inf;
			final double expInfSup = Math.exp(-inf * (sup - inf)) - 1;
			double z = 0;

			while (!stop) {
				z = Math.log(1 + r.nextDouble() * expInfSup);
				final double e = -Math.log(r.nextDouble());
				stop = twoInfSq * e > z * z;
			}

			return inf - z / inf;
		}

		// If a in the left tail a < xmin, use rejection algorithm with a
		// Gaussian proposal
		if (inf < xmin) {
			boolean stop = false;
			double rand = (inf + sup) / 2;

			while (!stop) {
				rand = r.nextGaussian();
				stop = rand >= inf && rand <= sup;
			}

			return rand;
		}

		return chopin(inf, sup, xmax);
	}
}