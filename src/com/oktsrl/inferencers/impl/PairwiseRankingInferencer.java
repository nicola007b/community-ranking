package com.oktsrl.inferencers.impl;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.distribution.NormalDistribution;

import com.oktsrl.Model;
import com.oktsrl.inferencers.BayesianInferencer;
import com.oktsrl.math.BuildMatrixFactoryOKT;
import com.oktsrl.math.Indices;
import com.oktsrl.math.MatrixFactoryOKT;
import com.oktsrl.math.MatrixOKT;
import com.oktsrl.math.Summation;
import com.oktsrl.models.impl.PairwiseRankingModel;
import com.oktsrl.utils.BidimensionalIndex;
import com.oktsrl.utils.Settings;
import com.oktsrl.utils.TruncatedNormalRandomGenerator;

public class PairwiseRankingInferencer implements BayesianInferencer {

	private static final int BINARY_SEARCH_ACTIVATION = 12;

	private static int binarySearch(int[] v, int key) {
		int low = 0;
		int high = v.length - 1;

		while (low <= high) {
			final int mid = low + high >>> 1;

			if (v[mid] < key)
				low = mid + 1;
			else if (v[mid] > key)
				high = mid - 1;
			else
				return mid; // key found
		}
		return -(low + 1); // key not found
	}

	private static void log(final String message, final Object... params) {
		System.out.println(String.format(message, params));
	}

	private static int min(int i, int j) {
		if (i < j)
			return i;

		return j;
	}

	private static String padd(double d, int i) {
		final long pow = Math.round(Math.pow(10, i));

		d *= pow;
		d = (double) Math.round(d) / pow;

		return String.valueOf(d);
	}

	private static void shuffle(int[] v, int elements, Random rnd) {
		for (int i = elements; i > 1; --i) {
			final int j = rnd.nextInt(i);
			final int tmp = v[i - i];
			v[i - 1] = v[j];
			v[j] = tmp;
		}
	}

	private static int sortedLinearSearch(int[] v, int key) {
		final int n = v.length;

		for (int i = 0; i < n; ++i) {
			if (v[i] > key)
				return -(i + 1);

			if (v[i] == key)
				return i;
		}

		return -(n + 1);
	}

	private static int sortedSearch(int[] v, int key) {
		if (v.length < BINARY_SEARCH_ACTIVATION)
			return sortedLinearSearch(v, key);
		else
			return binarySearch(v, key);
	}

	private static String timeToString(double time) {
		String measure = "secondi";
		if (time > 60d) {
			measure = "minuti";
			time /= 60d;
		} else if (time > 3600d) {
			measure = "ore";
			time /= 3600d;
		}
		return String.format("%.2f %s", time, measure);
	}

	protected int nUsers;
	protected int nItems;

	protected MatrixFactoryOKT factory;
	protected MatrixOKT socialNetwork;
	protected MatrixOKT preferenceMatrix;
	protected MatrixOKT unknownLinks;

	protected BidimensionalIndex index;

	protected MatrixOKT Theta;
	private MatrixOKT Omega;
	protected MatrixOKT[] ThetaAll;
	protected MatrixOKT[] OmegaAll;

	protected MatrixOKT Ze;
	protected MatrixOKT[] Zr;

	protected TruncatedNormalRandomGenerator tnrg;
	protected Random rg;

	// parametri configurabili
	protected int nFactors; // topics
	protected double alpha; // observation noise (precision)
	protected double beta; // observation noise (precision)
	protected double gamma1; // prior number of active elements
	protected double gamma2; // prior number of non-active elements
	protected int nTrials; // cicli di Gibbs
	protected int maxEpoch;
	protected int epochHistorySize; // numero di epoch da prendere in
	// considerazione
	protected double pairPercent; // rapporto di coppie da prendere in
	// considerazione per i confronti
	protected int minComparisons; // minimo numero di confronti di items per
	// item e per utente
	protected long startInferenceTime;

	// Hyper parameters
	private MatrixOKT muTheta, muOmega;
	private MatrixOKT sigmaThetaInv, sigmaOmegaInv;
	private MatrixOKT W0_theta, W0_omega;
	private MatrixOKT mu0_theta, mu0_omega;
	private int beta0_theta, beta0_omega;
	private int ni0_theta, ni0_omega;

	// support
	private HashMap<Integer, HashSet<Integer>>[] pairwiseComparisons;
	private HashMap<Integer, int[]> neighborsPerUser;
	private HashMap<Integer, int[]> usersPerItem;
	private HashMap<Integer, LinkedList<Integer>> Y;

	public PairwiseRankingInferencer(final Settings settings) {
		nFactors = settings.getInt("nTopics", 5);
		alpha = settings.getReal("alpha", 2d);
		beta = settings.getReal("beta", 2d);
		gamma1 = settings.getReal("gamma1", 2d);
		gamma2 = settings.getReal("gamma2", 2d);
		nTrials = settings.getInt("nTrials", 2);
		maxEpoch = settings.getInt("maxEpoch", 20);
		epochHistorySize = settings.getInt("epochHistorySize", 10);
		pairPercent = settings.getReal("pairPercent", 0.1);
		minComparisons = settings.getInt("minComparisons", 5);

		final String xFile = settings.getString("xFile");
		final String yuFile = settings.getString("yuFile");
		final String ncellFile = settings.getString("ncellFile");

		final int seed = settings.getInt("seed", 101);

		tnrg = new TruncatedNormalRandomGenerator(xFile, yuFile, ncellFile,
				seed);
		rg = new Random(seed);
	}

	@SuppressWarnings("all")
	private double computeAllInOneLogLikelihood_Alternative(int epoch) {
		final NormalDistribution nd = new NormalDistribution();

		final int limit = min(epoch, ThetaAll.length);
		final HashMap<Integer, HashMap<Integer, Double>>[] z_uvk = new HashMap[limit];

		System.out.println("Computing log-likelihood...");
		double oldPercent = 0;

		double log_Pr_ep = 0;
		double log_Pr_rp = 0;

		// Observed links & items
		for (int u = 0; u < nUsers; ++u) {
			final double percent = (double) u / nUsers;

			if (percent >= oldPercent + 0.1) {
				System.out.println("\tlog-likelihood: " + padd(percent * 50, 1)
						+ "%");
				oldPercent = percent;
			}

			// Links
			final int[] neighbors = neighborsPerUser.get(u);

			for (int l = -sortedSearch(neighbors, u) - 1, n = neighbors.length; l < n; ++l) {
				final int v = neighbors[l];
				double cumulative_z_uv = 0;

				for (int k = 0; k < limit; ++k) {
					final MatrixOKT thetaUT = ThetaAll[k].rows(u);
					final MatrixOKT thetaVT = OmegaAll[k].rows(v);

					// Anche senza trasposta il dot() funziona lo stesso
					final double phi = nd.cumulativeProbability(thetaVT
							.dot(thetaUT));

					cumulative_z_uv += phi;
				}

				log_Pr_ep += Math.log(cumulative_z_uv) - Math.log(limit);
			}

			// Items
			for (final Map.Entry<Integer, HashSet<Integer>> e : pairwiseComparisons[u]
					.entrySet()) {

				final int itemI = e.getKey();

				for (final int itemJ : e.getValue()) {
					double cumulative_w_uij = 0;

					for (int k = 0; k < limit; ++k) {
						final MatrixOKT omegaI = OmegaAll[k].rows(itemI);
						final MatrixOKT deltaOmegaT = omegaI.sub(OmegaAll[k]
								.rows(itemJ));
						final MatrixOKT thetaUT = ThetaAll[k].rows(u);

						// Anche senza trasposta il dot() funziona lo stesso
						final double phi = nd.cumulativeProbability(deltaOmegaT
								.dot(thetaUT));

						cumulative_w_uij += phi;
					}

					log_Pr_rp += 2 * (Math.log(cumulative_w_uij) - Math
							.log(limit));
				}
			}
		}

		// Negative links: we are assuming that all the selected unknown
		// links are negative XXX richiedi conferma a Beppe
		final Indices[] unknownIndices = unknownLinks
				.find(0, MatrixOKT.GREATER);
		final int n = unknownIndices[0].count();

		oldPercent = 0;

		for (int j = 0; j < n; j++) {
			final double percent = (double) j / n;

			if (percent >= oldPercent + 0.1) {
				System.out.println("\tlog-likelihood: "
						+ padd(percent * 50 + 50, 1) + "%");
				oldPercent = percent;
			}

			final int u = unknownIndices[0].get(j);
			final int v = unknownIndices[1].get(j);

			if (v <= u)
				continue;

			double cumulative_z_uv = 0;

			for (int k = 0; k < limit; ++k) {
				final MatrixOKT thetaUT = ThetaAll[k].rows(u);
				final MatrixOKT thetaVT = OmegaAll[k].rows(v);

				// Anche senza trasposta il dot() funziona lo stesso
				final double phi = nd.cumulativeProbability(thetaVT
						.dot(thetaUT));

				cumulative_z_uv += phi;
			}

			log_Pr_ep += Math.log(1 - cumulative_z_uv / limit);
		}

		System.out.println("\tlog-likelihood: 100.0% complete");

		return log_Pr_ep + log_Pr_rp;
	}

	private double computeIncrementalLogLikelihood(int epoch) {
		final NormalDistribution nd = new NormalDistribution();

		final MatrixOKT theta_k = ThetaAll[epoch];
		final MatrixOKT omega_k = OmegaAll[epoch];

		double log_Pr_ep = 0;
		double log_Pr_rp = 0;

		System.out.println("Computing log-likelihood...");
		double oldPercent = 0;

		// Observed links & items
		for (int u = 0; u < nUsers; ++u) {
			final double percent = (double) u / nUsers;

			if (percent >= oldPercent + 0.1) {
				System.out.println("\tlog-likelihood: " + padd(percent * 50, 1)
						+ "%");
				oldPercent = percent;
			}

			final MatrixOKT thetaUT = theta_k.rows(u);

			// Links
			final int[] neighbors = neighborsPerUser.get(u);

			for (int l = -sortedSearch(neighbors, u) - 1, n = neighbors.length; l < n; ++l) {
				final int v = neighbors[l];

				final MatrixOKT thetaVT = theta_k.rows(v);

				// Anche senza trasposta il dot() funziona lo stesso
				final double z_uv = thetaVT.dot(thetaUT);
				final double phi = nd.cumulativeProbability(z_uv);

				log_Pr_ep += Math.log(phi);
			}

			// Items
			for (final Map.Entry<Integer, HashSet<Integer>> e : pairwiseComparisons[u]
					.entrySet()) {

				final int itemI = e.getKey();
				final MatrixOKT omegaI = omega_k.rows(itemI);

				for (final int itemJ : e.getValue()) {
					final MatrixOKT deltaOmegaT = omegaI.sub(omega_k
							.rows(itemJ));

					// Anche senza trasposta il dot() funziona lo stesso
					final double w_uij = deltaOmegaT.dot(thetaUT);
					final double phi = nd.cumulativeProbability(w_uij);

					log_Pr_rp += 2 * Math.log(phi);
				}
			}
		}

		oldPercent = 0;
		int iter = 0;
		final int maxIter = Y.size();

		// Negative links
		for (final Map.Entry<Integer, LinkedList<Integer>> e : Y.entrySet()) {
			final double percent = (double) iter / maxIter;
			++iter;

			if (percent >= oldPercent + 0.1) {
				System.out.println("\tlog-likelihood: "
						+ padd(percent * 50 + 50, 1) + "%");
				oldPercent = percent;
			}

			final int u = e.getKey();
			final MatrixOKT thetaUT = theta_k.rows(u);

			for (final int v : e.getValue()) {
				if (v <= u)
					continue;

				final MatrixOKT thetaVT = theta_k.rows(v);

				// Anche senza trasposta il dot() funziona lo stesso
				final double z_uv = thetaVT.dot(thetaUT);
				final double phi = nd.cumulativeProbability(z_uv);

				log_Pr_ep += Math.log(1 - phi);
			}
		}

		System.out.println("\tlog-likelihood: 100.0% complete");

		return log_Pr_ep + log_Pr_rp;
	}

	private void hyperSample() {
		MatrixOKT[] hyperParams;

		hyperParams = samplingHyperParams(Theta, W0_theta, mu0_theta,
				beta0_theta, ni0_theta);

		muTheta = hyperParams[0];
		sigmaThetaInv = hyperParams[1].inverse();

		hyperParams = samplingHyperParams(Omega, W0_omega, mu0_omega,
				beta0_omega, ni0_omega);

		muOmega = hyperParams[0];
		sigmaOmegaInv = hyperParams[1].inverse();
	}

	private void initializeHyperParams() {

		// parameters of Inv-Whishart distribution for Theta
		W0_theta = factory.eye(nFactors);
		beta0_theta = 2;
		ni0_theta = nFactors;

		mu0_theta = factory.zeros(nFactors, 1);

		// parameters of Inv-Whishart distribution for Omega
		W0_omega = factory.eye(nFactors);
		beta0_omega = 2;
		ni0_omega = nFactors;

		mu0_omega = factory.zeros(nFactors, 1);
	}

	private void initializeParams() {
		Theta = factory.randn(nUsers, nFactors).mulMe(0.1);
		Omega = factory.randn(nItems, nFactors).mulMe(0.1);
		Y = new HashMap<Integer, LinkedList<Integer>>(nUsers);

		ThetaAll = new MatrixOKT[epochHistorySize];
		OmegaAll = new MatrixOKT[epochHistorySize];

		Ze = factory.sparse(nUsers, nUsers);
		Zr = new MatrixOKT[nUsers];

		for (int u = 0; u < nUsers; ++u)
			Zr[u] = factory.sparse(nItems, nItems);
	}

	@SuppressWarnings("unchecked")
	private void initializeSupportParams() {
		System.out.println("\tLoading item pairs for comparisons...");

		pairwiseComparisons = new HashMap[nUsers];

		final int[] copy = new int[nItems];

		double oldPercent = 0;

		for (int u = 0; u < nUsers; ++u) {
			final double percent = (double) u / nUsers;

			if (percent >= oldPercent + 0.1) {
				System.out.println("\t\tLoading: " + padd(percent * 100, 1)
						+ "%");
				oldPercent = percent;
			}

			final int[] items = preferenceMatrix.findColumnIndices(0, u,
					MatrixOKT.GREATER).toArray();

			pairwiseComparisons[u] = new HashMap<Integer, HashSet<Integer>>(
					items.length);

			System.arraycopy(items, 0, copy, 0, items.length);

			final int lowerbound = min(items.length, minComparisons);
			int n = (int) Math.ceil(items.length * pairPercent);

			if (n < lowerbound)
				n = lowerbound;

			for (final int i : items) {
				if (n > lowerbound)
					shuffle(copy, items.length, rg);

				for (int index = 0; index < n; ++index) {
					if (i == copy[index])
						continue;

					int min = i;
					int max = copy[index];

					if (min > max) {
						min = max;
						max = i;
					}

					HashSet<Integer> comparisons = pairwiseComparisons[u]
							.get(min);

					if (comparisons == null) {
						comparisons = new HashSet<Integer>(n);
						pairwiseComparisons[u].put(min, comparisons);
					}

					comparisons.add(max);
				}
			}
		}

		System.out.println("\t\tLoading: 100.0% complete");

		System.out.println("\tLoading neighbors per user...");

		neighborsPerUser = new HashMap<Integer, int[]>(nUsers);
		oldPercent = 0;

		for (int u = 0; u < nUsers; ++u) {
			final double percent = (double) u / nUsers;

			if (percent >= oldPercent + 0.1) {
				System.out.println("\t\tLoading: " + padd(percent * 100, 1)
						+ "%");
				oldPercent = percent;
			}

			neighborsPerUser.put(u,
					socialNetwork.findColumnIndices(0, u, MatrixOKT.NOT)
					.toArray());
		}

		System.out.println("\t\tLoading: 100.0% complete");

		System.out.println("\tLoading users per item...");

		usersPerItem = new HashMap<Integer, int[]>(nUsers);
		oldPercent = 0;

		for (int i = 0; i < nItems; ++i) {
			final double percent = (double) i / nItems;

			if (percent >= oldPercent + 0.1) {
				System.out.println("\t\tLoading: " + padd(percent * 100, 1)
						+ "%");
				oldPercent = percent;
			}

			usersPerItem.put(i,
					preferenceMatrix.findRowIndices(0, i, MatrixOKT.GREATER)
					.toArray());
		}

		System.out.println("\t\tLoading: 100.0% complete");
	}

	private void logNextStep(final int epoch) {
		if (epoch == 0)
			log("Step: %d di %d, tempo rimanente %s", epoch + 1, maxEpoch,
					"in valutazione");
		else {
			final double remaningTime = (System.currentTimeMillis() - startInferenceTime)
					/ epoch / 1000d * (maxEpoch - epoch);
			log("Step: %d di %d, tempo rimanente %s", epoch + 1, maxEpoch,
					timeToString(remaningTime));
		}
	}

	@Override
	public Model runInference() {
		nUsers = socialNetwork.rowsCount();
		nItems = preferenceMatrix.columnsCount();

		if (nUsers < 240)
			factory = BuildMatrixFactoryOKT
					.getInstance(BuildMatrixFactoryOKT.BLAS);
		else
			factory = BuildMatrixFactoryOKT
					.getInstance(BuildMatrixFactoryOKT.UJMP);

		log("Generating model: users(%d), items(%d), features(%d)", nUsers,
				nItems, nFactors);

		startInferenceTime = System.currentTimeMillis();

		initializeHyperParams();
		initializeParams();

		long time = System.currentTimeMillis();
		System.out.println("Building support parameters...");
		initializeSupportParams();
		System.out.println("... done. Elapsed time: "
				+ (System.currentTimeMillis() - time));

		double cumulativeLogLikelihood = 0;
		int k = 0;

		for (int epoch = 0; epoch < maxEpoch; epoch++) {
			logNextStep(epoch);

			// Sample hyperparams
			hyperSample();

			// Start doing Gibbs updates over user and movie feature vectors
			// given hyperparams. dimensione di R_train: n_users x n_items
			for (int gibbs = 0; gibbs < nTrials; gibbs++) {
				System.out.println("Gibbs' iteration: " + gibbs + " over "
						+ nTrials);

				System.out.println("Sampling Y");
				time = System.currentTimeMillis();
				samplingY();
				System.out.println("Elapsed time: "
						+ (System.currentTimeMillis() - time));

				System.out.println("Sampling Ze");
				time = System.currentTimeMillis();
				samplingZe();
				System.out.println("Elapsed time: "
						+ (System.currentTimeMillis() - time));

				System.out.println("Sampling Zr");
				time = System.currentTimeMillis();
				samplingZr();
				System.out.println("Elapsed time: "
						+ (System.currentTimeMillis() - time));

				System.out.println("Sampling Theta");
				time = System.currentTimeMillis();
				samplingTheta();
				System.out.println("Elapsed time: "
						+ (System.currentTimeMillis() - time));

				System.out.println("Sampling Omega");
				time = System.currentTimeMillis();
				samplingOmega();
				System.out.println("Elapsed time: "
						+ (System.currentTimeMillis() - time));
			}

			final int pos = epoch % epochHistorySize;

			ThetaAll[pos] = Theta.getCopy();
			OmegaAll[pos] = Omega.getCopy();

			if (epoch % 10 == 0) {
				cumulativeLogLikelihood += computeIncrementalLogLikelihood(pos);
				log("LogLike al passo %d:\t%lf", epoch, cumulativeLogLikelihood
						/ ++k);
			}
		}

		log("\nGenerazione modello completata: " + "tempo richiesto %s",
				timeToString((System.currentTimeMillis() - startInferenceTime) / 100.0));

		return new PairwiseRankingModel(ThetaAll, OmegaAll, index);
	}

	/**
	 * Sample from item hyperparams (see paper for details)
	 */
	protected MatrixOKT[] samplingHyperParams(final MatrixOKT x,
			final MatrixOKT W0, final MatrixOKT mu0, final int beta0,
			final int ni0) {

		final int n = x.rowsCount();

		final int ni0star = ni0 + n;
		final int beta0star = beta0 + n;

		final MatrixOKT xAvg = factory.mean(x).transpose(true);
		final MatrixOKT Sx = factory.covariance(x);
		final MatrixOKT mu0MinusxAvg = mu0.sub(xAvg);

		final Summation summation = factory.getSummation();
		summation.add(W0.inverse(), Sx.mulMe(n));

		summation.add(mu0MinusxAvg.mul(mu0MinusxAvg.transpose(true)).mulMe(
				n * beta0 / beta0star));

		final MatrixOKT W0star = summation.getResult().inverse();
		final MatrixOKT mu0star = mu0.mul(beta0).sumMe(xAvg.mul(n))
				.divMe(beta0star);

		final MatrixOKT sigma = factory.wishart(W0star, ni0star, null);
		final MatrixOKT lam = factory.cholesky(sigma.mul(beta0star).inverse())
				.transpose(true);

		final MatrixOKT mu = lam.mul(factory.randn(nFactors, 1)).sumMe(mu0star);

		return new MatrixOKT[] { mu, sigma };
	}

	protected void samplingOmega() {
		final Summation sigmaSummation = factory.getSummation();
		final Summation muSummation = factory.getSummation();

		final MatrixOKT sigmaOmegaInvMulMuOmega = sigmaOmegaInv.mul(muOmega);

		double oldPercent = 0;

		for (int i = 0; i < nItems; ++i) {
			final double percent = (double) i / nItems;

			if (percent >= oldPercent + 0.1) {
				System.out.println("\tSampling: " + padd(percent * 100, 1)
						+ "%");
				oldPercent = percent;
			}

			sigmaSummation.reset();
			sigmaSummation.add(sigmaOmegaInv);

			muSummation.reset();
			muSummation.add(sigmaOmegaInvMulMuOmega);

			for (final int u : usersPerItem.get(i)) {
				final MatrixOKT thetaUT = Theta.rows(u);
				final MatrixOKT thetaU = thetaUT.transpose(true);
				final MatrixOKT thetaUSq = thetaU.mul(thetaUT);

				final HashSet<Integer> comparisonsOfUI = pairwiseComparisons[u]
						.get(i);

				if (comparisonsOfUI != null) {
					for (final int j : comparisonsOfUI) {
						muSummation.add(thetaUSq.mul(Omega.rows(j).transpose(
								true)));
						muSummation.add(thetaU.mul(Zr[u].get(i, j)));
					}

					sigmaSummation.add(thetaUSq.mulMe(comparisonsOfUI.size()));
				}
			}

			final MatrixOKT sigmaStarOmegaInv = sigmaSummation.getResult();
			final MatrixOKT sigmaStarOmega = sigmaStarOmegaInv.inverse();
			final MatrixOKT muStarOmega = sigmaStarOmega.mul(muSummation
					.getResult());

			// campionamento con decomposizione di Cholesky
			final MatrixOKT lam = factory.cholesky(sigmaStarOmega).transpose(
					true);

			Omega.putRow(i,
					lam.mul(factory.randn(nFactors, 1)).sumMe(muStarOmega));
		}

		System.out.println("\tSampling: 100.0% complete");
	}

	protected void samplingTheta() {
		final Summation sigmaSummation = factory.getSummation();
		final Summation muSummation = factory.getSummation();
		final MatrixOKT sigmaThetaInvMulMuTheta = sigmaThetaInv.mul(muTheta);

		double oldPercent = 0;

		for (int u = 0; u < nUsers; ++u) {
			final double percent = (double) u / nUsers;

			if (percent >= oldPercent + 0.1) {
				System.out.println("\tSampling: " + padd(percent * 100, 1)
						+ "%");
				oldPercent = percent;
			}

			sigmaSummation.reset();
			sigmaSummation.add(sigmaThetaInv);

			muSummation.reset();
			muSummation.add(sigmaThetaInvMulMuTheta);

			for (final Map.Entry<Integer, HashSet<Integer>> e : pairwiseComparisons[u]
					.entrySet()) {

				final int itemI = e.getKey();
				final MatrixOKT omegaI = Omega.rows(itemI);

				for (final int itemJ : e.getValue()) {
					final MatrixOKT deltaOmegaT = omegaI.sub(Omega.rows(itemJ));

					sigmaSummation.add(deltaOmegaT.transpose(true).mul(
							deltaOmegaT));
					muSummation.add(deltaOmegaT.mulMe(Zr[u].get(itemI, itemJ)));
				}
			}

			for (final int v : neighborsPerUser.get(u)) {
				final MatrixOKT thetaVT = Theta.rows(v);

				sigmaSummation.add(thetaVT.transpose(true).mul(thetaVT));

				if (u < v)
					muSummation.add(thetaVT.mulMe(Ze.get(u, v)));
				else
					muSummation.add(thetaVT.mulMe(Ze.get(v, u)));
			}

			final LinkedList<Integer> neighbors = Y.get(u);

			if (neighbors != null)
				for (final int v : neighbors) {
					final MatrixOKT thetaVT = Theta.rows(v);

					sigmaSummation.add(thetaVT.transpose(true).mul(thetaVT));

					if (u < v)
						muSummation.add(thetaVT.mulMe(Ze.get(u, v)));
					else
						muSummation.add(thetaVT.mulMe(Ze.get(v, u)));

				}

			// found simgaStarTheta
			// getResult() transposes the matrix. Why? XXX
			final MatrixOKT sigmaStarThetaInv = sigmaSummation.getResult();
			final MatrixOKT sigmaStarTheta = sigmaStarThetaInv.inverse();

			// found muStarTheta
			final MatrixOKT muStarTheta = sigmaStarTheta.mul(muSummation
					.getResult());

			// campionamento con decomposizione di Cholesky
			final MatrixOKT lam = factory.cholesky(sigmaStarTheta).transpose(
					true);

			Theta.putRow(u,
					lam.mul(factory.randn(nFactors, 1)).sumMe(muStarTheta));
		}

		System.out.println("\tSampling: 100.0% complete");
	}

	private void samplingY() {

		// campioniamo solo quei valori di Y per i quali unknownLinks(u,v) > 0
		final Indices[] unknownIndices = unknownLinks
				.find(0, MatrixOKT.GREATER);

		final int n = unknownIndices[0].count();
		double oldPercent = 0;

		// Reset Y
		for (final LinkedList<Integer> list : Y.values())
			list.clear();

		Y.clear();

		for (int j = 0; j < n; j++) {
			final double percent = (double) j / n;

			if (percent >= oldPercent + 0.1) {
				System.out.println("\tSampling: " + padd(percent * 100, 1)
						+ "%");
				oldPercent = percent;
			}

			final int u = unknownIndices[0].get(j);
			final int v = unknownIndices[1].get(j);

			// Anche senza trasposta il dot() funziona lo stesso
			final double prod = Theta.rows(u).dot(Theta.rows(v));

			final double eta = Math.log(gamma1 / gamma2);

			// Sigmoid XXX non mi convice il segno, forse è -beta
			final double prY = 1d / (1d + Math.exp(beta * prod * prod / 2d
					+ eta));

			// scegliamo se assegnare l'arco o meno a Y_uv nel seguente modo:
			// generiamo un numero random x tra 0 e 1, se prY e' > x allora Y_uv
			// avrà un arco altrimenti no
			if (rg.nextDouble() < prY) {
				LinkedList<Integer> list = Y.get(u);

				if (list == null) {
					list = new LinkedList<Integer>();
					Y.put(u, list);
				}

				list.add(v);
			}
		}

		System.out.println("\tSampling: 100.0% complete");
	}

	private void samplingZe() {
		double oldPercent = 0;

		for (int u = 0; u < nUsers; ++u) {
			final double percent = (double) u / nUsers;

			if (percent >= oldPercent + 0.1) {
				System.out.println("\tSampling: " + padd(percent * 100, 1)
						+ "%");
				oldPercent = percent;
			}

			final MatrixOKT thetaUT = Theta.rows(u);

			for (final int v : neighborsPerUser.get(u)) {

				// Anche senza trasposta il dot() funziona lo stesso
				final double avg = Theta.rows(v).dot(thetaUT);

				if (u < v)
					Ze.set(u, v, tnrg.next(0, Double.POSITIVE_INFINITY, avg, 1));
				else
					Ze.set(v, u, tnrg.next(0, Double.POSITIVE_INFINITY, avg, 1));
			}

			for (final int v : Y.get(u)) {

				// Anche senza trasposta il dot() funziona lo stesso
				final double avg = Theta.rows(v).dot(thetaUT);

				if (u < v)
					Ze.set(u, v, tnrg.next(Double.NEGATIVE_INFINITY, 0, avg, 1));
				else
					Ze.set(v, u, tnrg.next(Double.NEGATIVE_INFINITY, 0, avg, 1));
			}
		}

		System.out.println("\tSampling: 100.0% complete");
	}

	protected void samplingZr() {
		double oldPercent = 0;

		for (int u = 0; u < nUsers; ++u) {
			final double percent = (double) u / nUsers;

			if (percent >= oldPercent + 0.1) {
				System.out.println("\tSampling: " + padd(percent * 100, 1)
						+ "%");
				oldPercent = percent;
				System.gc();
			}

			final MatrixOKT thetaU = Theta.rows(u).transpose(true);

			// Find u's items
			for (final Map.Entry<Integer, HashSet<Integer>> e : pairwiseComparisons[u]
					.entrySet()) {

				final int itemI = e.getKey();
				final MatrixOKT omegaIT = Omega.rows(itemI);

				for (final int itemJ : e.getValue()) {
					final double avg = omegaIT.sub(Omega.rows(itemJ)).dot(
							thetaU);

					Zr[u].set(
							itemI,
							itemJ,
							preferenceMatrix.get(u, itemI) > preferenceMatrix
							.get(u, itemJ) ? tnrg.next(0,
									Double.POSITIVE_INFINITY, avg, 1) : tnrg
									.next(Double.NEGATIVE_INFINITY, 0, avg, 1));
				}
			}
		}

		System.out.println("\tSampling: 100.0% complete");
	}

	public void setIndex(BidimensionalIndex index) {
		this.index = index;
	}

	@Override
	public void setSocialFoldMatrix(final MatrixOKT matrix) {
		unknownLinks = matrix;
	}

	@Override
	public void setSocialNetworkMatrix(final MatrixOKT matrix) {
		socialNetwork = matrix;
	}

	@Override
	public void setTrainingMatrix(final MatrixOKT matrix) {
		preferenceMatrix = matrix;
	}
}