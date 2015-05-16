package com.oktsrl.inferencers.impl;

import java.util.Random;

import com.oktsrl.Indices;
import com.oktsrl.MatrixFactoryOKT;
import com.oktsrl.MatrixOKT;
import com.oktsrl.Model;
import com.oktsrl.Summation;
import com.oktsrl.inferencers.BayesianInferencer;
import com.oktsrl.models.impl.GBModelImpl;
import com.oktsrl.utils.EdgeType;
import com.oktsrl.utils.Settings;
import com.oktsrl.utils.TruncatedNormalRandomGenerator;

/*
 * Permission is granted for anyone to copy, use, modify, or distribute this
 * program and accompanying programs and documents for any purpose, provided
 * this copyright notice is retained and prominently displayed, along with
 * a note saying that the original programs are available from our
 * web page.
 * The programs and documents are distributed without any warranty, express or
 * implied.  As the programs were written for research purposes only, they have
 * not been tested to the degree that would be advisable in any important
 * application.  All use of these programs is entirely at the user's own risk.
 */

/**
 * TODO che fa transpose(true)? Perch� il true? Qual � la differenza con
 * transopose() senza parametri?
 * 
 * TODO come funziona il factory.randn(num_feat, 1)?
 */
public class BayesianModelInferencerImpl implements BayesianInferencer {

	private static void log(final String message, final Object... params) {
		System.out.println(String.format(message, params));
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

	protected MatrixFactoryOKT factory;

	protected int nUsers;
	protected int nItems;

	protected double meanRating;
	protected MatrixOKT socialNetwork;
	protected MatrixOKT preferenceMatrix;
	protected MatrixOKT unknownLinks;

	protected MatrixOKT Theta, Omega, Y;
	protected MatrixOKT Ze;
	protected MatrixOKT[] Zr;

	protected MatrixOKT[] ThetaAll;
	protected MatrixOKT[] OmegaAll;

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
	protected int epochTopSize; // dimensionamente ultimi epoch validi

	protected long startInferenceTime;

	// Hyper parameters
	private MatrixOKT muTheta, muOmega;
	private MatrixOKT sigmaThetaInv, sigmaOmegaInv;
	private MatrixOKT W0_theta, W0_omega;
	private MatrixOKT mu0_theta, mu0_omega;
	private int beta0_theta, beta0_omega;
	private int ni0_theta, ni0_omega;

	public BayesianModelInferencerImpl(final Settings settings) {
		nFactors = settings.getInt("nTopics", 5);
		alpha = settings.getReal("alpha", 2d);
		beta = settings.getReal("beta", 2d);
		gamma1 = settings.getReal("gamma1", 2d);
		gamma2 = settings.getReal("gamma2", 2d);
		nTrials = settings.getInt("nTrials", 2);
		maxEpoch = settings.getInt("maxEpoch", 20);
		epochHistorySize = settings.getInt("epochHistorySize", 10);
		epochTopSize = settings.getInt("epochHistorySize", 10);

		final int seed = settings.getInt("seed", 101);

		tnrg = new TruncatedNormalRandomGenerator(seed);
		rg = new Random(seed);
	}

	public double getMeanRating() {
		return meanRating;
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

		// parameters of Inv-Whishart distribution for Theta (see paper for
		// details)
		W0_theta = factory.eye(nFactors);
		beta0_theta = 2;
		ni0_theta = nFactors;

		// XXX � riga o colonna? Deve essere colonna
		mu0_theta = factory.zeros(nFactors, 1);

		// parameters of Inv-Whishart distribution for Omega (see paper for
		// details)
		W0_omega = factory.eye(nFactors);
		beta0_omega = 2;
		ni0_omega = nFactors;

		// XXX � riga o colonna? Deve essere colonna
		mu0_omega = factory.zeros(nFactors, 1);
	}

	private void initializeParams() {
		Theta = factory.randn(nUsers, nFactors).mulMe(0.1);
		Omega = factory.randn(nItems, nFactors).mulMe(0.1);
		Y = factory.sparse(nUsers, nUsers);

		ThetaAll = new MatrixOKT[epochHistorySize];
		OmegaAll = new MatrixOKT[epochHistorySize];

		Ze = factory.sparse(nUsers, nUsers);
		Zr = new MatrixOKT[nUsers];

		for (int u = 0; u < nUsers; ++u)
			Zr[u] = factory.sparse(nItems, nItems);
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

	public Model runInference() {
		factory = socialNetwork.getFactory();
		nUsers = socialNetwork.rowsCount();
		nItems = preferenceMatrix.columnsCount();

		initializeHyperParams();
		initializeParams();

		// calcolo il rating medio XXX A che serve?
		meanRating = factory.mean(preferenceMatrix, 0, MatrixOKT.GREATER);

		log("Avvio generazione modello: users(%d), items(%d), feature(%d)",
				nUsers, nItems, nFactors);

		startInferenceTime = System.currentTimeMillis();

		double llk = 0;
		int k = 0;
		for (int epoch = 0; epoch < maxEpoch; epoch++) {
			logNextStep(epoch);

			// Sample hyperparams
			hyperSample();

			// Start doing Gibbs updates over user and movie feature vectors
			// given hyperparams. dimensione di R_train: n_users x n_items
			for (int gibbs = 0; gibbs < nTrials; gibbs++) {
				samplingY();
				samplingZe();
				samplingZr();
				samplingTheta();
				samplingOmega();
			}

			final int pos = epoch % epochTopSize;

			ThetaAll[pos] = Theta.getCopy();
			OmegaAll[pos] = Omega.getCopy();
			
			if (epoch % 10 == 0){
				llk += updateLLK(epoch);
				k = k+1;
				llk /= k;
				log("LogLike al passo %d:\t%lf",epoch,llk);

			}
			
		}

		final double time = (System.currentTimeMillis() - startInferenceTime) / 1000d;
		log("\nGenerazione modello completata: " + "tempo richiesto %s",
				timeToString(time));

		return new GBModelImpl(ThetaAll, OmegaAll, meanRating);
	}

	private double updateLLK(int pos) {
		// TODO Auto-generated method stub
		return 0;
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
		summation.add(mu0MinusxAvg.mulMe(mu0MinusxAvg.transpose(true)).mulMe(
				n * beta0 / beta0star));

		final MatrixOKT W0star = summation.getResult().inverse();
		final MatrixOKT mu0star = mu0.mul(beta0).sumMe(xAvg.mul(n)).divMe(
				beta0star);

		final MatrixOKT sigma = factory.wishart(W0star, ni0star, null);
		final MatrixOKT lam = factory.cholesky(sigma.mul(beta0star).inverse())
				.transpose(true);
		final MatrixOKT mu = lam.mulMe(factory.randn(nFactors, 1)).sumMe(
				mu0star);

		return new MatrixOKT[] { mu, sigma };
	}

	protected void samplingOmega() {
		final Summation sigmaSummation = factory.getSummation();
		final Summation muSummation = factory.getSummation();

		final MatrixOKT sigmaOmegaInvMulMuOmega = sigmaOmegaInv.mul(muOmega);

		for (int i = 0; i < nItems; ++i) {

			// restituisce in userIndex gli indici di riga per i quali gli
			// elementi della colonna i hanno un valore > 0:
			// tutti gli utenti che hanno espresso un rating > 0 per l'item i
			//
			final Indices userIndex = preferenceMatrix.findRowIndices(0, i,
					MatrixOKT.GREATER);

			sigmaSummation.reset();
			sigmaSummation.add(sigmaOmegaInv);

			muSummation.reset();
			muSummation.add(sigmaOmegaInvMulMuOmega);

			for (final int u : userIndex.toArray()) {
				final MatrixOKT thetaUT = Theta.rows(u);
				final MatrixOKT thetaU = thetaUT.transpose(true);
				final MatrixOKT thetaUSq = thetaU.mul(thetaUT);

				// Find u's items
				final Indices itemsOfU = preferenceMatrix.findColumnIndices(0,
						u, MatrixOKT.GREATER);

				// XXX ottimizzato se ordinato
				for (final int j : itemsOfU.toArray()) {
					if (j <= i)
						continue;

					muSummation
							.add(thetaUSq.mul(Omega.rows(j).transpose(true)));
					muSummation.add(thetaU.mul(Zr[u].get(i, j)));
				}

				sigmaSummation.add(thetaUSq.mulMe(itemsOfU.count()));
			}

			final MatrixOKT sigmaStarOmegaInv = sigmaSummation.getResult();

			final MatrixOKT sigmaStarOmega = sigmaStarOmegaInv.inverse();

			final MatrixOKT muStarOmega = sigmaStarOmega.mul(muSummation
					.getResult());

			// campionamento con decomposizione di Cholesky
			final MatrixOKT lam = factory.cholesky(sigmaStarOmega).transpose(
					true);

			Omega.putRow(i, lam.mul(factory.randn(nFactors, 1)).sumMe(
					muStarOmega));
		}
	}

	protected void samplingTheta() {
		final Summation sigmaSummation = factory.getSummation();
		final Summation muSummation = factory.getSummation();
		final MatrixOKT sigmaThetaInvMulMuTheta = sigmaThetaInv.mul(muTheta);

		for (int u = 0; u < nUsers; ++u) {
			sigmaSummation.reset();
			sigmaSummation.add(sigmaThetaInv);

			muSummation.reset();
			muSummation.add(sigmaThetaInvMulMuTheta);

			// restituisce in itemIndex gli indici di colonna per i quali gli
			// elementi della riga u hanno un valore > 0 item per i quali
			// l'utente u ha espresso un rating > 0
			final Indices itemIndex = preferenceMatrix.findColumnIndices(0, u,
					MatrixOKT.GREATER);

			final int[] itemIndexToArray = itemIndex.toArray();
			final int n = itemIndexToArray.length;

			for (int i = 0; i < n; ++i)
				for (int j = i + 1; j < n; ++j) {

					int itemI = itemIndexToArray[i];
					int itemJ = itemIndexToArray[j];

					if (itemI > itemJ) {
						final int tmp = itemI;
						itemI = itemJ;
						itemJ = tmp;
					}

					final MatrixOKT deltaOmegaT = Omega.rows(itemI).sub(
							Omega.rows(itemJ));

					muSummation.add(deltaOmegaT.mul(Zr[u].get(itemI, itemJ)));

					sigmaSummation.add(deltaOmegaT.transpose(true).mulMe(
							deltaOmegaT));
				}

			final Indices userIndex = socialNetwork.findColumnIndices(0, u,
					MatrixOKT.NOT);

			for (final int v : userIndex.toArray()) {
				if (v == u)
					continue;

				// XXX rows() resituisce una matrice nuova???
				final MatrixOKT thetaVT = Theta.rows(v);

				final double edge = socialNetwork.get(u, v);

				if (edge == EdgeType.DIRECTED_EDGE)
					muSummation.add(thetaVT.mul(Ze.get(u, v)));
				else if (edge == EdgeType.INVERSE_EDGE)
					muSummation.add(thetaVT.mul(Ze.get(v, u)));
				else if (edge == EdgeType.DOUBLE_EDGE) {
					muSummation.add(thetaVT.mul(Ze.get(u, v)));
					muSummation.add(thetaVT.mul(Ze.get(v, u)));
				}

				sigmaSummation.add(thetaVT.transpose(true).mul(thetaVT));
			}

			final Indices negativeIndex = Y.findColumnIndices(0, u,
					MatrixOKT.NOT);

			for (final int v : negativeIndex.toArray()) {
				if (v == u)
					continue;

				final MatrixOKT thetaVT = Theta.rows(v);

				final double edge = Y.get(u, v);

				if (edge == EdgeType.DIRECTED_EDGE)
					muSummation.add(thetaVT.mul(Ze.get(u, v)));
				else if (edge == EdgeType.INVERSE_EDGE)
					muSummation.add(thetaVT.mul(Ze.get(v, u)));
				else if (edge == EdgeType.DOUBLE_EDGE) {
					muSummation.add(thetaVT.mul(Ze.get(u, v)));
					muSummation.add(thetaVT.mul(Ze.get(v, u)));
				}

				sigmaSummation.add(thetaVT.transpose(true).mul(thetaVT));
			}

			// found simgaStarTheta
			final MatrixOKT simgaStarThetaInv = sigmaSummation.getResult();
			final MatrixOKT simgaStarTheta = simgaStarThetaInv.inverse();

			// found muStarTheta
			final MatrixOKT muStarTheta = simgaStarTheta.mul(muSummation
					.getResult().transpose(true));

			// campionamento con decomposizione di Cholesky
			final MatrixOKT lam = factory.cholesky(simgaStarTheta).transpose(
					true);

			Theta.putRow(u, lam.mul(factory.randn(nFactors, 1)).sumMe(
					muStarTheta));
		}
	}

	private void samplingY() {

		// campioniamo solo quei valori di Y per i quali unknownLinks(u,v) > 0
		final Indices[] unknownIndices = unknownLinks
				.find(0, MatrixOKT.GREATER);

		final int n = unknownIndices[0].count();

		for (int j = 0; j < n; j++) {
			final int u = unknownIndices[0].get(j);
			final int v = unknownIndices[1].get(j);

			final double prod = Theta.rows(u)
					.dot(Theta.rows(v).transpose(true));

			final double eta = Math.log(gamma1 / gamma2);
			// XXX Controlla che non sia -beta
			final double prY = factory.sigmoid(beta * prod * prod / 2d + eta);

			// scegliamo se assegnare l'arco o meno a Y_uv nel seguente modo:
			// generiamo un numero random x tra 0 e 1, se prY e' > x
			// allora Y_uv avr� un arco altrimenti no
			final boolean newEdge = rg.nextDouble() < prY;
			final double oldEdge = Y.get(u, v);

			if (newEdge) {
				if (oldEdge == EdgeType.NO_EDGE) {
					Y.set(u, v, EdgeType.DIRECTED_EDGE);
					Y.set(v, u, EdgeType.INVERSE_EDGE);
				} else if (oldEdge == EdgeType.INVERSE_EDGE) {
					Y.set(u, v, EdgeType.DOUBLE_EDGE);
					Y.set(v, u, EdgeType.DOUBLE_EDGE);
				}
			} else if (oldEdge == EdgeType.DIRECTED_EDGE) {
				Y.set(u, v, EdgeType.NO_EDGE);
				Y.set(v, u, EdgeType.NO_EDGE);
			} else if (oldEdge == EdgeType.DOUBLE_EDGE) {
				Y.set(u, v, EdgeType.INVERSE_EDGE);
				Y.set(v, u, EdgeType.DIRECTED_EDGE);
			}
		}
	}

	private void samplingZe() {
		for (int u = 0; u < nUsers; ++u) {

			// XXX cambiando la codifica EdgeType potrebbe non funzionare
			final Indices userIndex = socialNetwork.findColumnIndices(0, u,
					MatrixOKT.GREATER);
			final MatrixOKT thetaUT = Theta.rows(u);

			for (final int v : userIndex.toArray()) {
				if (v == u)
					continue;

				final double avg = thetaUT.dot(Theta.rows(v).transpose(true));
				Ze.set(u, v, tnrg.next(0, Double.POSITIVE_INFINITY, avg, 1));
			}

			// XXX cambiando la codifica EdgeType potrebbe non funzionare
			final Indices negativeIndex = Y.findColumnIndices(0, u,
					MatrixOKT.GREATER);

			for (final int v : negativeIndex.toArray()) {
				if (v == u)
					continue;

				final double avg = thetaUT.dot(Theta.rows(v).transpose(true));
				Ze.set(u, v, tnrg.next(Double.NEGATIVE_INFINITY, 0, avg, 1));
			}
		}
	}

	protected void samplingZr() {
		for (int u = 0; u < nUsers; ++u) {
			final MatrixOKT thetaUT = Theta.rows(u);

			// Find u's items
			final Indices itemsOfU = preferenceMatrix.findColumnIndices(0, u,
					MatrixOKT.GREATER);
			final int[] itemsOfUArray = itemsOfU.toArray();

			// XXX potresti ottimizzare se il vettore fosse ordinato
			for (final int i : itemsOfUArray) {
				final MatrixOKT omegaIT = Omega.rows(i);

				for (final int j : itemsOfUArray) {
					if (j <= i)
						continue;

					final double avg = thetaUT.dot(omegaIT.sub(Omega.rows(j))
							.transpose());
					Zr[u].set(i, j,
							preferenceMatrix.get(u, i) > preferenceMatrix.get(
									u, j) ? tnrg.next(0,
									Double.POSITIVE_INFINITY, avg, 1) : tnrg
									.next(Double.NEGATIVE_INFINITY, 0, avg, 1));
				}
			}
		}
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