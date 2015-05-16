package com.oktsrl.utils;

import java.io.FileInputStream;
import java.util.Properties;

import com.oktsrl.MatrixOKT;
import com.oktsrl.Model;
import com.oktsrl.inferencers.impl.PairwiseRankingInferencer;

public class BayesianLauncher {

	public static void main(String[] args) throws Exception {
		System.out.println("Reading parameters...");

		final String conf = args[0];

		Properties p = new Properties();
		p.load(new FileInputStream(conf));

		final String net = p.getProperty("net");
		final String rat = p.getProperty("rat");
		final String unk = p.getProperty("unk");
		final String ind = p.getProperty("ind");
		final String infConfFile = p.getProperty("infConfFile");
		final String modelOutputFile = p.getProperty("modelOutputFile");

		final SerializatorIOManager io = new SerializatorIOManager();

		System.out.println("Loading data...");
		final MatrixOKT socialNetwork = (MatrixOKT) io.loadObject(net);
		final MatrixOKT preferenceMatrix = (MatrixOKT) io.loadObject(rat);
		final MatrixOKT unknowLinks = (MatrixOKT) io.loadObject(unk);
		final BidimensionalIndex index = (BidimensionalIndex) io
				.loadObject(ind);

		p = new Properties();
		p.load(new FileInputStream(infConfFile));
		final Settings s = new Settings(p);

		System.out.println("Running inference...");
		final PairwiseRankingInferencer inf = new PairwiseRankingInferencer(
				s);

		inf.setSocialFoldMatrix(unknowLinks);
		inf.setIndex(index);
		inf.setSocialNetworkMatrix(socialNetwork);
		inf.setTrainingMatrix(preferenceMatrix);

		final Model m = inf.runInference();

		System.out.println("Storing model...");
		io.storeObject(m, modelOutputFile);
	}
}