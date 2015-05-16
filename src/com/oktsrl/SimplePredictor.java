package com.oktsrl;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Properties;
import java.util.StringTokenizer;

import com.oktsrl.models.impl.PairwiseRankingModel;
import com.oktsrl.utils.SerializatorIOManager;

public class SimplePredictor {

	public static void main(String[] args) throws Exception {
		final String handlerFileName = args[0];
		final String modelFile = args[1];
		final String outputFile = args[2];

		final SerializatorIOManager io = new SerializatorIOManager();

		System.out.println("Loading model...");
		final PairwiseRankingModel rm = (PairwiseRankingModel) io.loadObject(modelFile);
		System.out.println("... done");

		final Properties properties = new Properties();
		properties.load(new FileInputStream(handlerFileName));

		final String ratingFile = properties.getProperty("ratings");
		final boolean ratSkipFirstLine = properties.getProperty(
				"ratSkipFirstLine").equalsIgnoreCase("true");
		final int ratUser = Integer.parseInt(
				properties.getProperty("ratUserIndex"), 10);
		final int ratItem = Integer.parseInt(
				properties.getProperty("ratItemIndex"), 10);

		System.out.println("reading entries and computing scores...");

		final BufferedReader br = new BufferedReader(new FileReader(ratingFile));
		final PrintWriter pw = new PrintWriter(new FileWriter(outputFile));

		pw.println("user\titem\tscore");

		String line;
		StringTokenizer st;

		if (ratSkipFirstLine)
			br.readLine();

		int kkk = 0;

		while ((line = br.readLine()) != null) {
			st = new StringTokenizer(line, " \t\n\r\f,;");

			int jj = 0;
			int userId = 0;
			int itemId = 0;

			while (st.hasMoreTokens()) {
				final String s = st.nextToken();

				if (jj == ratUser)
					userId = Integer.parseInt(s, 10);

				if (jj == ratItem)
					itemId = Integer.parseInt(s, 10);

				++jj;
			}

			final double score = rm.getScore(userId, itemId);

			pw.println(userId + "\t" + itemId + "\t" + score);

			++kkk;

			if (kkk % 1000 == 0)
				System.out.println("# " + kkk + " scores computed");
		}

		br.close();
		pw.close();
	}
}