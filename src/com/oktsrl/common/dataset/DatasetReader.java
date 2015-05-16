package com.oktsrl.common.dataset;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Properties;
import java.util.Random;
import java.util.StringTokenizer;

import com.oktsrl.BuildMatrixFactoryOKT;
import com.oktsrl.MatrixFactoryOKT;
import com.oktsrl.MatrixOKT;
import com.oktsrl.utils.Action;
import com.oktsrl.utils.BidimensionalIndex;
import com.oktsrl.utils.EdgeType;
import com.oktsrl.utils.SerializatorIOManager;

public class DatasetReader {

	private static BidimensionalIndex buildFinalIdex(
			final BidimensionalIndex index,
			final HashMap<Integer, HashSet<Integer>> links,
			final HashMap<Integer, LinkedList<Action>> ratings) {

		final BidimensionalIndex newIndex = new BidimensionalIndex();

		for (final Map.Entry<Integer, HashSet<Integer>> e : links.entrySet()) {
			newIndex.addRow(index.rowId(e.getKey()));

			for (final Integer i : e.getValue())
				newIndex.addRow(index.rowId(i));
		}

		for (final Map.Entry<Integer, LinkedList<Action>> e : ratings
				.entrySet()) {
			newIndex.addRow(index.rowId(e.getKey()));

			for (final Action a : e.getValue())
				newIndex.addColumn(index.columnId(a.getItem()));
		}

		return newIndex;
	}

	private static MatrixOKT[] buildMatrixes(final BidimensionalIndex oldIndex,
			final BidimensionalIndex newIndex,
			final HashMap<Integer, HashSet<Integer>> links,
			final HashMap<Integer, LinkedList<Action>> ratings, final Random r) {

		final int nUsers = newIndex.nRows();
		final int nItems = newIndex.nColumns();

		MatrixFactoryOKT factory;

		if (nUsers < 240)
			factory = BuildMatrixFactoryOKT
					.getInstance(BuildMatrixFactoryOKT.BLAS);
		else
			factory = BuildMatrixFactoryOKT
					.getInstance(BuildMatrixFactoryOKT.UJMP);

		final MatrixOKT socialNetwork = factory.sparse(nUsers, nUsers);
		int nnz = 0;

		final int oldNUsers = oldIndex.nRows();

		for (int u = 0; u < oldNUsers; ++u) {
			final HashSet<Integer> followees = links.get(u);

			if (followees != null)
				for (final Integer v : followees) {
					final int newU = newIndex.rowIndex(oldIndex.rowId(u));
					final int newV = newIndex.rowIndex(oldIndex.rowId(v));
					final double oldEdge = socialNetwork.get(newU, newV);

					if (oldEdge == EdgeType.NO_EDGE) {
						socialNetwork.set(newU, newV, EdgeType.DIRECTED_EDGE);
						socialNetwork.set(newV, newU, EdgeType.INVERSE_EDGE);
					} else if (oldEdge == EdgeType.INVERSE_EDGE) {
						socialNetwork.set(newU, newV, EdgeType.DOUBLE_EDGE);
						socialNetwork.set(newV, newU, EdgeType.DOUBLE_EDGE);
					}

					++nnz;
				}
		}

		final MatrixOKT preferenceMatrix = factory.sparse(nUsers, nItems);

		for (final LinkedList<Action> actions : ratings.values())
			for (final Action action : actions)
				preferenceMatrix.set(newIndex.rowIndex(oldIndex.rowId(action
						.getUser())), newIndex.columnIndex(oldIndex
						.columnId(action.getItem())), action.getRating());

		int count = nnz << 1;
		final MatrixOKT unknownLinks = factory.sparse(nUsers, nUsers);

		while (count > 0) {
			final int row = r.nextInt(nUsers);
			final int column = r.nextInt(nUsers);

			if (unknownLinks.get(row, column) == EdgeType.NO_EDGE
					&& socialNetwork.get(row, column) == EdgeType.NO_EDGE) {
				unknownLinks.set(row, column, 1);
				count--;
			}
		}

		return new MatrixOKT[] { socialNetwork, preferenceMatrix, unknownLinks };
	}

	private static void filterMatrices(
			final HashMap<Integer, HashSet<Integer>> links,
			final HashMap<Integer, LinkedList<Action>> ratings,
			final int netMinEdges, final int ratMinRatings, final int nUsers) {

		int jj = 0;
		int removed = 0;
		final HashSet<Integer> removedList = new HashSet<Integer>();
		final HashSet<Integer> toRemoveList = new HashSet<Integer>();

		do {
			++jj;

			if (jj % 100 == 0)
				System.out.println("#clean: " + jj + ", total removed: "
						+ removed);

			toRemoveList.clear();

			for (int u = 0; u < nUsers; ++u) {
				if (removedList.contains(u))
					continue;

				final HashSet<Integer> followees = links.get(u);
				final LinkedList<Action> actions = ratings.get(u);

				if (followees == null || followees.size() < netMinEdges
						|| actions == null || actions.size() < ratMinRatings) {
					toRemoveList.add(u);
					removedList.add(u);
				}
			}

			removed += toRemoveList.size();
			removeUsers(toRemoveList, links, ratings, nUsers);
		} while (!toRemoveList.isEmpty());

		if (jj % 100 != 0)
			System.out.println("#clean: " + jj + ", total removed: " + removed);
	}

	public static void main(final String[] args) throws Exception {
		final String handlerFileName = args[0];

		final Properties properties = new Properties();
		properties.load(new FileInputStream(handlerFileName));

		final String netFile = properties.getProperty("network");
		final String ratingFile = properties.getProperty("ratings");
		final boolean netSkipFirstLine = properties.getProperty(
				"netSkipFirstLine").equalsIgnoreCase("true");
		final boolean ratSkipFirstLine = properties.getProperty(
				"ratSkipFirstLine").equalsIgnoreCase("true");
		final int netSource = Integer.parseInt(properties
				.getProperty("netSourceIndex"), 10);
		final int netDestination = Integer.parseInt(properties
				.getProperty("netDestinationIndex"), 10);
		final int netMinEdges = Integer.parseInt(properties
				.getProperty("netMinEdges"), 10);
		final int ratUser = Integer.parseInt(properties
				.getProperty("ratUserIndex"), 10);
		final int ratItem = Integer.parseInt(properties
				.getProperty("ratItemIndex"), 10);
		final int ratRating = Integer.parseInt(properties
				.getProperty("ratRatingIndex"), 10);
		final int ratMinRatings = Integer.parseInt(properties
				.getProperty("ratMinRatings"), 10);

		final String netOut = properties.getProperty("netOut");
		final String ratOut = properties.getProperty("ratOut");
		final String unkOut = properties.getProperty("unkOut");
		final String indOut = properties.getProperty("indOut");

		final long seed = Long.parseLong(properties.getProperty("seed"), 10);
		final Random r = new Random(seed);

		final BidimensionalIndex index = new BidimensionalIndex();

		System.out.println("reading network...");
		final HashMap<Integer, HashSet<Integer>> links = readNetwork(netFile,
				index, netSkipFirstLine, netSource, netDestination);

		System.out.println("reading ratings...");
		final HashMap<Integer, LinkedList<Action>> ratings = readRatings(
				ratingFile, index, ratSkipFirstLine, ratUser, ratItem,
				ratRating);

		System.out.println("initial");
		System.out.println("nUsers: " + index.nRows());
		System.out.println("nItems: " + index.nColumns());

		System.out.println("filtering...");
		filterMatrices(links, ratings, netMinEdges, ratMinRatings, index
				.nRows());

		System.out.println("building index...");
		final BidimensionalIndex newIndex = buildFinalIdex(index, links,
				ratings);

		System.out.println("building matrixes...");
		final MatrixOKT[] matrixes = buildMatrixes(index, newIndex, links,
				ratings, r);

		System.out.println("storing output...");
		storeInferenceInput(netOut, matrixes[0], ratOut, matrixes[1], unkOut,
				matrixes[2], indOut, newIndex);

		System.out.println("final");
		System.out.println("nUsers: " + newIndex.nRows());
		System.out.println("nItems: " + newIndex.nColumns());

		System.out.println("... end");
	}

	private static HashMap<Integer, HashSet<Integer>> readNetwork(
			final String netFile, final BidimensionalIndex index,
			final boolean netSkipFirstLine, final int netSource,
			final int netDestination) throws Exception {

		final BufferedReader br = new BufferedReader(new FileReader(netFile));
		String line;
		StringTokenizer st;

		if (netSkipFirstLine)
			br.readLine();

		final HashMap<Integer, HashSet<Integer>> links = new HashMap<Integer, HashSet<Integer>>(
				1024);

		while ((line = br.readLine()) != null) {
			st = new StringTokenizer(line, " \t\n\r\f,;");

			int jj = 0;
			int source = 0, destination = 0;

			while (st.hasMoreTokens()) {
				final int i = Integer.parseInt(st.nextToken(), 10);

				if (jj == netSource)
					source = i;

				if (jj == netDestination)
					destination = i;

				++jj;
			}

			final int sourceIndex = index.addRow(source);
			final int destinationIndex = index.addRow(destination);

			HashSet<Integer> list = links.get(sourceIndex);

			if (list == null) {
				list = new HashSet<Integer>();
				links.put(sourceIndex, list);
			}

			list.add(destinationIndex);
		}

		return links;
	}

	private static HashMap<Integer, LinkedList<Action>> readRatings(
			final String ratingFile, final BidimensionalIndex index,
			final boolean ratSkipFirstLine, final int ratUser,
			final int ratItem, final int ratRating) throws Exception {

		final BufferedReader br = new BufferedReader(new FileReader(ratingFile));
		String line;
		StringTokenizer st;

		if (ratSkipFirstLine)
			br.readLine();

		final HashMap<Integer, LinkedList<Action>> ratings = new HashMap<Integer, LinkedList<Action>>(
				1024);

		while ((line = br.readLine()) != null) {
			st = new StringTokenizer(line, " \t\n\r\f,;");

			int jj = 0;
			int user = 0;
			int item = 0;
			double rating = 0;

			while (st.hasMoreTokens()) {
				final String s = st.nextToken();

				if (jj == ratUser)
					user = Integer.parseInt(s, 10);

				if (jj == ratItem)
					item = Integer.parseInt(s, 10);

				if (jj == ratRating)
					rating = Double.parseDouble(s);

				++jj;
			}

			final int userIndex = index.addRow(user);
			final int itemIndex = index.addColumn(item);

			LinkedList<Action> list = ratings.get(userIndex);

			if (list == null) {
				list = new LinkedList<Action>();
				ratings.put(userIndex, list);
			}

			list.add(new Action(userIndex, itemIndex, rating));
		}

		return ratings;
	}

	private static void removeUsers(final HashSet<Integer> removedUsers,
			final HashMap<Integer, HashSet<Integer>> links,
			final HashMap<Integer, LinkedList<Action>> ratings, final int nUsers) {

		for (final Integer u : removedUsers) {
			ratings.remove(Integer.valueOf(u));
			links.remove(Integer.valueOf(u));
		}

		for (int v = 0; v < nUsers; ++v) {
			final HashSet<Integer> followees = links.get(v);

			if (followees == null)
				continue;

			followees.removeAll(removedUsers);
		}
	}

	private static void storeInferenceInput(final String netOut,
			final MatrixOKT net, final String ratOut, final MatrixOKT rat,
			final String unkOut, final MatrixOKT unk, final String indOut,
			final BidimensionalIndex index) throws Exception {

		final SerializatorIOManager m = new SerializatorIOManager();

		m.storeObject(net, netOut);
		m.storeObject(rat, ratOut);
		m.storeObject(unk, unkOut);
		m.storeObject(index, indOut);
	}
}