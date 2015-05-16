package com.oktsrl.models.impl;

import com.oktsrl.MatrixOKT;
import com.oktsrl.Model;
import com.oktsrl.utils.BidimensionalIndex;

public final class PairwiseRankingModel implements Model {

	private static final long serialVersionUID = -4432232006908441346l;

	private final MatrixOKT[] ThetaAll;
	private final MatrixOKT[] OmegaAll;
	private final BidimensionalIndex index;

	public PairwiseRankingModel(MatrixOKT[] Omega_All, MatrixOKT[] Theta_All,
			BidimensionalIndex index) {
		OmegaAll = Omega_All;
		ThetaAll = Theta_All;
		this.index = index;
	}

	public double getScore(int userId, int itemId) {
		final int userIndex = index.rowIndex(userId);
		final int itemIndex = index.columnIndex(itemId);

		double sum = 0;
		final int nEpoch = OmegaAll.length;

		for (int epoch = 0; epoch < nEpoch; epoch++)
			sum += ThetaAll[epoch].rows(userIndex).dot(
					OmegaAll[epoch].rows(itemIndex));

		return sum / nEpoch;
	}
}