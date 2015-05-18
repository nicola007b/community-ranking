package com.oktsrl.inferencers;

import com.oktsrl.Inferencer;
import com.oktsrl.math.MatrixOKT;

public interface BayesianInferencer extends Inferencer {

	public void setSocialFoldMatrix(MatrixOKT matrix);

	public void setSocialNetworkMatrix(MatrixOKT matrix);

	public void setTrainingMatrix(MatrixOKT matrix);
}