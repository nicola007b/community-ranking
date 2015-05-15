package com.oktsrl.inferencers;

import com.oktsrl.Inferencer;
import com.oktsrl.MatrixOKT;

public interface BayesianInferencer extends Inferencer {
    public void setSocialNetworkMatrix(MatrixOKT matrix);
    public void setTrainingMatrix(MatrixOKT matrix);
    public void setSocialFoldMatrix(MatrixOKT matrix);

}
