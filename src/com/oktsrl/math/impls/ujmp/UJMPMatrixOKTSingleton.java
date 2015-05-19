package com.oktsrl.math.impls.ujmp;

import org.ujmp.core.Matrix;

import com.oktsrl.math.MatrixFactoryOKT;

public class UJMPMatrixOKTSingleton extends UJMPMatrixOKT {

	private static final long serialVersionUID = 4751589345525613574L;

	public UJMPMatrixOKTSingleton(Matrix matrix, MatrixFactoryOKT factory) {
		super(null, factory);
	}

	public void setMatrix(Matrix matrix) {
		this.matrix = matrix;
	}
}