package com.oktsrl.math.impls;

import com.oktsrl.math.MatrixFactoryOKT;
import com.oktsrl.math.MatrixOKT;

public abstract class AbstractMatrixOKT implements MatrixOKT {

	private static final long serialVersionUID = -3726717911004981351L;

	protected MatrixFactoryOKT factory;

	protected AbstractMatrixOKT(MatrixFactoryOKT factory) {
		this.factory = factory;
	}

	@Override
	public MatrixFactoryOKT getFactory() {
		return factory;
	}
}