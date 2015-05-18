package com.oktsrl.math;

public final class BuildMatrixFactoryOKT {
	public static final String BLAS = "BLAS";
	public static final String UJMP = "UJMP";

	public static MatrixFactoryOKT getInstance(String library) {
		try {
			final String className = String.format(
					"com.oktsrl.math.impls.%s.%sMatrixFactoryOKT",
					library.toLowerCase(), library);
			return (MatrixFactoryOKT) Class.forName(className).newInstance();
		} catch (final InstantiationException e) {
			e.printStackTrace();
		} catch (final IllegalAccessException e) {
			e.printStackTrace();
		} catch (final ClassNotFoundException e) {
			e.printStackTrace();
		}

		return null;
	}
}