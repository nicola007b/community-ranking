package com.oktsrl;

public class PredictorEvaluator {
	private int numberOfOccurrences;
	private double mSE;
	private double mAE;
	private double mAPE;

	public PredictorEvaluator() {
		numberOfOccurrences = 0;
		mSE = 0.0d;
		mAE = 0.0d;
		mAPE = 0.0d;
	}

	public void addContribution(double realValue, double predictedValue) {
		numberOfOccurrences++;
		double delta = realValue - predictedValue;
		mSE += Math.pow(delta, 2);
		mAE += Math.abs(delta);
		if (realValue != 0)
			mAPE += (Math.abs(delta) / realValue);
	}

	public double getMAE() {
		return mAE / numberOfOccurrences;
	}

	public double getMAPE() {
		return mAPE / numberOfOccurrences;
	}

	public double getMSE() {
		return mSE / numberOfOccurrences;
	}

	public double getRMSE(){
		return Math.sqrt(mSE
				/ numberOfOccurrences);
	}
	
	public int getNumberOfOccurrences() {
		return numberOfOccurrences;
	}

	public void print() {

		String prints = "mse:"
				+ (mSE / numberOfOccurrences)
				+ "\nrmse:"
				+ (Math.sqrt(mSE
						/ numberOfOccurrences)) + "\nmae:"
				+ (mAE / numberOfOccurrences)
				+ "\nmape:" + (mAPE / numberOfOccurrences);

		System.out.println(prints);

	}

}

