package com.oktsrl.test;

import org.apache.commons.math3.distribution.NormalDistribution;

public class CumulativeNormalTest {

	public static void main(String[] args) {
		final NormalDistribution nd = new NormalDistribution();
		
		double i = -10;
		double delta = 0.1;
		do {
			System.out.println(i+":\t"+nd.cumulativeProbability(i)+"\t"+cdf(i)+"\t"+phi(i));
		i = i + delta;
		}
		while (i <= 10);
	}
	public static double cdf(double x){
		double XAbs = Math.abs(x);
		double Cumnorm;
		if (XAbs > 37)
				Cumnorm = 0;
		else {
			double Exponential = Math.exp(-(XAbs*XAbs)/ 2.0);
			double build;
			if (XAbs < 7.07106781186547) {
				build = 3.52624965998911E-02 * XAbs + 0.700383064443688;
				build = build * XAbs + 6.37396220353165;
				build = build * XAbs + 33.912866078383;
				build = build * XAbs + 112.079291497871;
				build = build * XAbs + 221.213596169931;
				build = build * XAbs + 220.206867912376;
				Cumnorm = Exponential * build;
				build = 8.83883476483184E-02 * XAbs + 1.75566716318264;
				build = build * XAbs + 16.064177579207;
				build = build * XAbs + 86.7807322029461;
				build = build * XAbs + 296.564248779674;
				build = build * XAbs + 637.333633378831;
				build = build * XAbs + 793.826512519948;
				build = build * XAbs + 440.413735824752;
				Cumnorm = Cumnorm / build;
			}
			else {
				build = XAbs + 0.65;
				build = XAbs + 4 / build;
				build = XAbs + 3 / build;
				build = XAbs + 2 / build;
				build = XAbs + 1 / build;
				Cumnorm = Exponential / build / 2.506628274631;
			}
			
		}
		if (x > 0) 
			Cumnorm = 1 - Cumnorm;
		return Cumnorm;
	}
	
	public static double phi(double x)
	{
	    // constants
	    double a1 =  0.254829592;
	    double a2 = -0.284496736;
	    double a3 =  1.421413741;
	    double a4 = -1.453152027;
	    double a5 =  1.061405429;
	    double p  =  0.3275911;
	 
	    // Save the sign of x
	    int sign = 1;
	    if (x < 0)
	        sign = -1;
	    x = Math.abs(x)/Math.sqrt(2.0);
	 
	    // A&S formula 7.1.26
	    double t = 1.0/(1.0 + p*x);
	    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*Math.exp(-x*x);
	 
	    return 0.5*(1.0 + sign*y);
	}
	 /*
	void testPhi()
	{
	    // Select a few input values
	    double x[] = 
	    {
	        -3, 
	        -1, 
	        0.0, 
	        0.5, 
	        2.1 
	    };
	 
	    // Output computed by Mathematica
	    // y = Phi[x]
	    double y[] = 
	    { 
	        0.00134989803163, 
	        0.158655253931, 
	        0.5, 
	        0.691462461274, 
	        0.982135579437 
	    };
	 
	        int numTests = sizeof(x)/sizeof(double);
	 
	    double maxError = 0.0;
	    for (int i = 0; i < numTests; ++i)
	    {
	        double error = fabs(y[i] - phi(x[i]));
	        if (error > maxError)
	            maxError = error;
	    }
	 
	        std::cout << "Maximum error: " << maxError << "\n";
	} 
	*/
}