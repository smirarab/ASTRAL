package phylonet.coalescent;
import cern.jet.stat.*;
import cern.jet.math.*;

public class Posterior {
	private double m1;
	private double m2;
	private double m3;
	private double n;
	private double posterior;
	public Posterior(){
	}
	public Posterior(double ft1, double ft2, double ft3, double nt){
	 	m1 = ft1*nt/(ft1+ft2+ft3);
		m2 = ft2*nt/(ft1+ft2+ft3);
		m3 = ft3*nt/(ft1+ft2+ft3);;
		n  = nt;
		posterior=post();
	}
	public double getPost(){
		return posterior;
	}
	public String toString(){
		return("Value of m1: "+ m1 +"\n Value of m2: "+ m2 +"\n Value of n: "+ n+ "\n Value of posterior: " + posterior+"\n" );
	}
	public double betaRatio(double alpha1, double beta1, double alpha2, double beta2){
		double g = Gamma.logGamma(alpha1)+Gamma.logGamma(beta1)+Gamma.logGamma(alpha2+beta2)-(Gamma.logGamma(alpha2)+Gamma.logGamma(beta2)+Gamma.logGamma(alpha1+beta1));
		return Functions.exp.apply(g);
	}
	public double G(double x, double nt){
		return (1-Gamma.incompleteBeta(x+1,n-x+1,1./3.));
	}
	public double post(){
		double g = G(m1,n)/(
				G(m1,n)+Functions.pow.apply(2.,m2-m1)*G(m2,n)*betaRatio(m2+1.,n-m2+1.,m1+1.,n-m1+1.)+Functions.pow.apply(2.,m3-m1)*G(m3,n)*betaRatio(m3+1,n-m3+1,m1+1,n-m1+1));
		if (Double.isNaN(g)){
			if (m1/n>1./3.) g=1.;
			else g=0;
		}
		if (g == 0){
			double v1 = G(m1,n)/(G(m1,n)+G(m2,n)*Functions.pow.apply(2.,m2-m1)*Gamma.beta(m2+1,n-m2+1)/Gamma.beta(m1+1,n-m1+1)+G(m3,n)*Functions.pow.apply(2.,m3-m1)*Gamma.beta(m3+1,n-m3+1)/Gamma.beta(m1+1,n-m1+1));
			if (Double.isNaN(v1)){
				if (m1/n>1./3.) g=1.;
				else g=0;
			}
		}
		posterior = g;
		return posterior;
	}
      	
}
