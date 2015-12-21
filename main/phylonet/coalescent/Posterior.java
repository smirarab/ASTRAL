package phylonet.coalescent;

import cern.jet.stat.*;

public class Posterior extends cern.jet.math.Constants{
	private double m1;
	private double m2;
	private double m3;
	private double n;
	private double posterior;
	private static double LOG2 = Math.log(2.);
	public Posterior(double ft1, double ft2, double ft3, double nt){
	 	m1 = ft1*nt/(ft1+ft2+ft3);
		m2 = ft2*nt/(ft1+ft2+ft3);
		m3 = nt - m1 - m2;
		n  = nt;
		posterior=post();
	}
	public double getPost(){
		return posterior;
	}
	public  String toString(){
		StringBuffer out = new StringBuffer();
		out.append(posterior);	
		return  out.toString();
	}	
	public double betaRatio(double alpha1, double beta1, double alpha2, double beta2){
		double g = Gamma.logGamma(alpha1)+Gamma.logGamma(beta1)-Gamma.logGamma(alpha2)-Gamma.logGamma(beta2);
		return g;
	}
	public double G(double x, double nt){
		double g = 1- Gamma.incompleteBeta(x+1,nt-x+1,1./3.);
		if (g<=MACHEP){
				g = 0.;
		}
		return g;
	}
	public double post(){
		double g2 = LOG2*(m2-m1)+betaRatio(m2+1,n-m2+1,m1+1,n-m1+1);
		double g3 = LOG2*(m3-m1)+betaRatio(m3+1,n-m3+1,m1+1,n-m1+1);
		g2 = Math.exp(g2)*G(m2,n);
		g3 = Math.exp(g3)*G(m3,n);

		double g = G(m1,n)/(G(m1,n)+g2+g3);
		if (Double.isNaN(g) || Double.isInfinite(g)){
			if (3*m1>n) g=1.;
			else g=0;
		}
		/*if (g == 0){
			double v1 = G(m1,n)/(G(m1,n)+G(m2,n)*Functions.pow.apply(2.,m2-m1)*Gamma.beta(m2+1,n-m2+1)/Gamma.beta(m1+1,n-m1+1)+G(m3,n)*Functions.pow.apply(2.,m3-m1)*Gamma.beta(m3+1,n-m3+1)/Gamma.beta(m1+1,n-m1+1));
			if (Double.isNaN(v1)){
				if (m1/n>1./3.) g=1.;
				else g=0;
			}
		}*/
		posterior = g;
		return posterior;
	}    

	 public static void main(String[] args) {
		double m1 = Double.parseDouble(args[0]);
		double m2 = Double.parseDouble(args[1]);
		double m3 = Double.parseDouble(args[2]);
		double n  = Double.parseDouble(args[3]);
		Posterior a = new Posterior(m1,m2,m3,n);
		System.out.println(a.toString());
        }
}
