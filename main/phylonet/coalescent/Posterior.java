package phylonet.coalescent;
import cern.jet.stat.*;

public class Posterior extends cern.jet.math.Constants{
	private double m1;
	private double m2;
	private double m3;
	private double n;
	private double posterior = -1;
	private static double LOG2 = Math.log(2.);
	final private static String MESSAGE = "This shouldn't haved happened. Please report error with the following numbers: ";
	private static final boolean DEBUG = true;
	
	public Posterior(double ft1, double ft2, double ft3, double nt){
	 	m1 = ft1*nt/(ft1+ft2+ft3);
		m2 = ft2*nt/(ft1+ft2+ft3);
		m3 = nt - m1 - m2;
		n  = nt;
		//posterior=post();
	}
	public double getPost(){
		if (this.posterior == -1) {
			this.posterior = post();
		}
		return posterior;
	}
	public  String toString(){
		StringBuffer out = new StringBuffer();
		out.append(this.getPost());	
		return  out.toString();
	}	
	public double betaRatio(double alpha1, double beta1, double alpha2, double beta2){
		return Gamma.logGamma(alpha1)+Gamma.logGamma(beta1)-Gamma.logGamma(alpha2)-Gamma.logGamma(beta2);
	}
	public double G(double x, double nt){
		double g = 1- Gamma.incompleteBeta(x+1,nt-x+1,1./3.);
		if (g<=MACHEP){
			g = 0.;
		}
		return g;
	}
	public double r(double mi){
		return Math.exp(LOG2*(mi-m1)+betaRatio(mi+1,n-mi+1,m1+1,n-m1+1));
	}
	
	public double rG(double mi) {
		double x = G(mi,n) * r (mi);
		if (Double.isNaN(x) || Double.isInfinite(x)) {
			if (mi*3<n) {
				return 0;
			} else if (mi*3==n) {
				if (m1 > mi) {
					return 0;
				} else {
					throw new RuntimeException(MESSAGE + "\n" + m1 +" "+ m2 +" "+ m3 +" ");
				}
			} else {
				return Double.POSITIVE_INFINITY;
			}
		}
		return x;
	}
	private double post(){
		
		if (this.DEBUG) {
		 System.out.println(m1 +" "+ m2 +" "+ m3);
		 System.out.println(G(m1,n));
		 System.out.println(G(m2,n));
		 System.out.println(G(m3,n));
		 System.out.println("r2: " + r(m2));
		 System.out.println("r3: " + r(m3));
		}
		
		double g2 = rG(m2);
		double g3 = rG(m3);

		double g = G(m1,n)/(G(m1,n)+g2+g3);
		
		if (Double.isInfinite(g)) {
			throw new RuntimeException(MESSAGE + "\n" + m1 +" "+ m2 +" "+ m3 +" ");
		}
		if (Double.isNaN(g)) {
			if (m1*3 < n) return g = 0;
			else throw new RuntimeException(MESSAGE + "\n" + m1 +" "+ m2 +" "+ m3 +" ");
		}
		/*	if (Double.isNaN(g) || Double.isInfinite(g)){
			if (m1>m2 & m1>m3) g=1.;
			else g=0;
		}*/
		//if (g == 0){
		//	double v1 = G(m1,n)/(G(m1,n)+G(m2,n)*Functions.pow.apply(2.,m2-m1)*Gamma.beta(m2+1,n-m2+1)/Gamma.beta(m1+1,n-m1+1)+G(m3,n)*Functions.pow.apply(2.,m3-m1)*Gamma.beta(m3+1,n-m3+1)/Gamma.beta(m1+1,n-m1+1));
		//	if (Double.isNaN(v1)){
		//		if (m1/n>1./3.) g=1.;
		//		else g=0;
		//	}
		//}
		posterior = g;
		return posterior;
	}    

	 public static void main(String[] args) {
	/*	double m1 = Double.parseDouble(args[0]);
		double m2 = Double.parseDouble(args[1]);
		double m3 = Double.parseDouble(args[2]);
		double n  = Double.parseDouble(args[3]);
		Posterior a = new Posterior(m1,m2,m3,n);
		System.out.println(a.toString());*/
		 double m1 = 100000-30000;
		 double m2 = 100000+5000;
		 double n =  300000;
		 double m3 = n-m1-m2;
		 Posterior p = new Posterior(10, 10, 10, n);
		 p.m1 = m1;p.m2 = m2;p.m3 = m3;
		 System.out.println(p.getPost());

	 }
}
