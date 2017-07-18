package phylonet.coalescent;

import cern.jet.stat.*;

/***
 * Computes branch length and posterior support
 * @author smirarab
 *
 */
public class Posterior extends cern.jet.math.Constants{
	private double f1;
	private double f2;
	private double f3;
	private double n;
	private double posterior = -1;
	private static double LOG2 = Math.log(2.);
	private double pValue = -1;
	final private static String MESSAGE = "This shouldn't have happened."
			+ " Maybe you set lambda too high or too low? "
			+ "Please report the error with the following numbers: ";
	private static final boolean DEBUG = false;
	private double lambda;
	
	public Posterior(double ft1, double ft2, double ft3, double nt, double lambda){
		this.f1 = ft1;
	 	this.f2 = ft2;
		this.f3 = ft3;
		this.n  = nt;
		this.lambda = lambda;
		//posterior=post();
	}
	public double getPvalue(){
		if (pValue == -1) {
			this.pValue = pvalue();
		}
		return pValue;
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
	public String toString2(){
		StringBuffer out = new StringBuffer();
		out.append(this.getPvalue());
		return out.toString();
	}
	public double betaRatio(double alpha1, double beta1, double alpha2, double beta2){
		double a = Gamma.logGamma(alpha1)+Gamma.logGamma(beta1)-Gamma.logGamma(alpha2)-Gamma.logGamma(beta2);
		return a;
	}
	public double G(double x, double nt){
		double g = 1- Gamma.incompleteBeta(x+1,nt-x+2*lambda,1./3.);
		if (g<=MACHEP){
			g = 0.;
		}
		return g;
	}
	public double r(double mi){
		double b = Math.exp(LOG2*(mi-f1)+betaRatio(mi+1,n-mi+2*lambda,f1+1,n-f1+2*lambda));
		return b;
	}
	
	public double rG(double mi) {
		double x = G(mi,n) * r (mi);
		if (Double.isNaN(x) || Double.isInfinite(x)) {
			if (mi*3<n) {
				return 0;
			} else if (mi*3>=n) {
				if (f1 > mi) {
					return 0;
				} else if (f1<mi){
					return Double.POSITIVE_INFINITY;
				} else {
					throw new RuntimeException(MESSAGE + "\n" + f1 +" "+ f2 +" "+ f3 +" "+lambda + " " + n);
				}
			} else {
				return Double.POSITIVE_INFINITY;
			}
		}
		return x;
	}
	
	public double branchLength(){
		Double bl = 0.;
		if ((f1/(n+2*lambda)) >1./3) {
			bl = -Math.log(1.5*(1.0-(f1/(n+2*lambda))));
		} else {
			bl = 0.0;
		}
		return bl;
	}
	private double pvalue(){
		double fThird = n/3.;
		double p;
		double x; 
		x=Math.pow((f1-fThird),2)/fThird+Math.pow((f2-fThird),2)/fThird+Math.pow((f3-fThird),2)/fThird;
		p = Probability.chiSquareComplemented(2,x);
		
		System.err.println( f1 + " " + f2 + " " + f3 + " " + n + " " + p );
		if (n<=15) {				
				System.err.println("Not enough genes");
				return -2;
		} else {
		}
		
		return p;
	}
	private double post(){
		
		if (this.DEBUG) {
		 System.out.println(f1 +" "+ f2 +" "+ f3);
		 System.out.println("G1: " + G(f1,n));
		 System.out.println("G2: " + G(f2,n));
		 System.out.println("G3: " + G(f3,n));
		 System.out.println("r2: " + r(f2));
		 System.out.println("r3: " + r(f3));
		}
		
		double g2 = rG(f2);
		double g3 = rG(f3);

		double g = G(f1,n)/(G(f1,n)+g2+g3);
		
		if (Double.isInfinite(g)) {
			throw new RuntimeException(MESSAGE + "\n" + f1 +" "+ f2 +" "+ f3 +" "+lambda + " " + n);
		}
		if (Double.isNaN(g)) {
			if (f1*3 < n) return g = 0;
			else throw new RuntimeException(MESSAGE + "\n" + f1 +" "+ f2 +" "+ f3 +" "+lambda + " " + n);
		}
		/*	if (Double.isNaN(g) || Double.isInfinite(g)){
			if (m1>m2 & m1>m3) g=1.;
			else g=0;
		}*/
		
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
		 double m1 = 5;
		 double m2 = 0;
		 double m3 = 0;
		 double n =  5;
		 Posterior p = new Posterior(m1, m2, m3, n, 0.5);
		 p.f1 = m1;p.f2 = m2;p.f3 = m3;
		 System.out.println(p.getPost());
		 System.out.println(p.getPvalue()); 

	 }
}
