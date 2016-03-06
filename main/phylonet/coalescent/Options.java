package phylonet.coalescent;

public class Options {
	private boolean rooted;
	private boolean extrarooted;
	private boolean exactSolution;
	private boolean duploss;
	private int alg;
	private int addExtra;
	private boolean outputCompletedGenes;
	private boolean outputSearchSpace;
	private boolean runSearch;
	private int branchannotation; 
	private double lambda;
	
	//OLD parameters
	private double DLbdWeigth;
	private double CS;
	private double CD;

	public Options(boolean rooted, boolean extrarooted, 
			boolean exactSolution, boolean duploss, int alg, int addExtra,
			boolean outputCompletedGenes, boolean outSearch, boolean run,
			int branchannotation, double lambda) {
		this.rooted = rooted;
		this.extrarooted = extrarooted;
		this.exactSolution = exactSolution;
		this.duploss = duploss;
		this.alg = alg;
		this.addExtra = addExtra;
		this.outputCompletedGenes = outputCompletedGenes;
		this.outputSearchSpace = outSearch;
		this.runSearch = run;
		this.branchannotation = branchannotation;
		this.setLambda(lambda);
	}

	public boolean isRooted() {
		return rooted;
	}

	public void setRooted(boolean rooted) {
		this.rooted = rooted;
	}

	public boolean isExtrarooted() {
		return extrarooted;
	}

	public void setExtrarooted(boolean extrarooted) {
		this.extrarooted = extrarooted;
	}

	public boolean isExactSolution() {
		return exactSolution;
	}

	public void setExactSolution(boolean exactSolution) {
		this.exactSolution = exactSolution;
	}

	public boolean isDuploss() {
		return duploss;
	}

	public void setDuploss(boolean duploss) {
		this.duploss = duploss;
	}

	public int getAlg() {
		return alg;
	}

	public void setAlg(int alg) {
		this.alg = alg;
	}

	public int getAddExtra() {
		return addExtra;
	}

	public void setAddExtra(int addExtra) {
		this.addExtra = addExtra;
	}

	public boolean isOutputCompletedGenes() {
		return outputCompletedGenes;
	}

	public void setOutputCompletedGenes(boolean outputCompletedGenes) {
		this.outputCompletedGenes = outputCompletedGenes;
	}

	public boolean isOutputSearchSpace() {
		return outputSearchSpace;
	}

	public void setOutputSearchSpace(boolean outSearch) {
		this.outputSearchSpace = outSearch;
	}

	public boolean isRunSearch() {
		return runSearch;
	}

	public void setRunSearch(boolean run) {
		this.runSearch = run;
	}

	public int getBranchannotation() {
		return branchannotation;
	}

	public void setBranchannotation(int branchannotation) {
		this.branchannotation = branchannotation;
	}

	public double getDLbdWeigth() {
		return DLbdWeigth;
	}

	public void setDLbdWeigth(double dLbdWeigth) {
		DLbdWeigth = dLbdWeigth;
	}

	public double getCS() {
		return CS;
	}

	public void setCS(double cS) {
		CS = cS;
	}

	public double getCD() {
		return CD;
	}

	public void setCD(double cD) {
		CD = cD;
	}

	public double getLambda() {
		return lambda;
	}

	public void setLambda(double lambda) {
		this.lambda = lambda;
	}
}