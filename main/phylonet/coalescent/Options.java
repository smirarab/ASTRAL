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
	private String outputFile;
	private int samplingrounds;
	private int polylimit;
	
	//OLD parameters
	private double DLbdWeigth;
	private double CS;
	private double CD;
	private double trim;

	public Options(boolean rooted, boolean extrarooted, 
			boolean exactSolution, boolean duploss, int alg, int addExtra,
			boolean outputCompletedGenes, boolean outSearch, boolean run,
			int branchannotation, double lambda, String outputFile, int samplingrounds,int polylimit, double trim) {
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
		this.setOutputFile(outputFile);
		this.setSamplingrounds(samplingrounds);
		this.setPolylimit(polylimit);
		this.trim = trim;
	}

	public boolean isRooted() {
		return rooted;
	}


	public boolean isExtrarooted() {
		return extrarooted;
	}


	public boolean isExactSolution() {
		return exactSolution;
	}

	public boolean isDuploss() {
		return duploss;
	}


	public double getTrim() {
		return trim;
	}

	public void setTrim(double trim) {
		this.trim = trim;
	}

	public int getAlg() {
		return alg;
	}

	public int getAddExtra() {
		return addExtra;
	}

	/*public void setAddExtra(int addExtra) {
		this.addExtra = addExtra;
	}
	public void setAlg(int alg) {
		this.alg = alg;
	}
	public void setDuploss(boolean duploss) {
		this.duploss = duploss;
	}
	public void setExactSolution(boolean exactSolution) {
		this.exactSolution = exactSolution;
	}
	public void setExtrarooted(boolean extrarooted) {
		this.extrarooted = extrarooted;
	}
	
	public void setRooted(boolean rooted) {
		this.rooted = rooted;
	}
	public void setOutputCompletedGenes(boolean outputCompletedGenes) {
		this.outputCompletedGenes = outputCompletedGenes;
	}
	public void setOutputSearchSpace(boolean outSearch) {
		this.outputSearchSpace = outSearch;
	}
	public void setRunSearch(boolean run) {
		this.runSearch = run;
	}

	public void setBranchannotation(int branchannotation) {
		this.branchannotation = branchannotation;
	}

	*/

	public boolean isOutputCompletedGenes() {
		return outputCompletedGenes;
	}


	public boolean isOutputSearchSpace() {
		return outputSearchSpace;
	}

	
	public boolean isRunSearch() {
		return runSearch;
	}

	
	public int getBranchannotation() {
		return branchannotation;
	}


	public double getDLbdWeigth() {
		return DLbdWeigth;
	}


	public double getCS() {
		return CS;
	}

	public void setCS(double cS) {
		CS = cS;
	}
	public void setDLbdWeigth(double dLbdWeigth) {
		DLbdWeigth = dLbdWeigth;
	}

	public void setCD(double cD) {
		CD = cD;
	}
	
	public double getCD() {
		return CD;
	}

	public double getLambda() {
		return lambda;
	}

	public void setLambda(double lambda) {
		this.lambda = lambda;
	}

	public String getOutputFile() {
		return outputFile;
	}

	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}

	public Integer getSamplingrounds() {
		return samplingrounds;
	}

	public void setSamplingrounds(Integer samplingrounds) {
		this.samplingrounds = samplingrounds;
	}

	public int getPolylimit() {
		return polylimit;
	}

	public void setPolylimit(int polylimit) {
		this.polylimit = polylimit;
	}
}