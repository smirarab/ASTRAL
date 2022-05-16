package phylonet.coalescent;

public class Options {
	private boolean rooted;
	private boolean extrarooted;
	private boolean exactSolution;
	private int alg;
	private int addExtra;
	private boolean outputCompletedGenes;
	private boolean outputCompatibledGenes;
	private boolean compatibleNorun;
	private boolean outputSearchSpace;
	private boolean runSearch;
	private int branchannotation; 
	private double lambda;
	private String outputFile;

	private int samplingrounds;
	private int polylimit;
	private String freqOutputPath;
	//OLD parameters
	private double DLbdWeigth;
	private double CS;
	private double CD;
	private Integer minLeaves;
	private Integer geneRepeat;
	private boolean removeExtraTree;
	private boolean subUnitBranchLength;

	private boolean ustarDist;

	public Options(boolean rooted, boolean extrarooted, 
			boolean exactSolution, int alg, int addExtra,
			boolean outputCompletedGenes, boolean outSearch, boolean run,
			int branchannotation, double lambda, String outputFile, int samplingrounds, int polylimit,
			double trim, String freqOutputPath, Integer minimumLeaves, Integer geneRepeat, boolean removeExtraTree, 
			boolean useInnerNodeDist, boolean subUnitBranchLength, 
			boolean outputCompatibledGenes, boolean compatibleNorun) {

		this.rooted = rooted;
		this.extrarooted = extrarooted;
		this.exactSolution = exactSolution;
		this.alg = alg;
		this.addExtra = addExtra;
		this.outputCompletedGenes = outputCompletedGenes;
		this.outputCompatibledGenes = outputCompatibledGenes;
		this.compatibleNorun = compatibleNorun;
		this.outputSearchSpace = outSearch;
		this.runSearch = run;
		this.branchannotation = branchannotation;
		this.setLambda(lambda);
		this.setOutputFile(outputFile);
		this.setSamplingrounds(samplingrounds);
		this.setPolylimit(polylimit);
		this.freqOutputPath = freqOutputPath;
		this.setMinLeaves(minimumLeaves);
		this.setGeneRepeat(geneRepeat);
		this.setUstarDist(useInnerNodeDist);
		this.setRemoveExtraTree(removeExtraTree);
		this.setSubUnitBranchLength(subUnitBranchLength);
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

	public boolean isOutputCompatibledGenes() {
		return outputCompatibledGenes;
	}
	
	public boolean isCompatibleNorun() {
		return compatibleNorun;
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

	public boolean isUstarDist() {
		return ustarDist;
	}

	public void setUstarDist(boolean ustarDist) {
		this.ustarDist = ustarDist;
	}

	public int getPolylimit() {
		return polylimit;
	}

	public void setPolylimit(int polylimit) {
		this.polylimit = polylimit;
	}

	public boolean isRemoveExtraTree() {
		return removeExtraTree;
	}

	public void setRemoveExtraTree(boolean removeExtraTree) {
		this.removeExtraTree = removeExtraTree;
	}

	public String getFreqOutputPath() {
		return freqOutputPath;
	}

	public void setFreqOutputPath(String freqOutputPath) {
		this.freqOutputPath = freqOutputPath;
	}

	public Integer getMinLeaves() {
		return minLeaves;
	}

	public void setMinLeaves(Integer minLeaves) {
		this.minLeaves = minLeaves;
	}

	public Integer getGeneRepeat() {
		return geneRepeat;
	}

	public void setGeneRepeat(Integer geneRepeat) {
		this.geneRepeat = geneRepeat;
	}

	public boolean isSubUnitBranchLength() {
		return subUnitBranchLength;
	}

	public void setSubUnitBranchLength(boolean subUnitBranchLength) {
		this.subUnitBranchLength = subUnitBranchLength;
	}
}
