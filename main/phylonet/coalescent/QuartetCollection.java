package phylonet.coalescent;

public class QuartetCollection {
	
//	public class Quartet {
//		int a,b,c,d;
//	}

	float [][] quartetWeights;
	float [][] quartSumWeights;
	
	int twoNmin3;
	
	public QuartetCollection(int n) {
		twoNmin3 = 2 * n - 3;
		quartetWeights = new float[n*(n-1)/2][n*(n-1)/2];
		quartSumWeights = new float[n*(n-1)/2][n];
	}
	
	private int index(int a, int b) {
		return(a*(twoNmin3 - a)>>2+b-1);
	}
	public float getQuartetWeight(int a, int b, int c, int d) {
		return quartetWeights[index(a,b)][index(c,d)];
	}
	
	public float getquartetSumWeigth(int a, int b, int c) {
		return quartSumWeights[index(a,b)][c];
	}
	
	public void setQuartetWeight (int a, int b, int c, int d, float w ) {
		quartetWeights[index(a,b)][index(c,d)] = w;
	}
	
	public void setSumQuartetWeight (int a, int b, int c, float w ) {
		quartSumWeights[index(a,b)][c] = w;
	}
}
