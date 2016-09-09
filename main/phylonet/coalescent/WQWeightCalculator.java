package phylonet.coalescent;


import java.util.Iterator;
import java.util.List;

import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;

class WQWeightCalculator extends AbstractWeightCalculator<Tripartition> {

	WQInference inference;
	private WQDataCollection dataCollection;
	private int n;

	public WQWeightCalculator(AbstractInference<Tripartition> inference) {
		super(false);
		this.dataCollection = (WQDataCollection) inference.dataCollection;
		this.inference = (WQInference) inference;
		this.n = GlobalMaps.taxonIdentifier.taxonCount();
	}

	class QuartetWeightTask implements ICalculateWeightTask<Tripartition>{

		private Tripartition trip;

		public QuartetWeightTask(Tripartition trip) {
			this.trip = trip;
		}

		long calculateMissingWeight() {
			// System.err.print("Calculating weight for: " + biggerSTB);
			long weight = 0l;
			for (int i=0; i < dataCollection.finalCounts.length; i++) {
				weight += sharedQuartetCount(
						trip,dataCollection.finalTripartitions[i]) * 
						dataCollection.finalCounts[i];
			}
			return weight;
		}

		private long F(int a,int b,int c) {
			if (a<0 || b<0 || c<0) {
				throw new RuntimeException("negative side not expected: "+a+" "+b+" "+c);
			}
			long ret = (a+b+c-3);
			ret *= a*b*c;
			return ret;
		}	

		long sharedQuartetCount(Tripartition that, Tripartition other) {

			int 
			I0 = that.cluster1.getBitSet().intersectionSize(other.cluster1.getBitSet()),
			I1 = that.cluster1.getBitSet().intersectionSize(other.cluster2.getBitSet()),
			I2 = that.cluster1.getBitSet().intersectionSize(other.cluster3.getBitSet()),
			I3 = that.cluster2.getBitSet().intersectionSize(other.cluster1.getBitSet()),
			I4 = that.cluster2.getBitSet().intersectionSize(other.cluster2.getBitSet()),
			I5 = that.cluster2.getBitSet().intersectionSize(other.cluster3.getBitSet()),
			I6 = that.cluster3.getBitSet().intersectionSize(other.cluster1.getBitSet()),
			I7 = that.cluster3.getBitSet().intersectionSize(other.cluster2.getBitSet()),
			I8 = that.cluster3.getBitSet().intersectionSize(other.cluster3.getBitSet());

			return  F(I0,I4,I8)+F(I0,I5,I7)+
					F(I1,I3,I8)+F(I1,I5,I6)+
					F(I2,I3,I7)+F(I2,I4,I6); 
		}

		/*
		 * The main function used for scoring a tripartition
		 *
		 */
		 long calculateWeightByTraversal() { 
			 long weight = 0;
			 int[]  allsides = null;
			 Iterator<STITreeCluster> tit = dataCollection.treeAllClusters.iterator();
			 boolean newTree = true;
			 int[][] children = new int[n+2][3];
			 int top = 0;
			 for (Integer gtb: dataCollection.geneTreesAsInts){
				 if (newTree) {
					 STITreeCluster all = tit.next();
					 allsides = new int [] {
							 trip.cluster1.getBitSet().intersectionSize(all.getBitSet()),
							 trip.cluster2.getBitSet().intersectionSize(all.getBitSet()),
							 trip.cluster3.getBitSet().intersectionSize(all.getBitSet())};
					 newTree = false;
				 }
				 if (gtb >= 0){
					 if (trip.cluster1.getBitSet().get(gtb)) {
						    children[top][0] = 1;
							children[top][1] = 0;
							children[top][2] = 0;
						} else if (trip.cluster2.getBitSet().get(gtb)) {
						    children[top][0] = 0;
							children[top][1] = 1;
							children[top][2] = 0;
						} else if (trip.cluster3.getBitSet().get(gtb)) {
						    children[top][0] = 0;
							children[top][1] = 0;
							children[top][2] = 1;
						} 
					 top++;
				 } else if (gtb == Integer.MIN_VALUE) {
					 top = 0;
					 newTree = true;
				 } else if (gtb == -2) {

					 top--;
					 
					 int newSides0 = children[top][0] + children[top-1][0];
					 int newSides1 = children[top][1] + children[top-1][1];
					 int newSides2 = children[top][2] + children[top-1][2];

					 int side3s0 = allsides[0] - newSides0;
					 int side3s1 = allsides[1] - newSides1;
					 int side3s2 = allsides[2] - newSides2;

					 weight += 
							 F(children[top-1][0],children[top][1],side3s2)+
							 F(children[top-1][0],children[top][2],side3s1)+
							 F(children[top-1][1],children[top][0],side3s2)+
							 F(children[top-1][1],children[top][2],side3s0)+
							 F(children[top-1][2],children[top][0],side3s1)+
							 F(children[top-1][2],children[top][1],side3s0);
					 children[top-1][0] = newSides0;
					 children[top-1][1] = newSides1;
					 children[top-1][2] = newSides2;

				 } else {  // The following case is relevant only for polytomies. 

					 int newSides0 = 0, newSides1 = 0, newSides2 = 0;
					 for (int i = top-1; i >= top + gtb ; i--) {
						 newSides0 += children[i][0]; 
						 newSides1 += children[i][1]; 
						 newSides2 += children[i][2];
					 }
					 children[top][0] = allsides[0] - newSides0;
					 children[top][1] = allsides[1] - newSides1;
					 children[top][2] = allsides[2] - newSides2;

					 for (int i = top; i >= top + gtb; i--) { 
						 if ( children[i][0] == 0)
							 continue;

						 for (int j = top; j >= top + gtb; j--) {                    
							 if (children[j][1] == 0 || i == j)
								 continue;

							 for (int k = top; k >= top + gtb; k--) {
								 if (children[k][2] == 0 || k == i || k == j) 
									 continue;

								 weight += F(children[i][0],children[j][1],children[k][2]);
							 }
						 }
					 }
					 top = top + gtb + 1;
					 children[top-1][0] = newSides0;
					 children[top-1][1] = newSides1;
					 children[top-1][2] = newSides2;
				 } // End of polytomy section
			 }
			 return weight;
		 }

		 public Long calculateWeight() {
			 Long r = dataCollection.geneTreesAsInts != null? 
					 calculateWeightByTraversal():
						 calculateMissingWeight();
					 return r;
		 }
	}

	public void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees) {
	}


	@Override
	public ICalculateWeightTask<Tripartition> getWeightCalculateTask(Tripartition t) {
		return new QuartetWeightTask(t);
	}


	@Override
	protected void prepareWeightTask(ICalculateWeightTask<Tripartition> weigthWork, AbstractComputeMinCostTask<Tripartition> task) {
	}

}
