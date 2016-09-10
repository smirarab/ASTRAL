package phylonet.coalescent;

import java.util.Iterator;
import java.util.List;

import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;

class WQWeightCalculator extends AbstractWeightCalculator<Tripartition> {

	WQInference inference;
	private WQDataCollection dataCollection;
	private WeightCalculatorAlgorithm algorithm;

	public WQWeightCalculator(AbstractInference<Tripartition> inference) {
		super(false);
		this.dataCollection = (WQDataCollection) inference.dataCollection;
		this.inference = (WQInference) inference;
		this.algorithm = new TraversalWeightCalculator();
	}

	abstract class WeightCalculatorAlgorithm {
		long F(int a, int b, int c) {
			if (a < 0 || b < 0 || c < 0) {
				throw new RuntimeException(
					"negative side not expected: " + a + " " + b + " " + c);
			}
			long ret = (a + b + c - 3);
			ret *= a * b * c;
			return ret;
		}

		abstract Long calculateWeight(Tripartition trip);
	}

	@Override
	Long calculateWeight(Tripartition t,
			AbstractComputeMinCostTask<Tripartition> minCostTask) {
		return this.algorithm.calculateWeight(t);
	}

	class TraversalWeightCalculator extends WeightCalculatorAlgorithm {
		
		int[][] stack = new int[GlobalMaps.taxonIdentifier.taxonCount() + 2][3];

		int[][] overlap = new int[3][GlobalMaps.taxonIdentifier.taxonCount() +1];
		int[][] overlapind = new int[3][GlobalMaps.taxonIdentifier.taxonCount() +1];

		
		Long calculateWeight(Tripartition trip) {

			long weight = 0;
			int[] allsides = null;
			Iterator<STITreeCluster> tit = dataCollection.treeAllClusters.iterator();
			boolean newTree = true;
			int top = 0; // The first empty place on stack (generally)
			for (Integer gtb : dataCollection.geneTreesAsInts) {
				if (newTree) {
					STITreeCluster all = tit.next();
					allsides = new int[] {
							trip.cluster1.getBitSet().intersectionSize(all.getBitSet()),
							trip.cluster2.getBitSet().intersectionSize(all.getBitSet()),
							trip.cluster3.getBitSet().intersectionSize(all.getBitSet())};
					newTree = false;
				}
				if (gtb >= 0) { // Leaf nodes
					if (trip.cluster1.getBitSet().get(gtb)) {
						stack[top][0] = 1;
						stack[top][1] = 0;
						stack[top][2] = 0;
					} else if (trip.cluster2.getBitSet().get(gtb)) {
						stack[top][0] = 0;
						stack[top][1] = 1;
						stack[top][2] = 0;
					} else if (trip.cluster3.getBitSet().get(gtb)) {
						stack[top][0] = 0;
						stack[top][1] = 0;
						stack[top][2] = 1;
					} else { // This can happen due to missing data
						stack[top][0] = 0;
						stack[top][1] = 0;
						stack[top][2] = 0;
					}
					top++;
				} else if (gtb == Integer.MIN_VALUE) { // delimiter between trees
					top = 0;
					newTree = true;
				} else if (gtb == -2) { // Internal nodes

					top--;
					int newSides0 = stack[top][0] + stack[top - 1][0];
					int newSides1 = stack[top][1] + stack[top - 1][1];
					int newSides2 = stack[top][2] + stack[top - 1][2];

					int side3s0 = allsides[0] - newSides0;
					int side3s1 = allsides[1] - newSides1;
					int side3s2 = allsides[2] - newSides2;

					weight += F(stack[top][0], stack[top - 1][1], side3s2)
							+ F(stack[top][0], stack[top - 1][2], side3s1)
							+ F(stack[top][1], stack[top - 1][0], side3s2)
							+ F(stack[top][1], stack[top - 1][2], side3s0)
							+ F(stack[top][2], stack[top - 1][0], side3s1)
							+ F(stack[top][2], stack[top - 1][1], side3s0);

					stack[top - 1][0] = newSides0;
					stack[top - 1][1] = newSides1;
					stack[top - 1][2] = newSides2;
				} else { // The following case is relevant only for polytomies.

					int [] nzc = {0,0,0};
					int [] newSides = {0,0,0};
					for (int side = 0; side < 3; side++) {
						for (int i = top - 1; i >= top + gtb; i--) {
							if (stack[i][side] > 0) {
								newSides[side] += stack[i][side];
								overlap[side][nzc[side]] = stack[i][side]; 
								overlapind[side][nzc[side]++] = i;
							}
						}					
						stack[top][side] = allsides[side] - newSides[side];

						if (stack[top][side] > 0) {
							overlap[side][nzc[side]] = stack[top][side]; 
							overlapind[side][nzc[side]++] = top;
						}
					}

					for (int i = nzc[0] - 1; i >= 0; i--) {
						for (int j = nzc[1] - 1; j >= 0; j--) {
							for (int k = nzc[2] - 1; k >= 0; k--) {
								if ( (overlapind[0][i] != overlapind[1][j]) &&
										(overlapind[0][i] != overlapind[2][k]) &&
										(overlapind[1][j] != overlapind[2][k]))
									weight += F(overlap[0][i], overlap[1][j], overlap[2][k]);
							}
						}
					}
					top = top + gtb + 1;

					stack[top - 1][0] = newSides[0];
					stack[top - 1][1] = newSides[1];
					stack[top - 1][2] = newSides[2];
				} // End of polytomy section

			}

			return weight;
		}
	}

	class SetWeightCalculator extends WeightCalculatorAlgorithm {
		Long calculateWeight(Tripartition trip) {
			long weight = 0l;
			for (int i = 0; i < dataCollection.finalCounts.length; i++) {
				weight += sharedQuartetCount(trip,
						dataCollection.finalTripartitions[i])
						* dataCollection.finalCounts[i];
			}
			return weight;
		}

		long sharedQuartetCount(Tripartition that, Tripartition other) {

			int I0 = that.cluster1.getBitSet().intersectionSize(
					other.cluster1.getBitSet()), I1 = that.cluster1.getBitSet()
					.intersectionSize(other.cluster2.getBitSet()), I2 = that.cluster1
					.getBitSet().intersectionSize(other.cluster3.getBitSet()), I3 = that.cluster2
					.getBitSet().intersectionSize(other.cluster1.getBitSet()), I4 = that.cluster2
					.getBitSet().intersectionSize(other.cluster2.getBitSet()), I5 = that.cluster2
					.getBitSet().intersectionSize(other.cluster3.getBitSet()), I6 = that.cluster3
					.getBitSet().intersectionSize(other.cluster1.getBitSet()), I7 = that.cluster3
					.getBitSet().intersectionSize(other.cluster2.getBitSet()), I8 = that.cluster3
					.getBitSet().intersectionSize(other.cluster3.getBitSet());

			return F(I0, I4, I8) + F(I0, I5, I7) + F(I1, I3, I8)
					+ F(I1, I5, I6) + F(I2, I3, I7) + F(I2, I4, I6);
		}
	}
	
	public void useSetWeights() {
		algorithm = new SetWeightCalculator();
	}

	public void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees) {
	}


}
