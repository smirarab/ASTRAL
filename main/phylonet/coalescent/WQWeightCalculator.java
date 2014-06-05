package phylonet.coalescent;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.Iterator;
import java.util.List;

import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;

class WQWeightCalculator extends WeightCalculator<Tripartition> {

	WQInference inference;
	private WQDataCollection dataCollection;

	public WQWeightCalculator(Inference<Tripartition> inference) {
		this.dataCollection = (WQDataCollection) inference.dataCollection;
		this.inference = (WQInference) inference;
	}
	
	class QuartetWeightTask implements CalculateWeightTask<Tripartition>{

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

		private int F(int a,int b,int c) {
			return a*b*c*(a+b+c-3);
		}	
		
		int sharedQuartetCount(Tripartition that, Tripartition other) {
			
			//int [] I = new int [9];
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

		Intersects getSide(int i) {
			if (trip.cluster1.getBitSet().get(i)) {
				return new Intersects(1,0,0);
			} else if (trip.cluster2.getBitSet().get(i)) {
				return new Intersects(0,1,0);
			} else {
				return  new Intersects(0,0,1);
			}
		}
		
		class Intersects {
			int s0;
			int s1;
			int s2;
			
			public Intersects(int s0, int s1, int s2) {
				this.s0 = s0;
				this.s1 = s1;
				this.s2 = s2;
			}
			
		}
		
		long calculateMissingWeight2() {
			long weight = 0;
			Intersects  allsides = null;
			Iterator<STITreeCluster> tit = dataCollection.treeAllClusters.iterator();
			boolean newTree = true;
			
			Deque<Intersects> stack = new ArrayDeque<WQWeightCalculator.QuartetWeightTask.Intersects>();
			for (Integer gtb: dataCollection.geneTreesAsInts){
				if (newTree) {
					STITreeCluster all = tit.next();
					allsides = new Intersects(
						trip.cluster1.getBitSet().intersectionSize(all.getBitSet()),
						trip.cluster2.getBitSet().intersectionSize(all.getBitSet()),
						trip.cluster3.getBitSet().intersectionSize(all.getBitSet()));
					newTree = false;
				}
				if (gtb >= 0){
					stack.push(getSide(gtb));
				} else if (gtb == Integer.MIN_VALUE) {
					stack.clear();
					newTree = true;
				} else {
					for (int i = 0; i > gtb; i--) {
						Intersects side1 = stack.pop();
						Intersects side2 = stack.pop();
						Intersects side = new Intersects(
								side1.s0+side2.s0,
								side1.s1+side2.s1,
								side1.s2+side2.s2);
						stack.push(side);
						Intersects side3 = new Intersects  (
								(allsides.s0 - side.s0 ),
								(allsides.s1 - side.s1 ),
								(allsides.s2 - side.s2 )
						);
						weight += F(side1.s0,side2.s1,side3.s2)+
								F(side1.s0,side2.s2,side3.s1)+
								F(side1.s1,side2.s0,side3.s2)+
								F(side1.s1,side2.s2,side3.s0)+
								F(side1.s2,side2.s0,side3.s1)+
								F(side1.s2,side2.s1,side3.s0);
					}
				}
			}
			return weight;
		}

		public Long calculateWeight() {
			Long r = dataCollection.geneTreesAsInts != null? 
					calculateMissingWeight2():
						calculateMissingWeight();
			weights.put(trip, r);
			if (weights.size() % 100000 == 0)
				System.err.println("Calculated "+weights.size()+" weights");
			return r;
		}
	}

	public void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees) {
	}


	@Override
	public CalculateWeightTask<Tripartition> getWeightCalculateTask(Tripartition t) {
		return new QuartetWeightTask(t);
	}
	
	
	@Override
	protected void prepareWeightTask(CalculateWeightTask<Tripartition> weigthWork, ComputeMinCostTask<Tripartition> task) {
	}

}