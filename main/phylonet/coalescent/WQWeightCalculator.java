package phylonet.coalescent;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
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

		int calculateMissingWeight() {
			// System.err.print("Calculating weight for: " + biggerSTB);
			int weight = 0;
			for (int i=0; i < dataCollection.finalCounts.length; i++) {
				weight += sharedQuartetCount(trip,dataCollection.finalTripartitions[i]) * dataCollection.finalCounts[i];
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

		Intersects subtractSide(Intersects allsides, Intersects side) {
			return new Intersects  (
					(allsides.s0 - side.s0 ),
					(allsides.s1 - side.s1 ),
					(allsides.s2 - side.s2 )
			);
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
		
		
		int countAll(Intersects side1, Intersects side2, Intersects side3){
			return
					F(side1.s0,side2.s1,side3.s2)+
					F(side1.s0,side2.s2,side3.s1)+
					F(side1.s1,side2.s0,side3.s2)+
					F(side1.s1,side2.s2,side3.s0)+
					F(side1.s2,side2.s0,side3.s1)+
					F(side1.s2,side2.s1,side3.s0);

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
		
		int calculateMissingWeight2() {
			// System.err.print("Calculating weight for: " + biggerSTB);
			int weight = 0;
			//Map<TNode,Side > nodeData = new HashMap<TNode, Side>();
			Iterator<STITreeCluster> tit = dataCollection.treeAllClusters.iterator();
			Stack<Intersects> stack = new Stack<WQWeightCalculator.QuartetWeightTask.Intersects>();
			
			for (Tree tree : inference.trees){
				stack.clear();
				STITreeCluster all = tit.next();
				Intersects  allsides = new Intersects(
					trip.cluster1.getBitSet().intersectionSize(all.getBitSet()),
					trip.cluster2.getBitSet().intersectionSize(all.getBitSet()),
					trip.cluster3.getBitSet().intersectionSize(all.getBitSet()));
				
				for (TNode tn: tree.postTraverse()) {
					if (tn.isLeaf()) {
						stack.push(getSide((Integer) ((STINode) tn).getData()));
					} else {
						Intersects side1 = stack.pop();
						Intersects side2 = stack.pop();
						Intersects side = new Intersects(
								side1.s0+side2.s0,
								side1.s1+side2.s1,
								side1.s2+side2.s2);
						stack.push(side);
						Intersects side3 = subtractSide(allsides,side);
						weight += countAll(side1, side2, side3);
					}
				}
			}
			return weight;
		}

		public Integer calculateWeight() {
			int r = calculateMissingWeight2();
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