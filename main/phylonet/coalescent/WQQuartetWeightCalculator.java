package phylonet.coalescent;

import java.util.List;

import phylonet.tree.model.Tree;

class WQQuartetWeightCalculator extends WeightCalculator<Tripartition> {

	private WQDataCollection dataCollection;

	public WQQuartetWeightCalculator(Inference<Tripartition> inference) {
		dataCollection = (WQDataCollection) inference.dataCollection;
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
				weight += this.trip.sharedQuartetCount(dataCollection.finalTripartitions[i]) * dataCollection.finalCounts[i];
			}
			weights.put(trip, weight);
			if (weights.size() % 100000 == 0)
				System.err.println("Calculated "+weights.size()+" weights");
			return weight;
		}

		public Integer calculateWeight() {
			return calculateMissingWeight();
		}
	}

	public void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees) {
	}



	@Override
	public CalculateWeightTask<Tripartition> getWeightCalculateTask(Tripartition t) {
		return new QuartetWeightTask(t);
	}

}