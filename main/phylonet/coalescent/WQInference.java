package phylonet.coalescent;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;

import phylonet.coalescent.BipartitionWeightCalculator.Quadrapartition;
import phylonet.coalescent.BipartitionWeightCalculator.Results;
import phylonet.coalescent.WQWeightCalculator.QuartetWeightTask.Intersects;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;

public class WQInference extends AbstractInference<Tripartition> {
	
	int forceAlg = -1;
	long maxpossible;
	
	public WQInference(boolean rooted, boolean extrarooted, List<Tree> trees,
			List<Tree> extraTrees, boolean exactSolution, boolean duploss, int alg, int addExtra,
			boolean outputCompletedGenes, boolean outSearch, boolean run) {
		super(rooted, extrarooted, trees, extraTrees, exactSolution, 
				addExtra, outputCompletedGenes, outSearch, run);
		this.forceAlg = alg;
	}

	
	public void scoreGeneTree(Tree st) {

		mapNames();

		IClusterCollection clusters = newClusterCollection();


		this.dataCollection = newCounter(clusters);
		weightCalculator = newWeightCalculator();

		WQDataCollection wqDataCollection = (WQDataCollection) this.dataCollection;
		wqDataCollection.preProcess(this);
		wqDataCollection.initializeWeightCalculator(this);
		this.maxpossible = wqDataCollection.calculateMaxPossible();
		//System.err.println("Number of quartet trees in the gene trees: "+this.maxpossible);

		System.err.println(this.maxpossible);
		
		Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
		long sum = 0l;
		long maxsum = 0l;
		for (TNode node: st.postTraverse()) {
			if (node.isLeaf()) {
				String nodeName = node.getName(); //GlobalMaps.TaxonNameMap.getSpeciesName(node.getName());

				STITreeCluster cluster = new STITreeCluster();
				Integer taxonID = GlobalMaps.taxonIdentifier.taxonId(nodeName);
				cluster.addLeaf(taxonID);

				stack.add(cluster);

			} else {
				ArrayList<STITreeCluster> childbslist = new ArrayList<STITreeCluster>();
				BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
				for (TNode child: node.getChildren()) {
					STITreeCluster pop = stack.pop();
					childbslist.add(pop);
					bs.or(pop.getBitSet());
				}
				
				STITreeCluster cluster = new STITreeCluster();
				cluster.setCluster((BitSet) bs.clone());

				//((STINode)node).setData(new GeneTreeBitset(node.isRoot()? -2: -1));
				stack.add(cluster);


				STITreeCluster remaining = cluster.complementaryCluster();
				if (remaining.getClusterSize() != 0) {
					childbslist.add(remaining);
				}
				for (int i = 0; i < childbslist.size(); i++) {
					for (int j = i+1; j < childbslist.size(); j++) {
						for (int k = j+1; k < childbslist.size(); k++) {
							Tripartition trip = new Tripartition(childbslist.get(i),  childbslist.get(j), childbslist.get(k));
							Long s = weightCalculator.getWeight(trip, null);
							sum += s;
							//Long m = this.dataCollection.maxPossibleScore(trip);
							//((STINode)node).setData(s*100l/m);
							//maxsum += m;
						}
					}					       
				}
				if (childbslist.size() > 3) {
					for (STITreeCluster chid :childbslist) {
						System.err.print(chid.getClusterSize()+" ");
					}
					System.err.println(" (polytomy)");
				}
			}
		}
/*		if (4l*this.maxpossible != maxsum) {
			throw new RuntimeException("Hmm... "+maxsum+" "+4l*this.maxpossible);
		}*/
		System.err.println("Quartet score is: " + sum/4l);
		System.err.println("Normalized quartet score is: "+ (sum/4l+0.)/this.maxpossible);
		//System.out.println(st.toNewickWD());
		this.scoreGeneTree2(st);
	}

	
	public void scoreGeneTree2(Tree st) {

		weightCalculator = new BipartitionWeightCalculator(this);
		BipartitionWeightCalculator weightCalculator2 = (BipartitionWeightCalculator) weightCalculator;
		WQDataCollection wqDataCollection = (WQDataCollection) this.dataCollection;
		//wqDataCollection.initializeWeightCalculator(this);
		
		Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;
			if (node.isLeaf()) {
				String nodeName = node.getName(); //GlobalMaps.TaxonNameMap.getSpeciesName(node.getName());

				STITreeCluster cluster = new STITreeCluster();
				Integer taxonID = GlobalMaps.taxonIdentifier.taxonId(nodeName);
				cluster.addLeaf(taxonID);

				stack.add(cluster);
				node.setData(cluster);

			} else {
				ArrayList<STITreeCluster> childbslist = new ArrayList<STITreeCluster>();
				BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
				for (TNode child: n.getChildren()) {
					STITreeCluster pop = stack.pop();
					childbslist.add(pop);
					bs.or(pop.getBitSet());
				}
				
				STITreeCluster cluster = new STITreeCluster();
				cluster.setCluster((BitSet) bs.clone());

				//((STINode)node).setData(new GeneTreeBitset(node.isRoot()? -2: -1));
				stack.add(cluster);
				node.setData(cluster);
			}
		}
		stack = new Stack<STITreeCluster>();
		List<Results> mainfreqs = new ArrayList<Results>();
		List<Results> alt1freqs = new ArrayList<Results>();
		List<Results> alt2freqs = new ArrayList<Results>();
		List<Long> quartcount = new ArrayList<Long>();
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;
			if (node.isLeaf()) {
				stack.push((STITreeCluster) node.getData());
			} else {
				STITreeCluster cluster = (STITreeCluster) node.getData();
				STITreeCluster c1 = stack.pop();
				STITreeCluster c2 = stack.pop();
				stack.push(cluster);
				
				if (node.isRoot() ) {
					continue;
				}
				
				if (node.getChildCount() > 2) {
					throw new RuntimeException("Polytomies not implemented yet. ");
				}
				STITreeCluster sister;
				STITreeCluster remaining;
				Iterator<STINode> pcit = node.getParent().getChildren().iterator();
				STINode pc = pcit.next();
				if ( pc == n ) pc = pcit.next(); 
				sister = (STITreeCluster)pc.getData();
				if (node.getParent().isRoot() && node.getParent().getChildCount() == 3) {
					pc = pcit.next();
					if (pc == n) pc = pcit.next(); 
					remaining = (STITreeCluster)pc.getData();;					
				} else {
					remaining = ((STITreeCluster)node.getParent().getData()).complementaryCluster();
				}
				Quadrapartition quad = weightCalculator2.new Quadrapartition
						(c1,  c2, sister, remaining);
				Results s = weightCalculator2.getWeight(quad);
				mainfreqs.add(s);	
				
				quad = weightCalculator2.new Quadrapartition
						(c1, sister, c2, remaining);
				s = weightCalculator2.getWeight(quad);
				alt1freqs.add(s);
				
				quad = weightCalculator2.new Quadrapartition
						(c1, remaining, c2, sister);
				s = weightCalculator2.getWeight(quad);
				alt2freqs.add(s);
				
				quartcount.add( (c1.getClusterSize()+0l)
						* (c2.getClusterSize()+0l)
						* (sister.getClusterSize()+0l)
						* (remaining.getClusterSize()+0l));
			}
		}
		int i = 0;
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;
			if (node.isLeaf() || node.isRoot() ||
					(node.getParent().isRoot() && node.getParent().getChildCount() ==2)) {
				node.setData(null);
			} else{
				Results p = mainfreqs.get(i);
				Results a1 = alt1freqs.get(i);
				Results a2 = alt2freqs.get(i);
				Long quarc = quartcount.get(i);
				if (this.getBranchAnnotation() == 2) {
					node.setData("[q1="+p.freq+";q2="+a1.freq+";q3="+a2.freq+
							";f1="+p.succ+";f2="+a1.succ+";f3="+a2.succ+
							";QC="+quarc+"]");
				} else if (this.getBranchAnnotation() == 1){
					node.setData(p.freq);
				}
				i++;
				Double bl = -Math.log(1.5*(1.0-p.freq/100.0));
				if (bl.isInfinite()) {
					bl = 10.;
				}
				node.setParentDistance(bl);
			} 
		}
		
	}


	@Override
	Long getTotalCost(Vertex all) {
		System.err.println("Normalized score (portion of input quartet trees satisfied): " + 
				all._max_score/4./this.maxpossible);
		return (long) (all._max_score/4l);
	}


	@Override
	AbstractComputeMinCostTask<Tripartition> newComputeMinCostTask(AbstractInference<Tripartition> dlInference,
			Vertex all, IClusterCollection clusters) {
		return new WQComputeMinCostTask( (WQInference) dlInference, all,  clusters);
	}

	IClusterCollection newClusterCollection() {
		return new WQClusterCollection(GlobalMaps.taxonIdentifier.taxonCount());
	}
	
	WQDataCollection newCounter(IClusterCollection clusters) {
		return new WQDataCollection((WQClusterCollection)clusters, this.forceAlg, this);
	}



	@Override
	AbstractWeightCalculator<Tripartition> newWeightCalculator() {
		return new WQWeightCalculator(this);
	}

}
