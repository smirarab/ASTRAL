package phylonet.coalescent;
import phylonet.coalescent.BipartitionWeightCalculator.Results;
import phylonet.coalescent.Posterior;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;

import phylonet.coalescent.BipartitionWeightCalculator.Quadrapartition;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
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
				if (childbslist.size() > 3) {
					for (STITreeCluster chid :childbslist) {
						System.err.print(chid.getClusterSize()+" ");
					}
					System.err.println(" (polytomy)");
					if (this.getBranchAnnotation() == 4) {
						continue;
					}
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
			}
		}
/*		if (4l*this.maxpossible != maxsum) {
			throw new RuntimeException("Hmm... "+maxsum+" "+4l*this.maxpossible);
		}*/
		System.err.println("Quartet score is: " + sum/4l);
		System.err.println("Normalized quartet score is: "+ (sum/4l+0.)/this.maxpossible);
		//System.out.println(st.toNewickWD());
		this.scoreBranches(st);
	}

	
	public void scoreBranches(Tree st) {

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
		List<Long> mainfreqs = new ArrayList<Long>();
		List<Long> alt1freqs = new ArrayList<Long>();
		List<Long> alt2freqs = new ArrayList<Long>();
		List<Long> quartcount = new ArrayList<Long>();
		List<Integer> effn = new ArrayList<Integer>();
		
		
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;
			if (node.isLeaf()) {
				stack.push((STITreeCluster) node.getData());
			} else {
				STITreeCluster cluster = (STITreeCluster) node.getData();				
				if (node.isRoot() || node.getChildCount() > 2 || (node.getParent().getChildCount() >3) ||
						(node.getParent().getChildCount() >2 && !node.getParent().isRoot()) ) {
					for (int i =0; i< node.getChildCount(); i++) {
						stack.pop();
					}
					stack.push(cluster);
					continue;
				}
				
				STITreeCluster c1 = stack.pop();
				STITreeCluster c2 = stack.pop();
				stack.push(cluster);
				
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
				mainfreqs.add(s.qs);
				effn.add(s.effn);
				//System.err.println(s.effn + " " + quad);
				
				quad = weightCalculator2.new Quadrapartition
						(c1, sister, c2, remaining);
				s = weightCalculator2.getWeight(quad);
				alt1freqs.add(s.qs);
				
				quad = weightCalculator2.new Quadrapartition
						(c1, remaining, c2, sister);
				s = weightCalculator2.getWeight(quad);
				alt2freqs.add(s.qs);
				
				quartcount.add( (c1.getClusterSize()+0l)
						* (c2.getClusterSize()+0l)
						* (sister.getClusterSize()+0l)
						* (remaining.getClusterSize()+0l));
			}
		}
		int i = 0;
		
		DecimalFormat df = new DecimalFormat();
		df.setMaximumFractionDigits(2);
		
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;
			if (node.isLeaf() || node.isRoot() ||
					(node.getParent().isRoot() && node.getParent().getChildCount() ==2)) {
				node.setData(null);
			} else{
				if (node.isRoot() || node.getChildCount() > 2 || (node.getParent().getChildCount() >3) ||
						(node.getParent().getChildCount() >2 && !node.getParent().isRoot()) ) {
					node.setData(null);
					continue;
				}
				Long p = mainfreqs.get(i);
				Long a1 = alt1freqs.get(i);
				Long a2 = alt2freqs.get(i);
				Long quarc = quartcount.get(i);
				Integer effni = effn.get(i);
				Long sum = p+a1+a2;
				
				Posterior bl_tmp =new Posterior((double)p,(double)a1,(double) a2,(double)effni);
				double bl = bl_tmp.branchLength();
				
				node.setParentDistance(bl);
				if (this.getBranchAnnotation() == 0){
					node.setData(null);
				} else if (this.getBranchAnnotation() == 1){
					node.setData(df.format((p+.0)/sum*100));
				} else {
					Posterior pst_tmp = new Posterior((double)p,(double)a1,(double)a2,(double)effni);
					double post_m = pst_tmp.getPost();
					
					if (this.getBranchAnnotation() == 3) {
						node.setData(df.format(post_m));
					} else if (this.getBranchAnnotation() % 2 == 0) {
						pst_tmp = new Posterior((double)a1,(double)p,(double)a2,(double)effni);
						double post_a1 = pst_tmp.getPost();
						//pst_tmp =  new Posterior((double)a2,(double)p,(double)a1,(double)numTrees);
						double post_a2 = 1.0 - post_a1 - post_m;
						
						if (this.getBranchAnnotation() == 2)
							node.setData("[q1="+(p+.0)/sum+";q2="+(a1+.0)/sum+";q3="+(a2+.0)/sum+
									 ";f1="+p+";f2="+a1+";f3="+a2+
									 ";pp1="+post_m+";pp2="+post_a1+";pp3="+post_a2+
									 ";QC="+quarc+";EN="+effni+"]");
						else
							node.setData("'[pp1="+df.format(post_m)+";pp2="+df.format(post_a1)+";pp3="+df.format(post_a2)+"]'");
					} 
				}
				i++;
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
