package phylonet.coalescent;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Stack;

import phylonet.coalescent.BipartitionWeightCalculator.Quadrapartition;
import phylonet.coalescent.BipartitionWeightCalculator.Results;
import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;

public class WQInferenceNoCalculations extends AbstractInferenceNoCalculations<Tripartition> {
	
	int forceAlg = -1;
	long maxpossible;
	public WQInferenceNoCalculations(Options inOptions, List<Tree> trees, List<Tree> extraTrees) {
		super(inOptions, trees, extraTrees);
		
		this.forceAlg = inOptions.getAlg();
	}

	
	public WQInferenceNoCalculations(AbstractInference semiDeepCopy) {
		super(semiDeepCopy);
	}


	public double scoreGeneTree(Tree st, boolean initialize) {

		if (initialize) {
			mapNames();
	
			IClusterCollection clusters = newClusterCollection();
	
	
			this.dataCollection = newCounter(clusters);
			weightCalculator = newWeightCalculator();
	
			WQDataCollection wqDataCollection = (WQDataCollection) this.dataCollection;
			wqDataCollection.preProcess(this);
			wqDataCollection.initializeWeightCalculator(this);
			this.maxpossible = wqDataCollection.calculateMaxPossible();
			System.err.println("Number of quartet trees in the gene trees: "+this.maxpossible);
	
			//System.err.println(this.maxpossible);
		}
		
		GlobalMaps.taxonNameMap.getSpeciesIdMapper().stToGt((MutableTree) st);
		
		Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
		long sum = 0l;
		//long maxsum = 0l;
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
					if (this.getBranchAnnotation() % 2 == 0) {
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
		
		double logscore = this.scoreBranches(st);
		
		GlobalMaps.taxonNameMap.getSpeciesIdMapper().gtToSt((MutableTree) st);
		
		if (this.getBranchAnnotation() % 12 == 0) {
			System.err.println("log local posterior: "+logscore);
			return logscore;
		} else {
			return (sum/4l+0.)/this.maxpossible;
		}
	}
	
	
	boolean skipNode (TNode node) {
		return 	node.isLeaf() || node.isRoot() || node.getChildCount() > 2 || 
				(node.getParent().getChildCount() >3) ||
				(node.getParent().getChildCount() >2 && !node.getParent().isRoot())||
				((node.getParent().isRoot() && node.getParent().getChildCount() == 2));
	}
	
	class NodeData {
		Double mainfreq, alt1freqs, alt2freqs;
		Long quartcount;
		Integer effn ;
		Quadrapartition [] quads;
		STBipartition[] bipartitions;
		
	}
	
	public double scoreBranches(Tree st) {

		double ret = 0;
		
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
		
		
		//List<Double> mainfreqs = new ArrayList<Double>();
		//List<Double> alt1freqs = new ArrayList<Double>();
		//List<Double> alt2freqs = new ArrayList<Double>();
		//List<Long> quartcount = new ArrayList<Long>();
		//List<Integer> effn = new ArrayList<Integer>();
		//List< Quadrapartition []> quads = new ArrayList<Quadrapartition[]>();
		//List<STBipartition[]> bipartitions = new ArrayList<STBipartition[]>();
		
		Queue<NodeData> nodeDataList = new LinkedList<NodeData>();
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;
			if (node.isLeaf()) {
				stack.push((STITreeCluster) node.getData());
			} else {
				NodeData nd = new NodeData();
				nodeDataList.add(nd);
				
				STITreeCluster cluster = (STITreeCluster) node.getData();				
				if (skipNode(node) ) {
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
				} /* else if (node.getParent().isRoot() && node.getParent().getChildCount() == 2) {
					continue;
				} */ 
				else {
					remaining = ((STITreeCluster)node.getParent().getData()).complementaryCluster();
				}
				Quadrapartition quad = weightCalculator2.new Quadrapartition
						(c1,  c2, sister, remaining);
				if (this.getBranchAnnotation() == 7){
					if (remaining.getClusterSize() != 0 && sister.getClusterSize() != 0 && c2.getClusterSize() != 0 && c1.getClusterSize() != 0 ){
						System.err.print(c1.toString()+c2.toString()+"|"+sister.toString()+remaining.toString()+"\n");
					}
				}
				Results s = weightCalculator2.getWeight(quad);
				nd.mainfreq = s.qs;
				nd.effn = s.effn;
				//System.err.println(s.effn + " " + quad);
				
				Quadrapartition[] threequads = new Quadrapartition [] {quad, null,null};
				
				quad = weightCalculator2.new Quadrapartition
						(c1, sister, c2, remaining);
				s = weightCalculator2.getWeight(quad);
				nd.alt1freqs=s.qs;
				threequads[1] = quad;
				
				quad = weightCalculator2.new Quadrapartition
						(c1, remaining, c2, sister);
				s = weightCalculator2.getWeight(quad);
				nd.alt2freqs=s.qs;
				threequads[2] = quad;
			
				nd.quartcount= (c1.getClusterSize()+0l)
						* (c2.getClusterSize()+0l)
						* (sister.getClusterSize()+0l)
						* (remaining.getClusterSize()+0l);

					
				if (this.getBranchAnnotation() == 6) {
					STITreeCluster c1plussis = new STITreeCluster();
					c1plussis.setCluster((BitSet) c1.getBitSet().clone());
					c1plussis.getBitSet().or(sister.getBitSet());
					STITreeCluster c1plusrem = new STITreeCluster();
					c1plusrem.setCluster((BitSet) c1.getBitSet().clone());
					c1plusrem.getBitSet().or(remaining.getBitSet());
					
					STBipartition bmain = new STBipartition(cluster, cluster.complementaryCluster());
					STBipartition b2 = new STBipartition(c1plussis, c1plussis.complementaryCluster());
					STBipartition b3 = new STBipartition(c1plusrem, c1plusrem.complementaryCluster());
	
					STBipartition[] biparts = new STBipartition[] {bmain, b2, b3};
					nd.quads = threequads;
					nd.bipartitions = biparts;
				}
				
			}
		}
		//int i = 0;
		NodeData nd = null;
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;
			
			if (!node.isLeaf()) {
				nd = nodeDataList.poll();
			}
			if (skipNode(node)) {
				node.setData(null);
				continue;
			} 
				
			Double f1 = nd.mainfreq;
			Double f2 = nd.alt1freqs;
			Double f3 = nd.alt2freqs;
			Long quarc = nd.quartcount;
			Integer effni = nd.effn;
			//Long sum = p+a1+a2;
			
			Posterior post = new Posterior(
					f1,f2,f3,(double)effni, options.getLambda());
			double bl = post.branchLength();
			
			node.setParentDistance(bl);
			if (this.getBranchAnnotation() == 0){
				node.setData(null);
			} else if (this.getBranchAnnotation() == 1){
				node.setData(df.format((f1+.0)/effni*100));
			} else if (this.getBranchAnnotation() == 10) {
				df.setMaximumFractionDigits(5);
				node.setData(df.format(post.getPvalue()));
			} else {
				double postQ1 = post.getPost();
				ret += Math.log(postQ1);
				
				if (this.getBranchAnnotation() == 3 || this.getBranchAnnotation() == 12) {
					node.setData(df.format(postQ1));
				} else if (this.getBranchAnnotation() % 2 == 0) {
					post = new Posterior(f2,f1,f3,(double)effni, options.getLambda());
					double postQ2 = post.getPost();
					post =  new Posterior(f3,f1,f2,(double)effni, options.getLambda());
					double postQ3 = post.getPost();
					
					if (this.getBranchAnnotation() == 2)
						node.setData(
								"'[q1="+(f1)/effni+";q2="+(f2)/effni+";q3="+(f3)/effni+
								 ";f1="+f1+";f2="+f2+";f3="+f3+
								 ";pp1="+postQ1+";pp2="+postQ2+";pp3="+postQ3+
								 ";QC="+quarc+";EN="+effni+"]'");
					else if (this.getBranchAnnotation() == 4) {
						node.setData("'[pp1="+df.format(postQ1)+";pp2="+df.format(postQ2)+";pp3="+df.format(postQ3)+"]'");
					} else if (this.getBranchAnnotation() == 6){
						node.setData(df.format(postQ1));
						Quadrapartition[] threequads = nd.quads;
						STBipartition[] biparts = nd.bipartitions;
						System.err.println(threequads[0] +
								" [" + biparts[0].toString2() +"] : "+postQ1 +" ** f1 = "+f1+
								" f2 = "+f2+" f3 = "+f3+" EN = "+ effni+" **");
						System.err.println(threequads[1] +
								" ["+biparts[1].toString2()+"] : "+postQ2+ " ** f1 = "+f2+
								" f2 = "+f1+" f3 = "+f3+" EN = "+ effni+" **");
						System.err.println(threequads[2] +
								" ["+biparts[2].toString2()+"] : "+postQ3+ " ** f1 = "+f3+
								" f2 = "+f1+" f3 = "+f2+" EN = "+ effni+" **");
					}  else if (this.getBranchAnnotation() == 8){
						node.setData(
								"'[q1="+df.format((f1)/effni)+
								 ";q2="+df.format((f2)/effni)+
								 ";q3="+df.format((f3)/effni)+"]'");
					}
				}
				//i++;
			} 
		}
		if (!nodeDataList.isEmpty())
			throw new RuntimeException("Hmm, this shouldn't happen; "+nodeDataList);
		
		return ret;
	}


	@Override
	Long getTotalCost(Vertex all) {
		System.err.println("Normalized score (portion of input quartet trees satisfied): " + 
				all._max_score/4./this.maxpossible);
		return (long) (all._max_score/4l);
	}


	@Override
	AbstractComputeMinCostTask<Tripartition> newComputeMinCostTask(AbstractInference<Tripartition> dlInference,
			Vertex all, IClusterCollection clusters, boolean isWriteToQueue) {
		return new WQComputeMinCostTask( (WQInferenceNoCalculations) dlInference, all,  clusters, isWriteToQueue);
	}

	IClusterCollection newClusterCollection() {
		return new WQClusterCollection(GlobalMaps.taxonIdentifier.taxonCount());
	}
	
	WQDataCollection newCounter(IClusterCollection clusters) {
		return new WQDataCollection((WQClusterCollection)clusters, this.forceAlg, this);
	}



	@Override
	AbstractWeightCalculator<Tripartition> newWeightCalculator() {
		return new WQWeightCalculator(this, super.queue2);
	}
	AbstractWeightCalculatorNoCalculations<Tripartition> newWeightCalculatorNoCalculations() {
		return new WQWeightCalculatorNoCalculations(this, super.queue1);
	}

}
