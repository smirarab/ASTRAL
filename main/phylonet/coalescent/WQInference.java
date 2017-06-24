package phylonet.coalescent;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Stack;

import phylonet.coalescent.BipartitionWeightCalculator.Quadrapartition;
import phylonet.coalescent.BipartitionWeightCalculator.Results;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;

public class WQInference extends AbstractInference<Tripartition> {
	
	int forceAlg = -1;
	long maxpossible;
	
	public WQInference(Options inOptions, List<Tree> trees, List<Tree> extraTrees) {
		super(inOptions, trees, extraTrees);
		
		this.forceAlg = inOptions.getAlg();
	}


	/**
	 * Calculates maximum possible score, to be used for normalization.
	 * @return
	 */
	long calculateMaxPossible() {
		long weight = 0;
		Integer  allsides = null;
		Iterator<STITreeCluster> tit = ((WQDataCollection)this.dataCollection).treeAllClusters.iterator();
		boolean newTree = true;
		
		Deque<Integer> stack = new ArrayDeque<Integer>();
		// TODO: this should not use private stuff from weight calculator. 
		//       redo to use tree objects. 
		for (Integer gtb: ((WQWeightCalculator)this.weightCalculator).geneTreesAsInts()){
			if (newTree) {
				allsides = tit.next().getBitSet().cardinality();
				newTree = false;
			}
			if (gtb >= 0){
				stack.push(1);
			} else if (gtb == Integer.MIN_VALUE) {
				stack.clear();
				newTree = true;
			}  else {
			    ArrayList<Integer> children = new ArrayList<Integer>();
			    Integer newSide = 0;
			    for (int i = gtb; i < 0 ; i++) {
			    	Integer pop = stack.pop();
			        children.add(pop);
			        newSide+=pop;
			    }
			    stack.push(newSide);
                Integer sideRemaining = allsides - newSide;
                if ( sideRemaining !=0) {
                    children.add(sideRemaining);
                }
                for (int i = 0; i < children.size(); i++) {
                	Long a = children.get(i) + 0l;
                    
                    for (int j = i+1; j < children.size(); j++) {
                    	Long b = children.get(j) + 0l;
                        /*if (children.size() > 5) {
                        	if ((side1.s0+side2.s0 == 0? 1 :0) +
                        			(side1.s1+side2.s1 == 0? 1 :0) + 
                        			(side1.s2+side2.s2 == 0? 1:0) > 1)
                        		continue;
                        }
                        */
                        for (int k = j+1; k < children.size(); k++) {
                        	Long c = children.get(k) + 0l;
                            weight += (a+b+c-3) *a*b*c;
                        }
                    }
                }
			}
		}
		return weight/4l;
	}
	
	void initializeWeightCalculator() {
		((WQWeightCalculator)this.weightCalculator).setupGeneTrees(this);
		if (this.forceAlg == 2) {
			((WQWeightCalculator)this.weightCalculator).useSetWeightsAlgorithm();
		} 

		this.weightCalculator.initializeWeightContainer(
				this.trees.size() *  GlobalMaps.taxonIdentifier.taxonCount() * 2);
	}
	
	/**
	 * Score first computes the quartet scors and the calls
	 * scoreBranches to annotate branches (if needed). 
	 */
	public double scoreSpeciesTreeWithGTLabels(Tree st, boolean initialize) {

		if (initialize) {
			mapNames();
	
			IClusterCollection clusters = newClusterCollection();
	
	
			this.dataCollection = newCounter(clusters);
			weightCalculator = newWeightCalculator();
	
			WQDataCollection wqDataCollection = (WQDataCollection) this.dataCollection;
			wqDataCollection.preProcess(this);
			this.initializeWeightCalculator();			
			this.maxpossible = this.calculateMaxPossible();
			System.err.println("Number of quartet trees in the gene trees: "+this.maxpossible);
	
			//System.err.println(this.maxpossible);
		}
		
		Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
		long sum = 0l;

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
						}
					}					       
				}
			}
		}
		

		System.err.println("Quartet score is: " + sum/4l);
		System.err.println("Normalized quartet score is: "+ (sum/4l+0.)/this.maxpossible);
		//System.out.println(st.toNewickWD());

		if (this.getBranchAnnotation() == 0){
			for (TNode n: st.postTraverse()) {
				((STINode) n).setData(null);
			}
		} else {
			double logscore = this.scoreBranches(st);
			
			if (this.getBranchAnnotation() % 12 == 0) {
				System.err.println("log local posterior: "+logscore);
				return logscore;
			}
		}
		return (sum/4l+0.)/this.maxpossible;
		
	}
	
	
	private boolean skipNode (TNode node) {
		return 	node.isLeaf() || node.isRoot() || node.getChildCount() > 2 || 
				(node.getParent().getChildCount() >3) ||
				(node.getParent().getChildCount() >2 && !node.getParent().isRoot())||
				((node.getParent().isRoot() && node.getParent().getChildCount() == 2));
	}
	
	private class NodeData {
		Double mainfreq, alt1freqs, alt2freqs;
		Long quartcount;
		Integer effn ;
		Quadrapartition [] quads;
		STBipartition[] bipartitions;
		
	}
	
	/**
	 * Annotates the species tree branches with support, branch length, etc. 
	 * @param st
	 * @return
	 */
	private double scoreBranches(Tree st) {

		double ret = 0;
		
		weightCalculator = new BipartitionWeightCalculator(this,((WQWeightCalculator)this.weightCalculator).geneTreesAsInts());
		
		BipartitionWeightCalculator weightCalculator2 = (BipartitionWeightCalculator) weightCalculator;
		WQDataCollection wqDataCollection = (WQDataCollection) this.dataCollection;
		//wqDataCollection.initializeWeightCalculator(this);
		
		/**
		 * Add bitsets to each node for all taxa under it. 
		 * Bitsets are saved in nodes "data" field
		 */
		Path freqWriter = null;
		Path Rscript = null;
		List<String> freqWriterLines = new ArrayList<String>();
		List<String> RscriptLines = new ArrayList<String>();;
		if (this.getBranchAnnotation() == 16) {
			String freqOutputPath = this.options.getFreqOutputPath();
			freqWriter = Paths.get(freqOutputPath + "/freqQuad.csv");
			Rscript = Paths.get(freqOutputPath + "/freqQuadVisualization.R");
		}
		
		Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
		int numNodes = 0;
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
				String ndName = "N" + Integer.toString(numNodes);
				numNodes += 1;
				node.setName(ndName);
			}
		}
		stack = new Stack<STITreeCluster>();
		
		
		/**
		 * For each node,
		 *   1. create three quadripartitoins for the edge above it
		 *   2. score the quadripartition
		 *   3. save the scores in a list for annotations in the next loop
		 */
		Queue<NodeData> nodeDataList = new LinkedList<NodeData>();
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;
			if (node.isLeaf()) {
				stack.push((STITreeCluster) node.getData());
			} else {
				/**
				 * 1. Create quadripartion
				 */
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
				
				/**
				 * 2. Scores all three quadripartitoins
				 */
				Results s = weightCalculator2.getWeight(quad);
				nd.mainfreq = s.qs;
				nd.effn = s.effn;
				
				
				if (nd.effn < 20) {
					if (!GlobalMaps.taxonNameMap.getSpeciesIdMapper().isSingleSP(cluster.getBitSet()))
						System.err.println("You may want to ignore posterior probabilities and other statistics related to the following "
								+ "branch branch because the effective number of genes impacting it is only "+ nd.effn +
							":\n\t" +
							GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTClusterForGeneCluster(cluster));
				}
				
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

					
				if (this.getBranchAnnotation() == 6 || this.getBranchAnnotation() == 16) {
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
		
		/**
		 * Annotate each branch by updating its data field
		 * according to scores and user's annotation preferences. 
		 */
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
			Double effni = nd.effn + 0.0;
			
			if ( Math.abs((f1+f2+f3) - effni) > 0.01 ) {
				//System.err.println("Adjusting effective N from\t" + effni + "\tto\t" + (f1 + f2 + f3) + ". This should only happen as a result of polytomies in gene trees.");
				effni = f1 + f2 + f3;
			}
			
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
					} else if (this.getBranchAnnotation() == 16) {
						node.setData("'[pp1="+df.format(postQ1)+";pp2="+df.format(postQ2)+";pp3="+df.format(postQ3)+"]'");	
						Quadrapartition[] threequads = nd.quads;
						STBipartition[] biparts = nd.bipartitions;
						
						String lineTmp = node.getName() + "\t" + "t1" + "\t" + threequads[0].toString2() + "\t" + 
										 Double.toString(postQ1) + "\t" + Double.toString(f1) +
										 "\t" + Double.toString(effni);
						
						freqWriterLines.add(lineTmp);
						
						lineTmp = node.getName() + "\t" + "t2" + "\t" + threequads[1].toString2() + "\t" + 
								  Double.toString(postQ2) + "\t" + Double.toString(f2) + 
								  "\t" + Double.toString(effni);
						
						freqWriterLines.add(lineTmp);
						
						lineTmp = node.getName() + "\t" + "t3" + "\t" + threequads[2].toString2() + "\t" +
								  Double.toString(postQ3) + "\t" + Double.toString(f3) +
								  "\t" + Double.toString(effni);
						freqWriterLines.add(lineTmp);
					}
				}
				//i++;
			} 
		}
		if (!nodeDataList.isEmpty())
			throw new RuntimeException("Hmm, this shouldn't happen; "+nodeDataList);
		if (this.getBranchAnnotation() == 16) {
			try {
				RscriptLines.add("#!/usr/bin/env Rscript");
				RscriptLines.add("red='#d53e4f';orange='#1d91c0';blue='#41b6c4';colormap = c(red,orange,blue)");
				RscriptLines.add("require(reshape2);require(ggplot2); sizes <- c(7.67,8);");
				RscriptLines.add("dirPath = '"+options.getFreqOutputPath()+"'; filePath = paste(dirPath"
						+ ",'/freqQuadCorrected.csv',sep=''); md<-read.csv(filePath,header=F,sep='\\t'); md$value = md$V5/md$V6;");
				RscriptLines.add("ggplot(data=md)+aes(x=V8,y=value,fill=V8)+"
						+ "geom_bar(stat='identity',color=1,width=0.8,position='dodge')+"
						+ "theme_bw()+theme(axis.text.x=element_text(angle=90))+scale_fill_manual"
						+ "(values=colormap,name='Topology')+geom_hline(yintercept=1/3,size=0.4,linetype=2)+"
						+ "ylab('quartet score')+facet_wrap(~V7,scales='free_x')+scale_x_discrete(breaks=c())+xlab('')");
				RscriptLines.add("pdfFile = paste(dirPath,'/relativeFreq.pdf',sep=''); ggsave(pdfFile,width = sizes[1], height= sizes[2]);");
				RscriptLines.add("ggplot(data=md)+aes(x=V8,y=V4,fill=V8)+"
						+ "geom_bar(stat='identity',color=1,position='dodge')+"
						+ "theme_bw()+theme(axis.text.x=element_text(angle=90))+scale_fill_manual"
						+ "(values=colormap,name='Topology')+geom_hline(yintercept=1/3,size=0.4,linetype=2)+"
						+ "ylab('quartet score')+facet_wrap(~V7,scales='free_x')+scale_x_discrete(breaks=c())+xlab('')");
				RscriptLines.add("pdfFile = paste(dirPath, '/quartetScore.pdf',sep=''); ggsave(pdfFile,width = sizes[1], height= sizes[2]);");
				RscriptLines.add("ggplot(data=md)+aes(x=V8,y=V5,fill=V8)+"
						+ "geom_bar(stat='identity',color=1,position='dodge')+"
						+ "theme_bw()+theme(axis.text.x=element_text(angle=90))+scale_fill_manual"
						+ "(values=colormap,name='Topology')+geom_hline(aes(yintercept=V6/3),size=0.4,linetype=2)+"
						+ "ylab('quartet score')+facet_wrap(~V7,scales='free_x')+scale_x_discrete(breaks=c())+xlab('')");
				RscriptLines.add("pdfFile = paste(dirPath, '/countTopology.pdf',sep=''); ggsave(pdfFile,width = sizes[1], height= sizes[2]);");
				Files.write(freqWriter, freqWriterLines, Charset.forName("UTF-8"));
				Files.write(Rscript, RscriptLines, Charset.forName("UTF-8"));
			} catch (IOException e) {
				throw new RuntimeException("Hmm, the Rscript and frequency of Quadripartition files cannot be created!");
			}
		}
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
			Vertex all, IClusterCollection clusters) {
		return new WQComputeMinCostTask( (WQInference) dlInference, all,  clusters);
	}

	IClusterCollection newClusterCollection() {
		return new WQClusterCollection(GlobalMaps.taxonIdentifier.taxonCount());
	}
	
	WQDataCollection newCounter(IClusterCollection clusters) {
		return new WQDataCollection((WQClusterCollection)clusters, this);
	}



	@Override
	AbstractWeightCalculator<Tripartition> newWeightCalculator() {
		return new WQWeightCalculator(this);
	}


	@Override
	void setupMisc() {
		this.maxpossible = this.calculateMaxPossible();
		System.err.println("Number of quartet trees in the gene trees: " +
				this.maxpossible);
		
	}

	/**
	 * obsolete
	 */
	private void automaticallyDecideAlgorithm(int geneTreeTripartitonCountSize, int k){
		if (this.forceAlg != -1) {
			return;
		}
		if (k <= 0 || geneTreeTripartitonCountSize <= 0) {
			throw new RuntimeException("gene tree tripartition size or k not set properly");
		}
		if (this.forceAlg == -1) {
			this.forceAlg = ( GlobalMaps.taxonIdentifier.taxonCount() <= 32 || (geneTreeTripartitonCountSize < k*6)) ? 2 : 1;
		} else {
			throw new RuntimeException("Algorithm already set");
		}
	}
}
