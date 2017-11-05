package phylonet.coalescent;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
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
		if (weightCalculator instanceof WQWeightCalculator
				&& ((WQWeightCalculator)weightCalculator).algorithm instanceof WQWeightCalculator.CondensedTraversalWeightCalculator){
			return ((WQWeightCalculator.CondensedTraversalWeightCalculator)
					((WQWeightCalculator)weightCalculator).algorithm).polytree.maxScore / 4L
					- unresolvableQuartets();
		}

		//TODO: MUTIND: In the multi individual case, some quartets can never be satisfied. 
		//      We should compute their number and substract that from maxpossible here. 
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
        return weight/4l - unresolvableQuartets();
	}


	private long unresolvableQuartets() {
		long ret = 0;
		for (STITreeCluster gtCL : ((WQDataCollection)this.dataCollection).treeAllClusters) {
			int[] counts = new int [GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSpeciesCount()];
			int size = gtCL.getClusterSize();
			BitSet bs = gtCL.getBitSet();
	        for (int i = bs.nextSetBit(0); i >=0 ; i = bs.nextSetBit(i+1)) {
	            counts[(GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSpeciesIdForTaxon(i))]++;
	        }
	        for (long count : counts) {
	        	ret += (count*(count-1)*(count-2))/6l*(size-count)+
	        			(count*(count-1)*(count-2)*(count-3))/24l;
	        }
		}
		return ret;
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
	 * This method first computes the quartet scores and then calls
	 * scoreBranches to annotate branches (if needed).
	 * The method assumes the input tree st has labels of individuals (not species). 
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
			//ASTRAL IV SPECIFIC
			this.maxpossible = this.calculateMaxPossible();
			System.err.println("Number of quartet trees in the gene trees: "+this.maxpossible);

			//System.err.println(this.maxpossible);
		}

		Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
		long sum = 0l;

		for (TNode node: st.postTraverse()) {
			if (node.isLeaf()) {
				String nodeName = node.getName(); //GlobalMaps.TaxonNameMap.getSpeciesName(node.getName());

				STITreeCluster cluster = GlobalMaps.taxonIdentifier.newCluster();
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

				STITreeCluster cluster = GlobalMaps.taxonIdentifier.newCluster();
				cluster.setCluster((BitSet) bs.clone());

				//((STINode)node).setData(new GeneTreeBitset(node.isRoot()? -2: -1));
				stack.add(cluster);


				STITreeCluster remaining = cluster.complementaryCluster();
				if (remaining.getClusterSize() != 0) {
					childbslist.add(remaining);
				}
				if (childbslist.size() > 3) {
					/*for (STITreeCluster chid :childbslist) {
						System.err.print(chid.getClusterSize()+" ");
					}
					System.err.println(" (polytomy)");*/
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


		System.err.println("Final quartet score is: " + sum/4l);
		System.err.println("Final normalized quartet score is: "+ (sum/4l+0.)/this.maxpossible);
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
		BufferedWriter freqWriter = null;
		BufferedWriter Rscript = null;
		//List<String> freqWriterLines = new ArrayList<String>();
		if (this.getBranchAnnotation() == 16) {
			String freqOutputPath = this.options.getFreqOutputPath();
			try {
				Rscript = new BufferedWriter(new FileWriter(freqOutputPath + File.separator+ "freqQuadVisualization.R"));
				freqWriter = new BufferedWriter(new FileWriter(freqOutputPath + File.separator+ "freqQuad.csv"));
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}

		Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
		int numNodes = 0;
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;
			if (node.isLeaf()) {
				String nodeName = node.getName(); //GlobalMaps.TaxonNameMap.getSpeciesName(node.getName());

				STITreeCluster cluster = GlobalMaps.taxonIdentifier.newCluster();
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

				STITreeCluster cluster = GlobalMaps.taxonIdentifier.newCluster();
				cluster.setCluster((BitSet) bs.clone());

				//((STINode)node).setData(new GeneTreeBitset(node.isRoot()? -2: -1));
				stack.add(cluster);
				node.setData(cluster);
				if (options.getBranchannotation() == 16) {
					String ndName = "N" + Integer.toString(numNodes);
					numNodes += 1;
					node.setName(ndName);
				}
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
				 * 1. Create quadripartition
				 */
				NodeData nd = null;

				STITreeCluster cluster = (STITreeCluster) node.getData();		
				STITreeCluster c1 = null, c2 = null;
				long cs = cluster.getClusterSize()+0l;
				
				for (int i =0; i< node.getChildCount(); i++) {
					if (c1 == null)
						c1 = stack.pop();
					else if (c2 == null)
						c2 = stack.pop();
					else
						stack.pop();
				}
				stack.push(cluster);

				
				if (cs > 1 && GlobalMaps.taxonNameMap.getSpeciesIdMapper().isSingleSP(cluster.getBitSet()))
				{
					STITreeCluster[] sisterRemaining = getSisterRemaining(node);
					STITreeCluster sister = sisterRemaining[0]; 
					STITreeCluster remaining = sisterRemaining[1];
					
					nd = new NodeData();
					nodeDataList.add(nd);
					nd.mainfreq = 0d;
					nd.alt1freqs = 0d;
					nd.alt2freqs = 0d;
					nd.effn = 0;
					
					BitSet bitSet = cluster.getBitSet();			
					for (int j = bitSet.nextSetBit(0); j >= 0; j = bitSet.nextSetBit(j + 1)) {
						c1 = new STITreeCluster(cluster);
						c1.getBitSet().clear(j);
						c2 = GlobalMaps.taxonIdentifier.newCluster();
						c2.getBitSet().set(j);
						Quadrapartition[] threequads = new Quadrapartition [] { 
								weightCalculator2.new Quadrapartition (c1,  c2, sister, remaining), 
								weightCalculator2.new Quadrapartition (c1, sister, c2, remaining),
								weightCalculator2.new Quadrapartition (c1, remaining, c2, sister)
						};
						
						Results s = weightCalculator2.getWeight(threequads);

						nd.mainfreq += s.qs[0];
						nd.alt1freqs += s.qs[1];
						nd.alt2freqs += s.qs[2];
						nd.effn += s.effn;
					}
					
					nd.mainfreq /= cs;
					nd.alt1freqs /= cs;
					nd.alt2freqs /= cs;
					nd.effn /= (int) cs;
					
					nd.quartcount =  (cs*(cs-1)/2)
							* (sister.getClusterSize()+0l)
							* (remaining.getClusterSize()+0l);
					
					
				} else if (! skipNode(node) ) {
					STITreeCluster[] sisterRemaining = getSisterRemaining(node);
					STITreeCluster sister = sisterRemaining[0]; 
					STITreeCluster remaining = sisterRemaining[1];
					
					Quadrapartition[] threequads = new Quadrapartition [] { 
							weightCalculator2.new Quadrapartition (c1,  c2, sister, remaining), 
							weightCalculator2.new Quadrapartition (c1, sister, c2, remaining),
							weightCalculator2.new Quadrapartition (c1, remaining, c2, sister)
					};

					/**
					 * 2. Scores all three quadripartitoins
					 */
					Results s = weightCalculator2.getWeight(threequads);
					nd = new NodeData();
					nodeDataList.add(nd);
					nd.mainfreq = s.qs[0];
					nd.alt1freqs=s.qs[1];
					nd.alt2freqs=s.qs[2];
					nd.effn = s.effn;

					nd.quartcount= (c1.getClusterSize()+0l)
							* (c2.getClusterSize()+0l)
							* (sister.getClusterSize()+0l)
							* (remaining.getClusterSize()+0l);


					if (this.getBranchAnnotation() == 7){
						if (remaining.getClusterSize() != 0 && sister.getClusterSize() != 0 && c2.getClusterSize() != 0 && c1.getClusterSize() != 0 ){
							System.err.print(c1.toString()+c2.toString()+"|"+sister.toString()+remaining.toString()+"\n");
						}
					}
					if (this.getBranchAnnotation() == 6 || this.getBranchAnnotation() == 16) {
						STITreeCluster c1plussis = GlobalMaps.taxonIdentifier.newCluster();

						c1plussis.setCluster((BitSet) c1.getBitSet().clone());
						c1plussis.getBitSet().or(sister.getBitSet());
						STITreeCluster c1plusrem = GlobalMaps.taxonIdentifier.newCluster();
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
				
				if (nd != null && nd.effn < 20) {
					System.err.println("You may want to ignore posterior probabilities and other statistics related to the following "
							+ "branch branch because the effective number of genes impacting it is only "+ nd.effn +
							":\n\t" +
							GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTClusterForGeneCluster(cluster));
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
			if (nd == null ) {
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
				double pval = post.getPvalue();
				if (pval < 0) {
					System.err.println(""
							+ "Cannot perform polytomy test with effective N (after polytomies) "+ effni +
							":\n\t" +
							node);
					node.setData("NA");
				} else {
					node.setData(df.format(pval));
				}
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
						//STBipartition[] biparts = nd.bipartitions;

						String lineTmp = node.getName() + "\t" + "t1" + "\t" + threequads[0].toString2() + "\t" + 
								Double.toString(postQ1) + "\t" + Double.toString(f1) +
								"\t" + Double.toString(effni);

						try {
							freqWriter.write(lineTmp + "\n");

							lineTmp = node.getName() + "\t" + "t2" + "\t" + threequads[1].toString2() + "\t" + 
									Double.toString(postQ2) + "\t" + Double.toString(f2) + 
									"\t" + Double.toString(effni);

							freqWriter.write(lineTmp + "\n");

							lineTmp = node.getName() + "\t" + "t3" + "\t" + threequads[2].toString2() + "\t" +
									Double.toString(postQ3) + "\t" + Double.toString(f3) +
									"\t" + Double.toString(effni);
							freqWriter.write(lineTmp + "\n");
						} catch (IOException e) {
							throw new RuntimeException(e);
						}

					}
				}
				//i++;
			} 
		}
		if (!nodeDataList.isEmpty())
			throw new RuntimeException("Hmm, this shouldn't happen; "+nodeDataList);
		if (this.getBranchAnnotation() == 16) {
			try {
				Rscript.write("#!/usr/bin/env Rscript\n");
				Rscript.write("red='#d53e4f';orange='#1d91c0';blue='#41b6c4';colormap = c(red,orange,blue)\n");
				Rscript.write("require(reshape2);require(ggplot2);\n");
				Rscript.write("dirPath = '.'; filePath = paste(dirPath"
						+ ",'/freqQuadCorrected.csv',sep=''); md<-read.csv(filePath,header=F,sep='\\t'); md$value = md$V5/md$V6;\n");
				Rscript.write("a<-length(levels(as.factor(md$V7)))*3.7; b<-4; sizes <- c(a,b);\n");
				Rscript.write("md$V8<-reorder(md$V8,-md$value)\n");
				Rscript.write("ggplot(data=md)+aes(x=V8,y=value,fill=V9)+"
						+ "geom_bar(stat='identity',color=1,width=0.8,position='dodge')+"
						+ "theme_bw()+theme(axis.text.x=element_text(angle=90))+scale_fill_manual"
						+ "(values=colormap,name='Topology')+geom_hline(yintercept=1/3,size=0.4,linetype=2)+"
						+ "ylab('relative freq.')+facet_wrap(~V7,scales='free_x')+xlab('')\n");
				Rscript.write("pdfFile = paste(dirPath,'/relativeFreq.pdf',sep=''); ggsave(pdfFile,width = sizes[1], height= sizes[2]);\n");
				Rscript.close();
				freqWriter.close();
			} catch (IOException e) {
				throw new RuntimeException("Hmm, the Rscript and frequency of Quadripartition files cannot be created!");
			}

		}
		System.err.println(st);
		return ret;
	}


	private STITreeCluster[] getSisterRemaining(STINode node) {
		STITreeCluster [] sisterRemaining = {null,null};
		Iterator<STINode> pcit = node.getParent().getChildren().iterator();
		STINode pc = pcit.next();
		if ( pc == node ) pc = pcit.next(); 
		sisterRemaining[0] = (STITreeCluster)pc.getData();
		if (node.getParent().isRoot() && node.getParent().getChildCount() == 3) {
			pc = pcit.next();
			if (pc == node) pc = pcit.next(); 
			sisterRemaining[1] = (STITreeCluster)pc.getData();;					
		}  else if (node.getParent().isRoot() && node.getParent().getChildCount() == 2) {
			if (pc.getChildCount() == 2) {
				Iterator<STINode> nieceIt = pc.getChildren().iterator();
				sisterRemaining[0] = (STITreeCluster) nieceIt.next().getData();
				sisterRemaining[1] = (STITreeCluster) nieceIt.next().getData();
			}
		} 
		else {
			sisterRemaining[1] = ((STITreeCluster)node.getParent().getData()).complementaryCluster();
		}
		return sisterRemaining;
	}


	@Override
	Long getTotalCost(Vertex all) {
		System.err.println("Normalized score (portion of input quartet trees satisfied before correcting for multiple individuals): " + 
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
