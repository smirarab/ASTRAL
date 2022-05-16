package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITree;

/**
 * Sets up the set X
 * 
 * @author smirarab
 * 
 */
public class WQDataCollectionConstrainedMP extends WQDataCollectionMP {
	
	@Override
	protected List<STITree> preProcessTreesBeforeAddingToX(STITree tre) {
		return TreeCompletion.treeCompletionRepeat((STITree)tre, (STITree)constraintTree.get(0),1);
	}

	@Override
	protected Collection<Tree> preProcessTreesBeforeAddingToX(Collection<Tree> trees) {
		Collection<Tree> out = new ArrayList<Tree>();
		for (Tree tree: trees)
			out.addAll(preProcessTreesBeforeAddingToX((STITree) tree));
		return out;
	}


	@Override
	protected void addFromConsensusTreesToX(Collection<Tree> allGreedies) {
		for (Tree tr: allGreedies){
			super.addBipartitionsFromSignleIndTreesToX(tr, null, GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTTaxonIdentifier(), false, true);
		}
	}

	@Override
	protected void addBipartitionsFromSignleIndTreesToX(Tree tr, Collection<Tree> baseTrees, TaxonIdentifier id, boolean polytomy, boolean saveClusteres) {
		for (STITree tre: TreeCompletion.treeCompletionRepeat((STITree)tr, (STITree)constraintTree.get(0),1))
			super.addBipartitionsFromSignleIndTreesToX(tre, baseTrees, id, polytomy, saveClusteres);
	}
	

	@Override
	protected boolean shouldDoQuadratic(int th, TNode greedyNode, int j) {
		return false; 
	}
	
	@Override
	public void addExtraBipartitionByDistance(){
		return;
	}

	private List<Tree> constraintTree;

	public WQDataCollectionConstrainedMP(HashClusterCollection clusters, AbstractInference<Tripartition> inference) {
		super(clusters, inference);
		this.constraintTree = inference.getConstraintTree();
//		if (multiind)
//			for(TNode l:this.constraintTree.get(0).postTraverse()){
//				if(l.isLeaf()){
//					GlobalMaps.taxonNameMap.getSpeciesIdMapper().getTaxaForSpecies(l.getID())
//					for()
//				}
//			}
//		for (int s = 0; s< GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSpeciesCount(); s++){
//    		List<Integer> stTaxa = GlobalMaps.taxonNameMap.getSpeciesIdMapper().getTaxaForSpecies(s);
//    		int tid = stTaxa.get(GlobalMaps.random.nextInt(stTaxa.size()));
//			//sampleSpecificTaxonIdentifier.taxonId(sampleNames.get(sampleNames.size()-1));
//    	}
	}
	
	
}
