package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Stack;
import java.util.Map.Entry;

import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.util.Trees;
import phylonet.util.BitSet;

public class SingleIndividualSample {
	
	private List<Integer> sampleOrigIDs;
	private List<String> sampleNames;
	private TaxonIdentifier tempTaxonId;
	private SimilarityMatrix similarityMatrix;
	private int sampleSize;

	public SingleIndividualSample(SpeciesMapper spm, SimilarityMatrix matrix) {
		sampleOrigIDs = new ArrayList<Integer>();
		sampleNames = new ArrayList<String>();
		tempTaxonId = new TaxonIdentifier();
    	for (int s = 0; s< spm.getSpeciesCount(); s++){
    		List<Integer> stTaxa = spm.getTaxaForSpecies(s);
    		int tid = stTaxa.get(GlobalMaps.random.nextInt(stTaxa.size()));
    		sampleOrigIDs.add(tid);
			sampleNames.add(GlobalMaps.taxonIdentifier.getTaxonName(tid));
			tempTaxonId.taxonId(sampleNames.get(sampleNames.size()-1));
    	}
		setSampleSize(sampleOrigIDs.size());
		
		this.similarityMatrix = matrix.getInducedMatrix(this.sampleOrigIDs);
	}
	
	

	public List<Tree> contractTrees(Iterable<Tree> intrees){
		List<Tree> outtrees = new ArrayList<Tree>();			
		for (Tree tr : intrees) { 
			STITree ntr = new STITree(tr);
			ntr.constrainByLeaves(sampleNames);
			outtrees.add(ntr);
		}
		return outtrees;
	}

	public TaxonIdentifier getTaxonIdentifier() {
		return this.tempTaxonId;
	}



	public SimilarityMatrix getSimilarityMatrix() {
		return this.similarityMatrix;
	}
	



	public int getSampleSize() {
		return sampleSize;
	}



	public void setSampleSize(int sampleSize) {
		this.sampleSize = sampleSize;
	}
	
	BitSet toOriginalBitSet(BitSet bs) {
		BitSet ret = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
		for (int j = bs.nextSetBit(0); 
				j >= 0; j = bs.nextSetBit(j+1)) {
			ret.set(this.sampleOrigIDs.get(j));
		}
		return ret;
	}
	
}
