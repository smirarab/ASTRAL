package phylonet.coalescent;

import java.util.ArrayList;
import java.util.List;

import phylonet.tree.model.MutableTree;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITree;
import phylonet.util.BitSet;

/**
 * This class keeps track of a single-individual sample
 * of a multi-individual dataset. For each species, 
 * we will include only one of its individuals in any instance of this class. 
 * @author smirarab
 *
 */
public class SingleIndividualSample {
	
	/**
	 * IDs of sampled individuals in the original (global) taxon identifier
	 */
	private List<Integer> sampleGlobalIDs;
	/**
	 * Name of sampled individuals
	 */
	private List<String> sampleNames;
	/**
	 * A taxon identifier specifc to this subsample. This 
	 * taxon identifer will only include the individuals sampled. 
	 */
	//private TaxonIdentifier sampleSpecificTaxonIdentifier;
	// TODO: maybe we should take the distance matrix out of this class.
	//       not sure why it's here. 
	//private Matrix similarityMatrix;
	/**
	 * Size of the sample. 
	 */
	private int sampleSize;

	public SingleIndividualSample(SpeciesMapper spm, Matrix matrix) {
		sampleGlobalIDs = new ArrayList<Integer>();
		sampleNames = new ArrayList<String>();
		//sampleSpecificTaxonIdentifier = new TaxonIdentifier();
		//sampleSpecificTaxonIdentifier = GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTTaxonIdentifier();
		/*
		 * TODO: check if other parts of code need any changes
		 */
    	for (int s = 0; s< spm.getSpeciesCount(); s++){
    		List<Integer> stTaxa = spm.getTaxaForSpecies(s);
    		int tid = stTaxa.get(GlobalMaps.random.nextInt(stTaxa.size()));
    		sampleGlobalIDs.add(tid);
			sampleNames.add(GlobalMaps.taxonIdentifier.getTaxonName(tid));
			//sampleSpecificTaxonIdentifier.taxonId(sampleNames.get(sampleNames.size()-1));
    	}
		setSampleSize(sampleGlobalIDs.size());
		
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
	
	public Tree contractTree(Tree intree){	
		
		STITree ntr = new STITree(intree);
		ntr.constrainByLeaves(sampleNames); // sampleNames : GlobalMaps.taxonIdentifier.getTaxonName: gene tree names
		GlobalMaps.taxonNameMap.getSpeciesIdMapper().gtToSt2((MutableTree)ntr);
		return ntr;
	}
	

//	public TaxonIdentifier getTaxonIdentifier() {
//		return this.sampleSpecificTaxonIdentifier;
//	}




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
			ret.set(this.sampleGlobalIDs.get(j));
		}
		return ret;
	}
	
}
