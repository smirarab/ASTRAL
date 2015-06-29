package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TMutableNode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

public class SpeciesMapper {

    private int [] taxonIdToSpeciesId;
    private ArrayList<List<Integer>> speciesIdtoTaxonId;
    private ArrayList<Integer> speciesIdtoLowestTaxonId;
    private TaxonIdentifier speciesNameIdMap;

    public SpeciesMapper(int taxonCount) {
        this.taxonIdToSpeciesId = new int[taxonCount];
        this.speciesNameIdMap = new TaxonIdentifier();
        this.speciesIdtoTaxonId = new ArrayList<List<Integer>>();
        this.speciesIdtoLowestTaxonId = new ArrayList<Integer>();
    }

    public TaxonIdentifier getSTTaxonIdentifier() {
        return this.speciesNameIdMap;
    }

    public int getSpeciesIdForTaxon(int id) {
        return this.taxonIdToSpeciesId[id];
    }

    public int getLowestIndexIndividual(int stID){
    	return this.speciesIdtoLowestTaxonId.get(stID);
    }
    
    private void setSpeciesIdForTaxon(int taxonId, int speciesId) {
        this.taxonIdToSpeciesId[taxonId] = speciesId;
        for (int i = this.speciesIdtoTaxonId.size(); i <= speciesId; i++) {
            this.speciesIdtoTaxonId.add(new ArrayList<Integer>());
            this.speciesIdtoLowestTaxonId.add(null);
        }
        this.speciesIdtoTaxonId.get(speciesId).add(taxonId);
        if (this.speciesIdtoLowestTaxonId.get(speciesId) == null ||
        		this.speciesIdtoLowestTaxonId.get(speciesId) > taxonId) {
        	this.speciesIdtoLowestTaxonId.set(speciesId,taxonId);
        }
        //System.err.println("Mapped taxon "+taxonId +" to species "+speciesId);
    }

    protected void setSpeciesIdForTaxon(int taxonId, String speciesName) {
        this.setSpeciesIdForTaxon(taxonId, this.speciesNameIdMap.taxonId(speciesName));
    }

    protected void setSpeciesIdForTaxon(String taxonName, String speciesName) {
        this.setSpeciesIdForTaxon(GlobalMaps.taxonIdentifier.taxonId(taxonName),
                this.speciesNameIdMap.taxonId(speciesName));
    }

    protected List<Integer> getTaxaForSpecies(Integer species){
        return this.speciesIdtoTaxonId.get(species);
    }

    public int getSpeciesCount() {
        return this.speciesNameIdMap.taxonCount();
    }
    
    public String getSpeciesNames() {
        if (this.taxonIdToSpeciesId.length == this.speciesNameIdMap.taxonCount()) {
            return Arrays.toString(this.speciesNameIdMap.getAllTaxonNames());
        } else {
            HashMap<String, String> stToGtNameMap = new HashMap<String, String>();
            int i = 0;
            for (List<Integer> set : this.speciesIdtoTaxonId) {
                ArrayList<String> gtNames = new ArrayList<String>();
                for (Integer gi : set) {
                    gtNames.add(GlobalMaps.taxonIdentifier.getTaxonName(gi));
                }
                stToGtNameMap.put(this.getSpeciesName(i++), gtNames.toString());
            }
            return stToGtNameMap.toString();
        }
    }

    public String getSpeciesName(int stId) {
        return this.speciesNameIdMap.getTaxonName(stId);
    }

    public Integer speciesId(String name) {
        return this.speciesNameIdMap.taxonId(name);
    }

    protected BitSet getSTBisetForGeneBitset(BitSet bs) {
        BitSet stbs = new BitSet(this.getSpeciesCount());
        for (int i = bs.nextSetBit(0); i >=0 ; i = bs.nextSetBit(i+1)) {
            stbs.set(this.getSpeciesIdForTaxon(i));
        }
        return stbs;
    }
    
    protected BitSet getGeneBisetForSTBitset(BitSet bs) {
        BitSet gtbs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
        for (int i = bs.nextSetBit(0); i >=0 ; i = bs.nextSetBit(i+1)) {
            for (int j : this.getTaxaForSpecies(i)) {
                gtbs.set(j);
            }
        }
        return gtbs;
    }
    
    public STITreeCluster getSTClusterForGeneCluster(STITreeCluster geneCluster) {
        return this.getSTClusterForGeneBitSet(geneCluster.getBitSet());
    }
    
    public STITreeCluster getSTClusterForGeneBitSet(BitSet geneBitSet) {
        STITreeCluster stCluster = new STITreeCluster(this.speciesNameIdMap);
        stCluster.setCluster(this.getSTBisetForGeneBitset(geneBitSet));
        return stCluster;
    }

    public STITreeCluster getGeneClusterForSTCluster(BitSet stBitset) {
        STITreeCluster geneCluster = GlobalMaps.taxonIdentifier.newCluster();
        geneCluster.setCluster(this.getGeneBisetForSTBitset(stBitset));
        return geneCluster;
    }
    public STITreeCluster getGeneClusterForSTCluster(STITreeCluster stCluster) {
        return this.getGeneClusterForSTCluster(stCluster.getBitSet());
    }
    
    public void addMissingIndividuals(BitSet geneBS) {
        BitSet stBS = this.getSTBisetForGeneBitset(geneBS);
        for (int i = stBS.nextSetBit(0); i >=0 ; i = stBS.nextSetBit(i+1)) {
            List<Integer> taxonsForSpecies = this.getTaxaForSpecies(i);
            for (Iterator<Integer> it = taxonsForSpecies.iterator(); it.hasNext();) {
                Integer gi = it.next();
                geneBS.set(gi);               
            }
        }
    }

    /*
     * return whether bs has only taxa from one species.
     */
    public boolean isSingleSP(BitSet bs) {
        int i = bs.nextSetBit(0);
        int prevId =  this.getSpeciesIdForTaxon(i);
        for (; i >=0 ; i = bs.nextSetBit(i+1)) {
            int stID = this.getSpeciesIdForTaxon(i);
            if (prevId != stID) {
                return false;
            }
            prevId = stID;
        }
        return true;
    }
    
    public double meanSampling() {
    	return (this.taxonIdToSpeciesId.length+0.0)/this.getSpeciesCount();
    }
    
    public List<String> getGeneNamesForSpeciesName(String species) {
        ArrayList<String> ret = new ArrayList<String>();
        for (Integer id: this.speciesIdtoTaxonId.get(this.speciesId(species))) {
            ret.add(GlobalMaps.taxonIdentifier.getTaxonName(id));
        }
        return ret;
    }
    
    public void stToGt(MutableTree gt) {
        for (String leave : gt.getLeaves()) {
            TMutableNode node = gt.getNode(leave);
            List<String> geneNamesForSpeciesName = this.getGeneNamesForSpeciesName(node.getName());
            if (geneNamesForSpeciesName.size() == 1) {
            	if (! node.getName().equals(geneNamesForSpeciesName.get(0)))
            		node.setName(geneNamesForSpeciesName.get(0));
            } else {
                for (String gene: geneNamesForSpeciesName) {
                    node.createChild(gene);
                }
                node.setName(""); 
            }
        }
    }
    
    SimilarityMatrix convertToSpeciesDistance(SimilarityMatrix matrix) {
		float [][] STsimMatrix = new float[this.getSpeciesCount()][this.getSpeciesCount()];
		float[][] denum = new float[this.getSpeciesCount()][this.getSpeciesCount()];
		int n = matrix.getSize();
		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				int stI =  this.getSpeciesIdForTaxon(i);
				int stJ =  this.getSpeciesIdForTaxon(j);
				STsimMatrix[stI][stJ] += matrix.get(i,j); 
				STsimMatrix[stJ][stI] = STsimMatrix[stI][stJ];
				denum[stI][stJ] ++;
				denum[stJ][stI] ++;
			}
		}
		for (int i = 0; i < this.getSpeciesCount(); i++) {
			for (int j = 0; j < this.getSpeciesCount(); j++) {
				STsimMatrix[i][j] = denum[i][j] == 0 ? 0 : 
					STsimMatrix[i][j] / denum[i][j];
			}
			STsimMatrix[i][i] = 1;
			//System.err.println(Arrays.toString(this.distSTMatrix[i]));
		}
		System.err.println("Species tree distances calculated ...");
		
		SimilarityMatrix ret = new SimilarityMatrix(STsimMatrix);
		
		return ret;
	}
}