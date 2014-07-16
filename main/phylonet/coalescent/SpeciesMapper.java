package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import phylonet.graph.model.MutableNode;
import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TMutableNode;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

public class SpeciesMapper {

    private int [] taxonIdToSpeciesId;
    private List<Set<Integer>> speciesIdtoTaxonId;
    private TaxonIdentifier speciesNameIdMap;

    public SpeciesMapper(int taxonCount) {
        this.taxonIdToSpeciesId = new int[taxonCount];
        this.speciesNameIdMap = new TaxonIdentifier();
        this.speciesIdtoTaxonId = new ArrayList<Set<Integer>>();
    }

    public TaxonIdentifier getSTTaxonIdentifier() {
        return this.speciesNameIdMap;
    }

    public int getSpeciesIdForTaxon(int id) {
        return this.taxonIdToSpeciesId[id];
    }

    private void setSpeciesIdForTaxon(int taxonId, int speciesId) {
        this.taxonIdToSpeciesId[taxonId] = speciesId;
        for (int i = this.speciesIdtoTaxonId.size(); i <= speciesId; i++) {
            this.speciesIdtoTaxonId.add(new TreeSet<Integer>());
        }
        this.speciesIdtoTaxonId.get(speciesId).add(taxonId);
        //System.err.println("Mapped taxon "+taxonId +" to species "+speciesId);
    }

    protected void setSpeciesIdForTaxon(int taxonId, String speciesName) {
        this.setSpeciesIdForTaxon(taxonId, this.speciesNameIdMap.taxonId(speciesName));
    }

    protected void setSpeciesIdForTaxon(String taxonName, String speciesName) {
        this.setSpeciesIdForTaxon(GlobalMaps.taxonIdentifier.taxonId(taxonName),
                this.speciesNameIdMap.taxonId(speciesName));
    }

    public Set<Integer> getTaxonsForSpecies(Integer species){
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
            for (Set<Integer> set : this.speciesIdtoTaxonId) {
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
            for (int j : this.getTaxonsForSpecies(i)) {
                gtbs.set(j);
            }
        }
        return gtbs;
    }
    
    public STITreeCluster getSTClusterForGeneCluster(STITreeCluster geneCluster) {
        STITreeCluster stCluster = new STITreeCluster(this.speciesNameIdMap);
        stCluster.setCluster(this.getSTBisetForGeneBitset(geneCluster.getBitSet()));
        return stCluster;
    }

    public STITreeCluster getGeneClusterForSTCluster(BitSet stBitset) {
        STITreeCluster geneCluster = new STITreeCluster();
        geneCluster.setCluster(this.getGeneBisetForSTBitset(stBitset));
        return geneCluster;
    }
    public STITreeCluster getGeneClusterForSTCluster(STITreeCluster stCluster) {
        return this.getGeneClusterForSTCluster(stCluster.getBitSet());
    }
    
    public void addMissingIndividuals(BitSet geneBS) {
        BitSet stBS = this.getSTBisetForGeneBitset(geneBS);
        for (int i = stBS.nextSetBit(0); i >=0 ; i = stBS.nextSetBit(i+1)) {
            Set<Integer> taxonsForSpecies = this.getTaxonsForSpecies(i);
            for (Iterator<Integer> it = taxonsForSpecies.iterator(); it.hasNext();) {
                Integer gi = it.next();
                geneBS.set(gi);               
            }
        }
    }

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
                node.setName(geneNamesForSpeciesName.get(0));
            } else {
                for (String gene: geneNamesForSpeciesName) {
                    node.createChild(gene);
                }
                node.setName(""); 
            }
        }
    }

}