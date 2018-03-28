package phylonet.coalescent;

import java.util.HashMap;
import java.util.List;

import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

public interface Matrix {

    int getSize();

    float get(int i, int j);
    
    boolean isDistance();
    
    int getBetterSideByFourPoint(int x, int a, int b, int c) ;

    Iterable<BitSet> getQuadraticBitsets();

    int getClosestPresentTaxonId(BitSet presentBS, int missingId);

    Matrix getInducedMatrix(HashMap<String, Integer> randomSample, TaxonIdentifier id);
    
    List<BitSet> inferTreeBitsets();
    
    Matrix populate (List<STITreeCluster> treeAllClusters, List<Tree> geneTrees, SpeciesMapper spm);
    
    List<BitSet> resolvePolytomy(List<BitSet> bsList, boolean original);

}