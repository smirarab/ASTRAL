package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeSet;

import phylonet.util.BitSet;

public abstract class AbstractMatrix implements Matrix{

    protected Integer n;
    protected float[][] matrix;
    List<TreeSet<Integer>> orderedTaxonBySimilarity;
    
    abstract int compareTwoValues(float a, float b);
    
    abstract Matrix factory(float [][] a);

    public AbstractMatrix() {
        super();
    }

    @Override
    public int getSize() {
        return n;
    }

    @Override
    public float get(int i, int j) {
        return this.matrix[i][j];
    }
    
    void assureOrderedTaxa () {
        if (this.orderedTaxonBySimilarity == null) {
            this.orderedTaxonBySimilarity = this.sortByDistance(this.matrix);
        }
    }
    

    /**
     * Returns the id of the closest taxon that is either present in presentBS or
     * has a smaller id than mssingId
     * @param presentBS
     * @param missingId
     * @return
     */
    @Override
    public int getClosestPresentTaxonId(BitSet presentBS, int missingId) {
        this.assureOrderedTaxa();
        int closestId = -1;
        for (Integer other: this.orderedTaxonBySimilarity.get(missingId)){
            if ( missingId > other // other is already added
                    || presentBS.get(other) // other was in original tree
                    ) {
                closestId = other;
                break;
            }
        }
        if (closestId == -1) {
            throw new RuntimeException("Bug: this should not be reached");
        }
        return closestId;
    }
    
    List<TreeSet<Integer>> sortByDistance(float[][] refMatrix) {
        List<TreeSet<Integer>> ret = new ArrayList<TreeSet<Integer>>(n);
        List<Integer> range = Utils.getRange(n);
        for (int i = 0; i < n; i++) {
            final float[] js = refMatrix[i];
            TreeSet<Integer> indices = sortColumn(range, js);
            ret.add(indices);
        }
        return ret;
    }
    
    TreeSet<Integer> sortColumn(List<Integer> range, final float[] js) {
        TreeSet<Integer> indices = new TreeSet<Integer>(new Comparator<Integer>() {

            @Override
            public int compare(Integer o1, Integer o2) {
                if (o1 == o2) {
                    return 0;
                }
                int comp = compareTwoValues(js[o1], js[o2]) ;
                return  comp == 0 ? - o1.compareTo(o2) : comp;
            }
        });
        indices.addAll(range);
        return indices;
    }

  //TODO: generate iterable, not list
    @Override
    public Iterable<BitSet> getQuadraticBitsets() {
        List<BitSet> newBitSets = new ArrayList<BitSet>();
        ArrayList<Integer> inds = new ArrayList<Integer> (n);
        for (int i = 0; i < n; i++) {
            inds.add(i);
        }
        for (final float[] fs : this.matrix) {
            Collections.sort(inds, new Comparator<Integer>() {
    
                @Override
                public int compare(Integer i1, Integer i2) {
                    if (i1 == i2) {
                        return 0;
                    }
                    int vc = compareTwoValues(fs[i1],fs[i2]);
                    if (vc != 0) {
                        return vc;
                    }
                    return i1 > i2 ? 1 : -1;
                }
            });
            BitSet stBS = new BitSet(n);
            //float previous = fs[inds.get(1)];
            //float lastStep = 0;
            for (int sp : inds) {
                stBS.set(sp);
                /*if (previous - fs[sp] < 0) {
                        continue;
                    }*/
                newBitSets.add((BitSet) stBS.clone());
                //lastStep = previous - fs[sp];
                //previous = fs[sp];
            }
            //System.err.println(this.clusters.getClusterCount());
        }
        return newBitSets;
    }
    
    
    
    /**
     * 
     * @param randomSample a map from species names to positions in the new induced matrix. 
     * @param id
     * @return
     */
    @Override
    public Matrix getInducedMatrix(HashMap<String, Integer> randomSample, TaxonIdentifier id) {

        int sampleSize = randomSample.size();
        float[][] matrix = new float[sampleSize][sampleSize];

        for (Entry<String, Integer> row : randomSample.entrySet()) {
            int rowI = id.taxonId(row.getKey());
            int i = row.getValue();
            for (Entry<String, Integer> col : randomSample.entrySet()) {
                int colJ = id.taxonId(col.getKey());
                matrix[i][col.getValue()] = this.matrix[rowI][colJ];
            }
        }

        return factory(matrix);
    }

    /*
     * 
     * TODO check what is this
     */
    public Matrix getInducedMatrix(List<Integer> sampleOrigIDs) {

        int sampleSize = sampleOrigIDs.size();
        float[][] similarityMatrixnew = new float [sampleSize][sampleSize];

        int i = 0;
        for (Integer rowI : sampleOrigIDs) {
            int j = 0;
            for (Integer colJ : sampleOrigIDs) {
                similarityMatrixnew[i][j] = this.matrix[rowI][colJ];
                j++;
            }
            i++;
        }
        Matrix ret = Factory.instance.newSimilarityMatrix(similarityMatrixnew);
        return ret;
    }
    
    public static float[][] convertType(AbstractMatrix in) {
        float[][] SM = new float[in.n][in.n];
        for (int i = 0; i < in.n; i++) {
            for (int j = 0; j < in.n; j++) {
                if (i == j) {
                    SM[i][j] = 0;
                } else if (in.matrix[i][j] == -99.0) {
                    SM[i][j] = -99;
                } else {
                    SM[i][j] = GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSpeciesCount() - in.matrix[i][j];
                }
            }
        }
        return SM;
    }

}