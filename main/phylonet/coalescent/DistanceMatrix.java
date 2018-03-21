package phylonet.coalescent;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import phylonet.tree.io.NewickReader;
import phylonet.tree.io.ParseException;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

public class DistanceMatrix extends AbstractMatrix implements Matrix {

    public DistanceMatrix(int n) {
        this.n = n;
    }

    public DistanceMatrix(float[][] from) {
        this.n = from.length;
        this.matrix = from;
    }

    int compareTwoValues(float f1, float f2) {
        int vc = Float.compare(f1,f2);
        return vc;
    }
    
    public int getBetterSideByFourPoint(int x, int a, int b, int c) {
        double xa = this.matrix[x][a];
        double xb = this.matrix[x][b];
        double xc = this.matrix[x][c];
        double ab = this.matrix[a][b];
        double ac = this.matrix[a][c];
        double bc = this.matrix[b][c];
        double ascore = (xb + ac) - (xa + bc) ; 
        double bscore = (xa + bc) - (xb + ac); 
        double cscore = (xb + ac) - (xc + ab) ; 
        return ascore >= bscore ?
                ascore >= cscore ? a : c :
                    bscore >= cscore ? b : c;   
    }
    

    List<BitSet> PhyDstar(SpeciesMapper spm) {
        /* Write the distance matrix to a file, then use the file as input for PhyD*.
         * Call PhyD* to construct a tree. */
        File matrix;
        try {
            matrix = File.createTempFile("distance", ".matrix");
            BufferedWriter out = new BufferedWriter(new FileWriter(matrix));
            out.write(n + "\n");
            for (int i = 0; i < n; i++) {
                out.write(spm.getSpeciesName(i) + " ");
                for (int j = 0; j < n; j++) {
                    out.write(this.matrix[i][j] + " ");
                }
                out.write("\n");
            }
            out.close();
            String[] arg = new String[]{"java","-jar","PhyDstar.java","-i",matrix.getCanonicalPath()};
            System.err.println("Invoking PhyDstar with " + Arrays.toString(arg));
            PhyDstar.main(arg);
        } catch (IOException e) { throw new RuntimeException(); }
    
        Tree phyDtree = generatePhyDstarTree(matrix);
    
        List<BitSet> ret = new ArrayList<BitSet>();
        /* Extract list of BitSets of leaf-sets */
        for (TNode node : phyDtree.postTraverse()) {
            if (node.isRoot()) {
                continue;
            } else if (node.isLeaf()) {
                BitSet leaf = new BitSet();
                leaf.set(GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTTaxonIdentifier().taxonId(node.getName()));
                ((STINode)node).setData(leaf);
            } else {
                BitSet newbs = new BitSet(n);
                for (TNode cn : node.getChildren()) {
                    BitSet c = (BitSet) ((STINode)cn).getData();
                    newbs.or(c);
                }
                ((STINode)node).setData(newbs);
                ret.add(newbs);
            }
        }
        return ret;
    }

    List<BitSet> resolveByPhyDstar(List<BitSet> bsList, boolean original) {
        int size = bsList.size();
        float[][] internalMatrix = new float[size][size];
        HashMap<Integer, BitSet> bitMap = new HashMap<Integer, BitSet>(size);
    
        for (int n = 0; n < 10; n++) {
            int[][] pairLeaves = new int[size][2];
            List<int[]> ones = new ArrayList<int[]>();
    
            for (int i = 0; i < size; i++) {
                bitMap.put(i, bsList.get(i));
                BitSet bsI = bsList.get(i); 
                int[] bsIone = new int[bsI.cardinality()];
                int c = 0;
                for (int k = bsI.nextSetBit(0); k >= 0; k = bsI.nextSetBit(k + 1)) {
                    bsIone[c++] = k;
                }
                ones.add(bsIone);
                pairLeaves[i][0] = GlobalMaps.random.nextInt(bsI.cardinality());
                pairLeaves[i][1] = GlobalMaps.random.nextInt(bsI.cardinality());
            }
    
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    if (i == j) {
                        internalMatrix[i][j] = 0;
                        continue;
                    }
                    
                    float dI1I2 = this.matrix[ones.get(i)[pairLeaves[i][0]]][ones.get(i)[pairLeaves[i][1]]];
                    float dJ1J2 = this.matrix[ones.get(j)[pairLeaves[j][0]]][ones.get(j)[pairLeaves[j][1]]];
                    float dI1J1 = this.matrix[ones.get(i)[pairLeaves[i][0]]][ones.get(j)[pairLeaves[j][0]]];
                    float dI2J2 = this.matrix[ones.get(i)[pairLeaves[i][1]]][ones.get(j)[pairLeaves[j][1]]];
    
                    internalMatrix[i][j] += (dI1J1 + dI2J2 - dI1I2 - dJ1J2) / 2;
                }
            }
        }
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                internalMatrix[i][j] /= 10;
            }
        }
    
        File matrix;
        try {
            matrix = File.createTempFile("internal", ".matrix");
            BufferedWriter out = new BufferedWriter(new FileWriter(matrix));
            out.write(size + "\n");
            for (int i = 0; i < size; i++) {
                out.write(i + " ");
                for (int j = 0; j < size; j++) {
                    out.write(internalMatrix[i][j] + " ");
                }
                out.write("\n");
            }
            out.close();
    
            String[] arg = new String[]{"java","-jar","PhyDstar.java","-i",matrix.getCanonicalPath()};
            PhyDstar.main(arg);
        } catch (IOException e) { throw new RuntimeException(e); }
    
        Tree phyDtree = generatePhyDstarTree(matrix);
    
        List<BitSet> ret = new ArrayList<BitSet>();
    
        /* Extract list of BitSets of leaf-sets */
        for (TNode node : phyDtree.postTraverse()) {
            if (node.isRoot() ) {
                continue;
            } else if (node.isLeaf()) {
                if (original) {
                    BitSet leaf = bitMap.get(Integer.parseInt(node.getName()));
                    ((STINode)node).setData(leaf);
                } else {
                    BitSet leaf = new BitSet(size);
                    leaf.set(Integer.parseInt(node.getName()));
                    ((STINode)node).setData(leaf);
                }
            } else {
                BitSet newbs = new BitSet(n);
                for (TNode cn : node.getChildren()) {
                    BitSet c = (BitSet) ((STINode)cn).getData();
                    newbs.or(c);
                }
                ((STINode)node).setData(newbs);
                ret.add(newbs);
            }
        }
        return ret;
    }
    


    DistanceMatrix matricesByBranchDistance(List<Tree> geneTrees, SpeciesMapper spm) {
        this.matrix = new float[n][n]; // master gene matrix
        float[][] speciesSimilarityMatrix = new float[spm.getSpeciesCount()][spm.getSpeciesCount()]; // master species matrix
        float[][] pairNumMatrix = new float[n][n];
        float[][] speciesPN = new float[spm.getSpeciesCount()][spm.getSpeciesCount()];

        for (Tree tree : geneTrees) {
            HashMap<Integer, Integer> distanceMap = new HashMap<Integer, Integer>();
            float[][] geneSimilarityMatrix = new float[n][n];

            for (TNode node : tree.postTraverse()) {
                if (node.isLeaf()) {
                    BitSet tmp = new BitSet(n);
                    tmp.set(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
                    ((STINode)node).setData(tmp);

                    distanceMap.put(GlobalMaps.taxonIdentifier.taxonId(node.getName()), 0);
                } else {
                    List<BitSet> children = new ArrayList<BitSet>();
                    BitSet newbs = new BitSet(n);
                    for (TNode cn: node.getChildren()) {
                        BitSet c = (BitSet) ((STINode)cn).getData();
                        children.add(c);
                        newbs.or(c);
                    }
                    ((STINode)node).setData(newbs);

                    for (int i = 0; i < children.size() - 1; i++) {
                        BitSet left = children.get(i);
                        for (int j = i + 1; j < children.size(); j++) {
                            BitSet right = children.get(j);
                            for (int k = left.nextSetBit(0); k >= 0; k = left.nextSetBit(k + 1)) {
                                for (int l = right.nextSetBit(0); l >= 0; l = right.nextSetBit(l + 1)) {
                                    int d = distanceMap.get(k) + distanceMap.get(l) + 2;
                                    this.matrix[k][l] += d;
                                    this.matrix[l][k] = this.matrix[k][l];
                                    geneSimilarityMatrix[k][l] =  d;
                                    geneSimilarityMatrix[l][k] = geneSimilarityMatrix[k][l];
                                    pairNumMatrix[k][l] += 1;
                                    pairNumMatrix[l][k] = pairNumMatrix[k][l];
                                }
                            }
                        }
                    }

                    for (int index = 0; index < newbs.length(); index++) {
                        if (newbs.get(index)) {
                            distanceMap.put(index, distanceMap.get(index) + 1);
                        }
                    }
                }
            }

            SimilarityMatrix tmp = convertToSpeciesDistance(spm, geneSimilarityMatrix);
            for (int i = 0; i < spm.getSpeciesCount(); i++) {
                for (int j = i; j < spm.getSpeciesCount(); j++) {
                    speciesSimilarityMatrix[i][j] += tmp.matrix[i][j];
                    speciesSimilarityMatrix[j][i] = speciesSimilarityMatrix[i][j];
                    if (tmp.matrix[i][j] != 0) {
                        speciesPN[i][j] += 1;
                        speciesPN[j][i] = speciesPN[i][j];
                    }
                }
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    this.matrix[i][j] = 0;
                } else {
                    if (pairNumMatrix[i][j] != 0) {
                        this.matrix[i][j] /= pairNumMatrix[i][j];
                        this.matrix[i][j] =  this.matrix[i][j];
                    } else {
                        this.matrix[i][j] = -99;
                    }
                } 
            }
        }

        for (int i = 0; i < spm.getSpeciesCount(); i++) {
            for (int j = 0; j < spm.getSpeciesCount(); j++) {
                if (i == j) {
                    speciesSimilarityMatrix[i][j] = 0;
                } else {
                    if (speciesPN[i][j] != 0) {
                        speciesSimilarityMatrix[i][j] /= speciesPN[i][j];
                    } else {
                        speciesSimilarityMatrix[i][j] = -99;
                    }
                } 
            }
        }

        DistanceMatrix ret = new DistanceMatrix(speciesSimilarityMatrix);

        return ret;
    }

    @Override
    public boolean isDistance() {
        return true;
    }


    
    @Override
    public List<BitSet> inferTreeBitsets() {
       return PhyDstar(GlobalMaps.taxonNameMap.getSpeciesIdMapper());
    }

    @Override
    public Matrix populate(List<STITreeCluster> treeAllClusters, List<Tree> geneTrees, SpeciesMapper spm) {
       return matricesByBranchDistance(geneTrees, spm);
    }

    @Override
    public List<BitSet> resolvePolytomy(List<BitSet> bsList, boolean original) {
        return resolveByPhyDstar(bsList, original);
    }

    @Override
    Matrix factory(float[][] from) {
       return new DistanceMatrix(from);
    }

    private Tree generatePhyDstarTree(File matrix) {
    
            /* Write PhyD* tree back into ASTRAL */
            String newick;
            try {
                File phyDtree = new File(matrix.getCanonicalPath() + "_bionj.t");
                BufferedReader in = new BufferedReader(new FileReader(phyDtree));
                newick = in.readLine();
                in.close();
                matrix.delete();
                phyDtree.delete(); 
            } catch (IOException e) { throw new RuntimeException("Cannot find file: " + e); }
    
            /* Read the newick tree as an actual tree */
            Tree phyDstar_t = null;
    
            try {
    //          newick = newick.replaceAll("\\)[^,);]*", ")");
                NewickReader nr = new NewickReader(new StringReader(newick));
                phyDstar_t = nr.readTree();
    
            } catch (ParseException e) {
                throw new RuntimeException("Failed to Parse Tree: " , e);
            } catch (IOException e) {
                throw new RuntimeException();
            }
            return phyDstar_t;
        }
    

    // 339023690, 11085
    private void pupulateByBranchDistanceBS(List<Tree> geneTrees) {
        this.matrix = new float[n][n];
        float[][] pairNumMatrix = new float[n][n];
    
        for (Tree tree : geneTrees) {
            boolean[][] indicatorMatrix = new boolean[n][n];
            HashMap<Integer, Integer> distanceMap = new HashMap<Integer, Integer>();
    
            for (TNode node : tree.postTraverse()) {
                if (node.isLeaf()) {
    
                    // Setup BitSet
                    BitSet tmp = new BitSet(n);
                    tmp.set(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
                    ((STINode)node).setData(tmp);
    
                    // Update Map
                    distanceMap.put(GlobalMaps.taxonIdentifier.taxonId(node.getName()), 0);
                } else {
    
                    // Setup Bitset
                    List<BitSet> children = new ArrayList<BitSet>();
                    BitSet newbs = new BitSet(n);
                    for (TNode cn: node.getChildren()) {
                        BitSet c = (BitSet) ((STINode)cn).getData();
                        children.add(c);
                        newbs.or(c);
                    }
                    ((STINode)node).setData(newbs);
    
                    // Update similarity matrix
                    for (int i = 0; i < children.size() - 1; i++) {
                        BitSet left = children.get(i);
                        for (int j = i + 1; j < children.size(); j++) {
                            BitSet right = children.get(j);
                            for (int k = 0; k < left.length(); k++) {
                                if (left.get(k)) {
                                    for (int l = 0; l < right.length(); l++) {
                                        if (right.get(l)) {
                                            matrix[k][l] += distanceMap.get(k) + distanceMap.get(l) + 2;
                                            matrix[l][k] = matrix[k][l];
                                            pairNumMatrix[k][l] += 1;
                                            pairNumMatrix[l][k] = pairNumMatrix[k][l];
                                        }
                                    }
                                }
                            }
                        }
                    }
    
    
                    for (int index = 0; index < newbs.length(); index++) {
                        if (newbs.get(index)) {
                            distanceMap.put(index, distanceMap.get(index) + 1);
                        }
                    }
                }
            }           
        }
    
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    matrix[i][j] = 0;
                } else {
                    if (pairNumMatrix[i][j] != 0) {
                        matrix[i][j] /= pairNumMatrix[i][j];
                        matrix[i][j] = n - matrix[i][j];
                    } else {
                        matrix[i][j] = -99; 
                    }
                }
            }
        }
    }

    private SimilarityMatrix convertToSpeciesDistance(SpeciesMapper spm, float[][] SM) {
        float [][] STsimMatrix = new float[spm.getSpeciesCount()][spm.getSpeciesCount()];
        float[][] denum = new float[spm.getSpeciesCount()][spm.getSpeciesCount()];
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                int stI =  spm.getSpeciesIdForTaxon(i);
                int stJ =  spm.getSpeciesIdForTaxon(j);
                STsimMatrix[stI][stJ] += SM[i][j]; 
                STsimMatrix[stJ][stI] = STsimMatrix[stI][stJ];
                denum[stI][stJ] ++;
                denum[stJ][stI] ++;
            }
        }
    
        for (int i = 0; i < spm.getSpeciesCount(); i++) {
            for (int j = 0; j < spm.getSpeciesCount(); j++) {
                STsimMatrix[i][j] = denum[i][j] == 0 ? 0 : 
                    STsimMatrix[i][j] / denum[i][j];
            }
            STsimMatrix[i][i] = 0;
        }
        SimilarityMatrix ret = new SimilarityMatrix(STsimMatrix);
    
        return ret;
    }


}
