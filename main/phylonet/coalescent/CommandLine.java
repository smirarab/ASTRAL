package phylonet.coalescent;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import phylonet.coalescent.GlobalMaps.TaxonNameMap;
import phylonet.tree.io.NewickReader;
import phylonet.tree.io.ParseException;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITree;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import com.martiansoftware.jsap.SimpleJSAP;
import com.martiansoftware.jsap.Switch;
import com.martiansoftware.jsap.stringparsers.FileStringParser;

public class CommandLine {
	
    protected static String _versinon = "4.4.4";


    private static void exitWithErr(String extraMessage, SimpleJSAP jsap) {
        System.err.println();
        System.err.println(extraMessage);
        System.err.println();
        System.err.println("Usage: java -jar astral."+_versinon+".jar "+ jsap.getUsage());
        System.err.println();
        System.err.println(jsap.getHelp());
        System.exit( 1 );
    }


    private static SimpleJSAP getJSAP() throws JSAPException {
        return new SimpleJSAP(
                "ASTRAL (version" + _versinon + ")",
                "species tree inference from unrooted gene trees. "
                + "The ASTRAL algorithm maximizes the number of shared quartet trees with"
                + " the collection of all gene trees. The result of this optimization problem"
                + " is statistically consistent under the multi-species coalesent model."
                + " This software can also solve MGD and MGDL problems (see options) instead of ASTRAL.",
                    
                new Parameter[] {
                    
                    new FlaggedOption("input file", 
                            FileStringParser.getParser().setMustExist(true), null, JSAP.REQUIRED, 'i', "input",
                            "a file containing input gene trees in newick format. (required)"),
                            
                    new FlaggedOption( "output file",
                            FileStringParser.getParser(), null, JSAP.NOT_REQUIRED, 'o', "output",
                            "a filename for storing the output species tree. Defaults to outputting to stdout."),
                            
                    new Switch("exact",
                            'x', "exact",
                            "find the exact solution by looking at all clusters - recommended only for small (<18) numer of taxa."),
   
                    new FlaggedOption("extra trees", 
                            FileStringParser.getParser().setMustExist(true), null, JSAP.NOT_REQUIRED, 'e', "extra",
                            "provide extra trees used to enrich the set of clusters searched"),
                            
                    new FlaggedOption("bootstraps", 
                            FileStringParser.getParser().setMustExist(true), null, JSAP.NOT_REQUIRED, 'b', "bootstraps",
                            "perform multi-locus bootstrapping using input bootstrap replicate files (use --rep to change the number of replications). "
                            + "The file given with this option should have a list of the gene tree bootstrap files, one per line, and each line corresponding to one gene. "
                            + "By default performs site-only resampling, but gene/site resampling can also be used. "),

                    new FlaggedOption("replicates", 
                            JSAP.INTEGER_PARSER, "100", JSAP.NOT_REQUIRED, 'r', "reps",
                            "Set the number of bootstrap replicates done in multi-locus bootstrapping. "),

                    new FlaggedOption("seed", 
                            JSAP.LONG_PARSER, "692", JSAP.NOT_REQUIRED, 's', "seed",
                            "Set the seed number used in multi-locus bootstrapping. "),

                    new Switch("gene-sampling",
                            'g', "gene-resampling",
                            "perform gene tree resampling in addition to site resampling. Useful only with the -b option."),
       

                    new Switch( "duplication",
                            'd', "dup",
                            "Solves MGD problem. Minimizes the number duplications required to explain "
                            + "gene trees using DynaDup algorithm (Bayzid, 2011). Note that with this option, "
                            + "DynaDyp would be used *instead of* ASTRAL."),
                            
                    new FlaggedOption( "duploss weight",
                            JSAP.STRING_PARSER, null, JSAP.NOT_REQUIRED, 'l', "duploss",
                            "Solves MGDL problem. Minimizes the number duplication and losses required"
                            + " to explain gene trees using DynaDup algorithm. Note that with this option, "
                            + "DynaDyp would be used *instead of* ASTRAL. "
                            + "Use -l 0 for standard (homomorphic) definition, and -l 1 for our new bd definition. "
                            + "Any value in between weights the impact of missing taxa somewhere between these two extremes. "
                            + "-l auto will automaticaly pick this weight. "), 
                });
    }


    public static void main(String[] args) throws Exception{
		
	long startTime = System.currentTimeMillis();
	
        SimpleJSAP jsap;		
        JSAPResult config;
        int criterion = 2; // 2 for ASTRAL, 0 for dup, 1 for duploss
		boolean rooted = false;
		boolean extrarooted = false;		
		Map<String, String> taxonMap = null;
		String replace = null;
		String pattern = null;
		List<Tree> trees = new ArrayList<Tree>();
		List<List<Tree>> bstrees = new ArrayList<List<Tree>>();
		List<List<Tree>> inputSets = new ArrayList<List<Tree>>();
		List<Tree> extraTrees = new ArrayList<Tree>();
		// STITree scorest = null;
		// boolean explore = false;
		// double proportion = 0.0D;
		// boolean exhaust = false;
		// double bootstrap = 1.0D;
	    // boolean unresolved = false;   
        // double time = -1.0D;
		// double wd = 1.0D;
		double cs = 1.0D;
		double cd = 1.0D;
		double wh = 1.0D;
		boolean exact = false;
		int k = 0;
        BufferedWriter outbuffer;
		
        jsap = getJSAP();     
        config = jsap.parse(args);  
        if ( jsap.messagePrinted() ) {
            exitWithErr("",jsap);
        }
        
        File outfile = config.getFile("output file");  
        if (outfile == null) {
            outbuffer = new BufferedWriter(new OutputStreamWriter(System.out));
        } else {
            outbuffer = new BufferedWriter(new FileWriter(outfile));
        }
        
        exact = config.getBoolean("exact");
        
        if (config.getBoolean("duplication") && config.contains("duploss weight")) {
            exitWithErr("dup and duploss options cannot be used together. Choose only one. ",jsap);
        }
        if (config.getBoolean("duplication")) {
            criterion = 0;
            rooted = true;
            extrarooted = true;
            System.err.println("Using DynaDup application, minimizing MGD (not ASTRAL).");
        }
        if (config.contains("duploss weight")) {
            criterion = 1;
            rooted = true;
            extrarooted = true;
            String v = config.getString("duploss weight");
            if (v.equals("auto")) {
                wh = -1;
            } else {
                wh = Double.parseDouble(v);
                if (wh < 0.0D || wh > 1.0D) {
                    exitWithErr("duploss weight has to be between 0 and 1",jsap);
                };
            }
            System.err.println("Using DynaDup application, minimizing MGDL (not ASTRAL).");
        }
        
        System.err.println("\n================== ASTRAL ===================== \n" );
        System.err.println("This is ASTRAL version " + _versinon);

        System.err.println("Gene trees are treated as " + (rooted ? "rooted" : "unrooted"));

        Random random = new Random(config.getLong("seed"));
        
        try {

            readInputTrees(new BufferedReader(new FileReader(config.getFile("input file"))),
                    rooted, trees, true);			
            k = trees.size();
            System.err.println(k+" trees read from " + config.getFile("input file"));
            
            GlobalMaps.taxonIdentifier.lock();  

        } catch (IOException e) {
            System.err.println("Error when reading trees.");
            System.err.println(e.getMessage());
            e.printStackTrace();
            return;
        }      
        
        if (trees.size() == 0) {
            System.err.println("Empty list of trees. The function exits.");
            return;
        }
        
        String outgroup = GlobalMaps.taxonIdentifier.getTaxonName(0);
        System.err.println("All outputted trees will be *arbitrarily* rooted at "+outgroup);
            
        try {
            
            if (config.getFile("extra trees") != null) {
                readInputTrees(new BufferedReader(new FileReader(config.getFile("extra trees"))), 
                        extrarooted, extraTrees, false);
                System.err.println(extraTrees.size() + " extra trees read from "
                        + config.getFile("extra trees"));
            }
            
            if (config.getFile("bootstraps") != null) {
                String line;
                BufferedReader rebuff = new BufferedReader(new FileReader(config.getFile("bootstraps")));
                while ((line = rebuff.readLine()) != null) {
                    ArrayList<Tree> g = new ArrayList<Tree>();
                    readInputTrees(new BufferedReader(new FileReader(line)), rooted, g, true);
                    Collections.shuffle(g, random);
                    bstrees.add(g);
                }
                rebuff.close();
            }

        } catch (IOException e) {
            System.err.println("Error when reading extra trees.");
            System.err.println(e.getMessage());
            e.printStackTrace();
            return;
        } 

		if (taxonMap != null) {
			GlobalMaps.taxonNameMap = new TaxonNameMap(taxonMap);
		} else if (replace != null) {	
			GlobalMaps.taxonNameMap = new TaxonNameMap (pattern, replace);
		}
	
		if (config.getFile("bootstraps") != null) {
	        System.err.println("Bootstrapping with seed "+config.getLong("seed"));
		    for (int i = 0; i < config.getInt("replicates"); i++) {
		        List<Tree> input = new ArrayList<Tree>();
		        inputSets.add(input );   
		        try {
    		        if (config.getBoolean("gene-sampling")) {
    		            for (int j = 0; j < k; j++) {
                            input.add(bstrees.get(random.nextInt(k)).remove(0));                 
                        }
    		        } else {   		        
    		            for (List<Tree> gene : bstrees) {
    		                input.add(gene.get(i));
    		            }
    		        }
		        } catch (IndexOutOfBoundsException e) {
		            exitWithErr("Error: You seem to have asked for "+config.getInt("replicates")+
		                    " but only "+ i +" replicaes could be created.\n" + 
		                    " Note that for gene resampling, you need more input bootstrap" +
		                    " replicates than the number of species tee replicates.", jsap);
		        }
		    }
		}

        int j = 0;
        List<Tree> bootstraps = new ArrayList<Tree>();
		for ( List<Tree> input : inputSets) {  
	        System.err.println("\n======== Running bootstrap replicate " + j++);
	        bootstraps.add(runOnOneInput(criterion, 
	                rooted, extrarooted, extraTrees, cs, cd,
                    wh, exact, outbuffer, input, null, outgroup));
		}
	    
		if (bootstraps != null && bootstraps.size() != 0) {
		    STITree<Double> cons = (STITree<Double>) Utils.greedyConsensus(bootstraps);
            cons.rerootTreeAtNode(cons.getNode(GlobalMaps.taxonIdentifier.getTaxonName(0)));
            Utils.computeEdgeSupports(cons, bootstraps);
            writeTreeToFile(outbuffer, cons);
		}
		
        System.err.println("\n======== Running the main analysis");
        runOnOneInput(criterion, rooted, extrarooted, extraTrees, cs, cd,
                wh, exact, outbuffer, trees, bootstraps, outgroup);
           
		outbuffer.close();
		
	    System.err.println("ASTRAL finished in "  + 
	            (System.currentTimeMillis() - startTime) / 1000.0D + " secs");
	}


    private static Tree runOnOneInput(int criterion, boolean rooted,
            boolean extrarooted, List<Tree> extraTrees, double cs, double cd,
            double wh, boolean exact, BufferedWriter outbuffer, List<Tree> input, 
            Iterable<Tree> bootstraps, String outgroup) {
        long startTime;
        startTime = System.currentTimeMillis();
        
        Inference inference =
            initializeInference(criterion, rooted, extrarooted, input,
                    extraTrees, cs, cd, wh, exact);
        
        /* if (scorest != null) {inference.scoreGeneTree(scorest);System.exit(0); }*/
        List<Solution> solutions = inference.inferSpeciesTree();
   
        System.err.println("Optimal tree inferred in "
        		+ (System.currentTimeMillis() - startTime) / 1000.0D + " secs");
        
        Tree st = solutions.get(0)._st;
        st.rerootTreeAtNode(st.getNode(outgroup));
   
        if ((bootstraps != null) && (bootstraps.iterator().hasNext())) {
            for (Solution solution : solutions) {
                Utils.computeEdgeSupports((STITree<Double>) solution._st, bootstraps);
            }
        }
        writeTreeToFile(outbuffer, solutions.get(0)._st);
        
        return st;
    }

    private static Inference initializeInference(int criterion, boolean rooted,
            boolean extrarooted, List<Tree> trees, List<Tree> extraTrees,
            double cs, double cd, double wh, boolean exact) {
        Inference inference;		
		if (criterion == 1 || criterion == 0) {
			inference = new DLInference(rooted, extrarooted, 
					trees, extraTrees,exact ,criterion > 0);			
		} else if (criterion == 2) {
			inference = new WQInference(rooted, extrarooted, 
					trees, extraTrees, exact,criterion > 0, -1);
		} else {
			throw new RuntimeException("criterion not set?");
		}		
		inference.setDLbdWeigth(wh); 
		inference.setCS(cs);
		inference.setCD(cd);
        return inference;
    }


    private static void readInputTrees(BufferedReader treeBufferReader, 
            boolean rooted, List<Tree> trees, boolean checkCompleteness)
            throws FileNotFoundException, IOException {
        String line;
        int l = 0;			
        try {
           while ((line = treeBufferReader.readLine()) != null) {
        		l++;
        		Set<String> previousTreeTaxa = new HashSet<String>();
        		if (line.length() > 0) {
        			line = line.replaceAll("\\)[^,);]*", ")");
        			NewickReader nr = new NewickReader(new StringReader(line));
        			if (rooted) {
        				STITree<Double> gt = new STITree<Double>(true);
        				nr.readTree(gt);
        				if (checkCompleteness) {
            				if (previousTreeTaxa.isEmpty()) {
            					previousTreeTaxa.addAll(Arrays.asList(gt
            							.getLeaves()));
            				} else {
            					if (!previousTreeTaxa.containsAll(Arrays.asList(gt
            							.getLeaves()))) {
            					    treeBufferReader.close();
            						throw new RuntimeException(
            								"Not all trees are on the same set of taxa: "
            										+ gt.getLeaves() + "\n"
            										+ previousTreeTaxa);
            					}
            				}
        				}
        				trees.add(gt);
        			} else {						
        				Tree tr = nr.readTree();
        				trees.add(tr);
        				String[] leaves = tr.getLeaves();
                        for (int i = 0; i < leaves.length; i++) {
                            //if (!stLablel) {
                                GlobalMaps.taxonIdentifier.taxonId(leaves[i]);
                            //} else {
                             //   GlobalMaps.taxonNameMap.getSpeciesIdMapper().speciesId(leaves[i]);
                            //}
                        }
        			}
        		}
        	}
        	treeBufferReader.close();
        } catch (ParseException e) {
        	treeBufferReader.close();
        	throw new RuntimeException("Failed to Parse Tree number: " + l ,e);
        }
    }


    private static void writeTreeToFile(BufferedWriter outbuffer, Tree t) {
        try {
		    outbuffer.write(t.toStringWD()+ " \n");
		    outbuffer.flush();
		} catch (IOException e) {
		    System.err.println("Error when writing the species tree");
		    System.err.println(e.getMessage());
		    e.printStackTrace();
		}
    }


    /*      System.out.println("\t-st species tree file: The file containing a species tree to be scored.\n" +
    "\t                       If this option is provided the software only scores the species tree.");
System.out.println("\t-a mapping file: The file containing the mapping from alleles to speceis if multiple alleles sampled.\n" +
    "\t                 Alternatively, two reqular expressions for automatic name conversion (optional)");
*/
//System.out.println("\t-u treat input gene trees as unrooted (Not implemented!)");
//System.out.println("\t-xu treat extra trees input gene trees as unrooted (Not implemented!)");
/*System.out.println("\t-cs and -cd: these two EXPERIMENTAL options set two parameters (cs and cd) to a value between 0 and 1. \n" +
   "\t    For any cluster C if |C| >= cs*|taxa|, we add complementary clusters (with respect to C) of all subclusters of C\n" +
   "\t    if size of the subcluster is >= cd*|C|.\n" +
   "\t    By default cs = cd = 1; so no extra clusters are added. Lower cs and cd values could result in better scores\n" +
   "\t    (especially when gene trees have missing data) but can also increase the running time quite substantially.");
*/
//System.out.println("\t-f perform fast and less-accurate subtree-bipartition based search (Not implemented!).");

/*      if (option[0].equals("-st")) {
    if (option.length != 2) {
        printUsage();
        return;
    }                   
    BufferedReader tmp = new BufferedReader(new FileReader(
            option[1]));
    line = tmp.readLine();
    NewickReader nr = new NewickReader(new StringReader(line));
    scorest = new STITree(true);
    nr.readTree(scorest);               
} else if (option[0].equals("-a")) {
    if ( (option.length != 2) && (option.length != 3)) {
        printUsage();
        return;
    }
    if (option.length == 2) {
        BufferedReader br = new BufferedReader(new FileReader(
                option[1]));

        taxonMap = new HashMap<String, String>();
        while ((line = br.readLine()) != null) {
            // String line;
            String[] mapString = line.trim().split(";");
            for (String s : mapString) {
                String species = s.substring(0, s.indexOf(":"))
                        .trim();
                s = s.substring(s.indexOf(":") + 1);
                String[] alleles = s.split(",");
                for (String allele : alleles) {
                    allele = allele.trim();
                    if (taxonMap.containsKey(allele)) {
                        System.err
                                .println("The input file is not in correct format");
                        System.err
                                .println("Any gene name can only map to one species");
                        System.exit(-1);
                    } else {
                        taxonMap.put(allele, species);
                    }
                }
            }
        }
        br.close();
    } else {
        pattern = option[1];
        rep = option [2];
        if (rep.equals("/delete/")) {
            rep = "";
        }
    }
} else if (option[0].equals("-cs")) {
    if (option.length != 2) {
        printUsage();
        return;
    }
    try {
        cs = Double.parseDouble(option[1]);
        if ((cs <= 1.0D) && (cs >= 0.0D))
            continue;
        printUsage();
        return;
    } catch (NumberFormatException e) {
        System.err.println("Error in reading parameter");
        printUsage();
        return;
    }

} else if (option[0].equals("-cd")) {
    if (option.length != 2) {
        printUsage();
        return;
    }
    try {
        cd = Double.parseDouble(option[1]);
        if ((cd <= 1.0D) && (cd >= 0.0D))
            continue;
        printUsage();
        return;
    } catch (NumberFormatException e) {
        System.err.println("Error in reading parameter");
        printUsage();
        return;
    }

} 
}

*/

}
