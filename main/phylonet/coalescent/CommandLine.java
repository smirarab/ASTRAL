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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import com.martiansoftware.jsap.SimpleJSAP;
import com.martiansoftware.jsap.Switch;
import com.martiansoftware.jsap.stringparsers.FileStringParser;

import phylonet.tree.io.NewickReader;
import phylonet.tree.io.ParseException;
import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TMutableNode;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.util.Trees;

public class CommandLine{

	protected String ASTRAL;
	protected String _version = "5.17.0";

	protected SimpleJSAP jsap;

	protected CommandLine () {
		Factory.instance = new FactoryAstral3();
		this.ASTRAL = "ASTRAL-III";
	}

	protected void exitWithErr(String extraMessage) {
		Logging.log("");
		Logging.log(extraMessage);
		Logging.log("");
		Logging.log("Usage: java -jar astral."+_version+".jar "+ jsap.getUsage());
		Logging.log("");
		Logging.log(jsap.getHelp());
		System.exit( 1 );
	}


	protected SimpleJSAP getJSAP() throws JSAPException {
		return new SimpleJSAP(
				ASTRAL  + "(version" + _version + ")",
				"species tree inference from unrooted gene trees. "
						+ "The ASTRAL algorithm maximizes the number of shared quartet trees with"
						+ " the collection of all gene trees. The result of this optimization problem"
						+ " is statistically consistent under the multi-species coalescent model."
						+ " This software can also solve MGD and MGDL problems (see options) instead of ASTRAL.",
						getAstralParameters());
	}


	protected Parameter[] getAstralParameters() {
		return new Parameter[] {

				new FlaggedOption("input file", 
						FileStringParser.getParser().setMustExist(true), null, JSAP.REQUIRED, 
						'i', "input",
						"a file containing input gene trees in newick format. (required)"),

				new FlaggedOption("output file",
						FileStringParser.getParser(), null, JSAP.NOT_REQUIRED,
						'o', "output",
						"a filename for storing the output species tree. Defaults to outputting to stdout."),

				new Switch("internode-dist", 'A', "internode",
						"USe NJst-like internode distances instead of quartet distance for building the search space (X). Unpublished work. "),

				new FlaggedOption("score species trees", 
						FileStringParser.getParser().setMustExist(true), null, JSAP.NOT_REQUIRED, 
						'q', "score-tree",
						"score the provided species tree and exit"),

				new FlaggedOption("branch annotation level", 
						JSAP.INTEGER_PARSER, "3", JSAP.NOT_REQUIRED,
						't', "branch-annotate",
						"How much annotations should be added to each branch: 0, 1, or 2. \n"
								+ "0: no annotations. \n"
								+ "1: only the quartet support for the main resolution. \n"
								+ "2: full annotation (quartet support, quartet frequency, and posterior probability for all three alternatives, "
								+ "plus total number of quartets around the branch and effective number of genes).\n"
								+ "3 (default): only the posterior probability for the main resolution.\n"
								+ "4: three alternative posterior probabilities.\n"
								+ "8: three alternative quartet scores.\n"
								+ "16/32: hidden commands useful to create a file called freqQuad.csv.\n"
								+ "10: p-values of a polytomy null hypothesis test."),
				
				new Switch("subunit",
						'u', "subunit",
						"Output branch lengths computed by taking the median of input gene trees (often in substitution units) "),

				new FlaggedOption("bootstraps", 
						FileStringParser.getParser().setMustExist(true), null, JSAP.NOT_REQUIRED,
						'b', "bootstraps",
						"perform multi-locus bootstrapping using input bootstrap replicate files (use --rep to change the number of replications). "
								+ "The file given with this option should have a list of the gene tree bootstrap files, one per line, and each line corresponding to one gene. "
								+ "By default performs site-only resampling, but gene/site resampling can also be used. "),

				new FlaggedOption("replicates", 
						JSAP.INTEGER_PARSER, "100", JSAP.NOT_REQUIRED, 
						'r', "reps",
						"Set the number of bootstrap replicates done in multi-locus bootstrapping. "),

				new FlaggedOption("seed", 
						JSAP.LONG_PARSER, "692", JSAP.NOT_REQUIRED,
						's', "seed",
						"Set the seed number used in multi-locus bootstrapping. "),

				new Switch("gene-sampling",
						'g', "gene-resampling",
						"perform gene tree resampling in addition to site resampling. Useful only with the -b option."),

				new Switch("gene-only",
						JSAP.NO_SHORTFLAG, "gene-only",
						"perform bootstrapping but only with gene tree resampling. Should not be used with the -b option."),    

				new FlaggedOption("keep", 
						JSAP.STRING_PARSER, null, JSAP.NOT_REQUIRED, 
						'k', "keep",
						" -k completed: outputs completed gene trees (i.e. after adding missing taxa) to a file called [output file name].completed_gene_trees.\n"
								+ " -k completed_norun: outputs completed gene trees (i.e. after adding missing taxa) to a file called [output file name].completed_gene_trees.\n"
								+ " -k bootstraps: outputs individual bootstrap replicates to a file called [output file name].[i].bs\n"
								+ " -k bootstraps_norun: just like -k bootstraps, but exits after outputting bootstraps.\n"
								+ " -k searchspace_norun: outputs the search space and exits; use -k searchspace to continue the run after outputting the search space.\n"
								+ " -k compatible_trees: outputs gene trees made compatible with constraint tree to a file [output file name].compatible_gene_trees\n"
								+ " -k compatible_trees_norun: outputs gene trees made compatible with constraint tree to a file and exit [output file name].compatible_gene_trees\n"	                 
								+ "When -k option is used, -o option needs to be given. "
								+ "The file name specified using -o is used as the prefix for the name of the extra output files.").setAllowMultipleDeclarations(true),

				new FlaggedOption("outgroup", 
						JSAP.STRING_PARSER, null, JSAP.NOT_REQUIRED, 
						JSAP.NO_SHORTFLAG, "outgroup",
						" choose a single species to be used as outgroup FOR DISPLAY PUROPSES ONLY (has no effect on the actual unrooted tree inferred) "),

				new FlaggedOption("lambda", 
						JSAP.DOUBLE_PARSER, "0.5", JSAP.NOT_REQUIRED,
						'c', "lambda",
						"Set the lambda parameter for the Yule prior used in the calculations"
								+ " of branch lengths and posterior probabilities. Set to zero to get ML branch "
								+ "lengths instead of MAP."
								+ " Higher values tend to shorten estimated branch lengths and very"
								+ " high values can give inaccurate results (or even result in underflow)."),

				new FlaggedOption("mapping file", 
						FileStringParser.getParser().setMustExist(true), null, JSAP.NOT_REQUIRED, 
						'a', "namemapfile",
						"a file containing the mapping between names in gene tree and names in the species tree. "
								+ "The mapping file has one line per species, with one of two formats:\n"
								+ " species: gene1,gene2,gene3,gene4\n"
								+ " species 4 gene1 gene2 gene3 gene4\n"),

                new FlaggedOption("constraint tree", 
                        FileStringParser.getParser().setMustExist(true), null, JSAP.NOT_REQUIRED, 
                        'j', "constrain",
                        "provide constraint tree to be enforced on the search space and hence the final tree"),

				new FlaggedOption("minleaves", 
						JSAP.INTEGER_PARSER, null, JSAP.NOT_REQUIRED, 
						'm', "minleaves",
						"Remove genes with less than specified number of leaves "),

				new FlaggedOption("samplingrounds", 
						JSAP.INTEGER_PARSER, null, JSAP.NOT_REQUIRED, 
						JSAP.NO_SHORTFLAG, "samplingrounds",
						"For multi-individual datasets, perform these many rounds of individual sampling for"
								+ " building the set X. The program"
								+ " automatically picks this parameter if not provided or if below one."),

				new FlaggedOption("gene repetition", 
						JSAP.INTEGER_PARSER, "1", JSAP.NOT_REQUIRED,
						'w', "generepeat",
						"the number of trees sampled for each locus. "),

				new FlaggedOption("polylimit", 
						JSAP.INTEGER_PARSER, null, JSAP.NOT_REQUIRED, 
						JSAP.NO_SHORTFLAG, "polylimit",
						"Sets a limit for size of polytomies in greedy consensus trees where O(n) number"
								+ " of new  resolutions are added. ASTRAL-III sets automatic limits to guarantee polynomial"
								+ " time running time."),

				new Switch("exact",
						'x', "exact",
						"find the exact solution by looking at all clusters - recommended only for small (<18) number of taxa."),

				new Switch("rename",
						'R', "rename",
						"Simply rename gene trees according to species names given with -a; using the output can save memory as opposed to using the original file."),

				new FlaggedOption("extraLevel",
						JSAP.INTEGER_PARSER, "1", JSAP.NOT_REQUIRED,
						'p', "extraLevel",
						"How much extra bipartitions should be added: 0, 1, 2, or 3. "
								+ "0: adds nothing extra. "
								+ "1 (default): adds to X but not excessively (greedy resolutions). "
								+ "2: adds a potentially large number and therefore can be slow (quadratic distance-based)."
								+ "3: similar to default, but instead of completing input gene trees, it uses -f and -e as complted gene trees."),

				new FlaggedOption("extra trees", 
						FileStringParser.getParser().setMustExist(true), null, JSAP.NOT_REQUIRED, 
						'e', "extra",
						"provide extra trees (with gene labels) used to enrich the set of clusters searched"),

				new FlaggedOption("extra species trees", 
						FileStringParser.getParser().setMustExist(true), null, JSAP.NOT_REQUIRED, 
						'f', "extra-species",
						"provide extra trees (with species labels) used to enrich the set of clusters searched"),

				new FlaggedOption("remove extra tree bipartitions", 
						FileStringParser.getParser().setMustExist(true), null, JSAP.NOT_REQUIRED, 
						JSAP.NO_SHORTFLAG, "remove-bipartitions",
						"removes bipartitions of the provided extra trees (with species labels)"),

				new FlaggedOption("trimming threshold", 
						JSAP.DOUBLE_PARSER, "0", JSAP.NOT_REQUIRED,
						'd', "trimming",
						"trimming threshold is user's estimate on normalized score; the closer user's estimate is, the faster astral runs."),

		};
	}


	protected Options readOptions( boolean rooted, boolean extrarooted, double wh,
			JSAPResult config, List<Tree> mainTrees, List<List<String>> bootstrapInputSets) 
					throws JSAPException, IOException {

		Map<String, String> taxonMap = null;
		String replace = null;
		String pattern = null;
		Integer minleaves = null;
		Integer samplingrounds = null;
		Integer polylimit = null;
		String outfileName = null;
		Set<String> keepOptions = new HashSet<String>();
		String freqPath = null;        
		List<List<String>> bstrees = new ArrayList<List<String>>();
		int k = 0;


		File outfile = config.getFile("output file");  

		if (config.getBoolean("gene-only") && config.getFile("bootstraps") != null) {
			exitWithErr("--gene-only and -b cannot be used together");
		}

		if (outfile == null) {
			if (config.getInt("branch annotation level") % 16 == 0) {
				File extraTreeFile = config.getFile("score species trees");
				freqPath = extraTreeFile.getAbsoluteFile().getParentFile().getAbsolutePath();
			}
		} else {
			if (config.getInt("branch annotation level") % 16 == 0) {
				freqPath = outfile.getAbsoluteFile().getParentFile().getAbsolutePath();
			}
			outfileName = config.getFile("output file") == null ? null
					: config.getFile("output file").getCanonicalPath();
		}

		setupComputing(config);

		Logging.startLogger();

		if (config.getFile("mapping file") != null) {
			BufferedReader br = new BufferedReader(new FileReader(
					config.getFile("mapping file")));

			taxonMap = new HashMap<String, String>();
			String s;
			try {
				while ((s = br.readLine()) != null) {
					s = s.trim();
					if ("".equals(s)) {
						continue;
					}
					String species;
					String[] alleles;
					if ("".equals(s.trim()))  
						continue;
					if (s.indexOf(":") != -1) {
						species = s.substring(0, s.indexOf(":")).trim();
						s = s.substring(s.indexOf(":") + 1);
						alleles = s.split(",");
					} else {
						alleles = s.split(" ",3);
						species = alleles[0];
						alleles = alleles[2].split(" ");
					}
					for (String allele : alleles) {
						allele = allele.trim();
						if (taxonMap.containsKey(allele)) {
							Logging.log("The name mapping file is not in the correct format");
							Logging.log("A gene name can map to one only species name; check: " + allele + " which seems to appear at least twice: " +
									taxonMap.get(allele) + " & " + species);
							System.exit(-1);
						} else if (alleles.length > 1 && allele.equals(species)) {
							Logging.log("Error: The species name cannot be identical to gene names when "
									+ "multiple alleles exist for the same gene: "+ allele);
							System.exit(-1);
						}
						taxonMap.put(allele, species);
					}
				}
			} catch (Exception e) {
				br.close();
				throw new RuntimeException("\n** Error **: Your name mapping file looks incorrect.\n   Carefully check its format. ", e);
			}
			br.close();
		}

		minleaves = config.contains("minleaves")? config.getInt("minleaves"):null;      
		samplingrounds = config.contains("samplingrounds")? config.getInt("samplingrounds"):null;        
		polylimit = config.contains("polylimit")? config.getInt("polylimit"):null;

		try {

			//GlobalMaps.taxonIdentifier.taxonId("0");

			//Logging.log("Main input file: "+config.getFile("input file"));
			readInputTrees(mainTrees,
					readTreeFileAsString(config.getFile("input file")),
					rooted, true, false, minleaves, 
					config.getInt("branch annotation level"), null);			
			Logging.log(mainTrees.size() + " trees read from " + config.getFile("input file"));

			GlobalMaps.taxonIdentifier.lock();  
			Logging.logTimeMessage("");
		} catch (IOException e) {
			Logging.log("Error when reading trees.");
			Logging.log(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		} 

		if (mainTrees == null || mainTrees.size() == 0) {
			Logging.log("Empty list of trees. The function exits.");
			System.exit(1);
		} else {
			k = mainTrees.size();
		}

		if (taxonMap != null) {
			GlobalMaps.taxonNameMap = new TaxonNameMap(taxonMap);
		} else if (replace != null) {   
			GlobalMaps.taxonNameMap = new TaxonNameMap (pattern, replace);
		} else {
			GlobalMaps.taxonNameMap = new TaxonNameMap();
		}

		if (config.getStringArray("keep") != null && config.getStringArray("keep").length != 0) {
			if (outfileName == null) {
				throw new JSAPException("When -k option is used, -o is also needed.");
			}
			for (String koption : config.getStringArray("keep")) {
				if ("completed".equals(koption) ||
						"bootstraps".equals(koption) ||
						"bootstraps_norun".equals(koption)||
						"searchspace_norun".equals(koption)||
						"completed_norun".equals(koption) ||
						"searchspace".equals(koption) ||
						"compatible_trees".equals(koption)||
						"compatible_trees_norun".equals(koption)) {
					keepOptions.add(koption);
				} else {
					throw new JSAPException("-k "+koption+" not recognized.");
				}
			}
		}

		try {           
			if (config.getFile("bootstraps") != null) {
				String line;
				BufferedReader rebuff = new BufferedReader(new FileReader(config.getFile("bootstraps")));
				while ((line = rebuff.readLine()) != null) {
					List<String> g = readTreeFileAsString(new File(line));
					Collections.shuffle(g, GlobalMaps.random);
					bstrees.add(g);
				}
				rebuff.close();
			}

		} catch (IOException e) {
			Logging.log("Error when reading bootstrap trees.");
			Logging.log(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		} 


		if (config.getFile("bootstraps") != null || config.getBoolean("gene-only")) {
			Logging.log("Bootstrapping with seed " + config.getLong("seed"));
			for (int i = 0; i < config.getInt("replicates"); i++) {
				List<String> input = new ArrayList<String>();
				bootstrapInputSets.add(input);   
				try {
					if (config.getBoolean("gene-sampling")) {
						for (int j = 0; j < k; j++) {
							input.add(bstrees.get(GlobalMaps.random.nextInt(k)).remove(0));                 
						}
					} else if (config.getBoolean("gene-only")) { 
						for (int j = 0; j < k; j++) {
							input.add(mainTrees.get(GlobalMaps.random.nextInt(k)).toString());                 
						}	
					} else {
						for (List<String> gene : bstrees) {
							input.add(gene.get(i));
						}
					}
				} catch (IndexOutOfBoundsException e) {
					exitWithErr("Error: You seem to have asked for " + config.getInt("replicates") +
							" but only " + i + " replicates could be created.\n" + 
							" Note that for gene resampling, you need more input bootstrap" +
							" replicates than the number of species tree replicates.");
				}
				if (keepOptions.contains("bootstraps_norun") ||
						keepOptions.contains("bootstraps")) {
					String bsfn = outfile + ( "." + i + ".bs" );
					BufferedWriter bsoutbuffer = new BufferedWriter(new FileWriter(bsfn));
					for (String tree: input) {
						bsoutbuffer.write(tree + " \n");
					}
					bsoutbuffer.close();
				}
			}
			if (keepOptions.contains("bootstraps_norun") ||
					keepOptions.contains("bootstraps")) {
				Logging.log("bootstrap files written to files "+ outfile + ("." + 0 + ".bs") + 
						" to "+outfile + ( "." + config.getInt("replicates") + ".bs" ));
			}
			if (keepOptions.contains("bootstraps_norun")) {
				Logging.log("Exiting after outputting the bootstrap files");
				System.exit(0);
			}
		}


		Options options = new Options(rooted, extrarooted, 
				config.getBoolean("exact"), 1, 
				config.getInt("extraLevel"),
				keepOptions.contains("completed") || keepOptions.contains("completed_norun"), 
				keepOptions.contains("searchspace_norun") || keepOptions.contains("searchspace"), 
				!keepOptions.contains("searchspace_norun") && !keepOptions.contains("completed_norun"),
				config.getInt("branch annotation level"), 
				config.getDouble("lambda"),
				outfileName, samplingrounds == null ? -1 : samplingrounds, polylimit == null ? -1 : polylimit,
						config.getDouble("trimming threshold"), freqPath, minleaves,
						config.getInt("gene repetition"), 
						config.contains("remove extra tree bipartitions"),
						config.getBoolean("internode-dist"),
						config.getBoolean("subunit"), 
						keepOptions.contains("compatible_trees"), keepOptions.contains("compatible_trees_norun"));
		options.setDLbdWeigth(wh); 
		options.setCS(1d);
		options.setCD(1d);

		return options;
	}

	protected void setupComputing(JSAPResult config) {
	}


	public void process(String[] args) throws Exception{
		try {
			long startTime = System.currentTimeMillis();


			JSAPResult config;
			boolean rooted = false;
			boolean extrarooted = false;
			double wh = 1.0D;

			List<Tree> mainTrees = new ArrayList<Tree>();
			List<List<String>> bootstrapInputSets = new ArrayList<List<String>>();
			BufferedWriter outbuffer;

			Logging.initalize(Factory.instance.newLogger());
			Logging.log("\n================== "+ ASTRAL +" ===================== \n" );

			Logging.log("This is "+ ASTRAL +" version " + _version);

			checkLibraries();

			jsap = getJSAP();     
			config = jsap.parse(args);  
			if ( jsap.messagePrinted() ) {
				exitWithErr("");
			}

			Logging.log("Gene trees are treated as " + (rooted ? "rooted" : "unrooted"));


			GlobalMaps.random = new Random(config.getLong("seed"));
			GlobalMaps.taxonIdentifier = Factory.instance.newTaxonIdentifier();

			Options options = readOptions( rooted, extrarooted, wh, config,
					mainTrees, bootstrapInputSets);

			File outfile = config.getFile("output file");  
			if (outfile == null) {
				outbuffer = new BufferedWriter(new OutputStreamWriter(System.out));

			} else {

				outbuffer = new BufferedWriter(new FileWriter(outfile));
			}

			GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTTaxonIdentifier().lock();
			String outgroup = config.getString("outgroup") == null ? GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSpeciesName(0): config.getString("outgroup");
			//Logging.log("index"+ GlobalMaps.taxonNameMap.getSpeciesIdMapper().speciesId(outgroup));

			List<String> toScore = null;

			if (config.getBoolean("rename")) {
				renmaeFromGTtoST(mainTrees, outbuffer);
			} else if (config.getFile("score species trees") != null) {
	        	if(config.getFile("constraint tree") != null)
	        		Logging.log("!!Constraint tree will be ignored!!");
				Logging.log("Scoring " + config.getFile("score species trees"));
				toScore = readTreeFileAsString(config.getFile("score species trees"));
				runScore(rooted, mainTrees, outbuffer,
						options, outgroup, toScore);
			} else {

				runInference(config, rooted, extrarooted, 
						mainTrees, outbuffer, bootstrapInputSets,  options, outgroup);
			}
			// TODO: debug info
			Logging.log("Weight calculation took " + PolytreeA3.getTime() / 1000000000.0D + " secs");

			Logging.log("ASTRAL finished in " + 
					(System.currentTimeMillis() - startTime) / 1000.0D + " secs");
		} catch (Exception e) {
			throw (e);
		}
	}


	protected void checkLibraries() {
	}

	protected void runScore(boolean rooted,
			List<Tree> mainTrees,
			BufferedWriter outbuffer, Options options, String outgroup,
			List<String> toScore) throws FileNotFoundException, IOException {
		Logging.log("Scoring: " + toScore.size() +" trees");

		AbstractInference inference = Factory.instance.newInference(options, mainTrees, new ArrayList<Tree>(), new ArrayList<Tree>());
		double score = Double.NEGATIVE_INFINITY;
		List<Tree> bestTree = new ArrayList<Tree>(); 
		for (String trs : toScore) {   
			List<Tree> trees = new ArrayList<Tree>();
			readInputTrees(trees, Arrays.asList(new String[]{trs}),
					rooted, true, true, null, 1, false? 
							outgroup: null);
			Tree tr = trees.get(0);

			double nscore = inference.scoreSpeciesTreeWithGTLabels(tr, true);

			if (nscore > score) {
				score = nscore;
				bestTree.clear();
				bestTree.add(tr);
			} else if (nscore == score) {
				bestTree.add(tr);
			}

			if (!GlobalMaps.taxonNameMap.getSpeciesIdMapper().isSingleIndividual()) {
				Logging.log("Scored tree with gene names:\n"+tr.toNewickWD());
			}

			GlobalMaps.taxonNameMap.getSpeciesIdMapper().gtToSt((MutableTree) tr);

			if (options.getBranchannotation() != 12) {
				writeTreeToFile(outbuffer, tr);
			} 
		}
		if (options.getBranchannotation() == 12) {
			for (Tree bt: bestTree)
				writeTreeToFile(outbuffer, bt);
		}

		outbuffer.close();
	}


	protected void runInference(JSAPResult config,
			boolean rooted, boolean extrarooted,
			List<Tree> mainTrees, BufferedWriter outbuffer,
			List<List<String>> bootstrapInputSets, 
			Options options, String outgroup) throws JSAPException, IOException,
	FileNotFoundException {

		Logging.log("All output trees will be *arbitrarily* rooted at "+outgroup);

		List<Tree> extraTrees = new ArrayList<Tree>();
		List<Tree> constraintTree = new ArrayList<Tree>();
		//List<Tree> completedTrees = new ArrayList<Tree>();
		List<Tree> toRemoveExtraTrees = new ArrayList<Tree>();

		if (config.getFile("constraint tree") != null) {
			if (config.getFile("mapping file") != null) {
				Logging.log("Constrained version can not be used with Multi-individual version of ASTRAL for now. \n"+
						     "Consider providing your species definitions as a constraint to normal ASTRAL." );
				System.exit(0);
			}
			// TODO: -q: warn that -j is ignored but proceed. 
			// TODO: error on -x
			
		}
		try {

		    if (config.getFile("constraint tree") != null) {
		    
		    	readInputTrees(constraintTree, 
		        	readTreeFileAsString(config.getFile("constraint tree")), 
		                false, false, false, null, 1, null);///rooting
		        Logging.log(constraintTree.size() + " constraint tree read from "
		                + config.getFile("constraint tree"));
		        
	        	Logging.log("All gene trees are converted to be compatible with species tree.");
	        	
		    }
		    
			if (config.getFile("extra trees") != null) {
				readInputTrees(extraTrees, 
						readTreeFileAsString(config.getFile("extra trees")), 
						extrarooted, true, false, null, 1, null);
				Logging.log(extraTrees.size() + " extra trees read from "
						+ config.getFile("extra trees"));
			}

			if (config.getFile("extra species trees") != null) {
				readInputTrees(extraTrees,
						readTreeFileAsString(config.getFile("extra species trees")), 
						extrarooted, true, true, null, 1, null);
				Logging.log(extraTrees.size() + " extra trees read from "
						+ config.getFile("extra trees"));
			}

			if (config.getFile("remove extra tree bipartitions") != null) {
				readInputTrees(toRemoveExtraTrees,
						readTreeFileAsString(config.getFile("remove extra tree bipartitions")), 
						true, true, true, null, 1, null);
				Logging.log(toRemoveExtraTrees.size() + " extra trees to remove from search space read from "
						+ config.getFile("remove extra tree bipartitions"));
			}

		} catch (IOException e) {
			Logging.log("Error when reading extra trees.");
			Logging.log(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		} 


		int j = 0;
		List<Tree> bootstraps = new ArrayList<Tree>();
		for ( List<String> input : bootstrapInputSets) {  
			Logging.log("\n======== Running bootstrap replicate " + j++);
			List<Tree> trees = new ArrayList<Tree>();
			readInputTrees(trees, input, rooted, false, false, options.getMinLeaves(),
					config.getInt("branch annotation level"), null);
			bootstraps.add(runOnOneInput( extraTrees,toRemoveExtraTrees, outbuffer, trees, null, outgroup, options, constraintTree));
		}

		if (bootstraps != null && bootstraps.size() != 0) {
			STITree<Double> cons = (STITree<Double>) Factory.instance.greedyCons().greedyConsensus(bootstraps, false, 
					GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTTaxonIdentifier(), false);
			cons.rerootTreeAtNode(cons.getNode(outgroup));
			Trees.removeBinaryNodes(cons);
			Utils.computeEdgeSupports(cons, bootstraps);
			writeTreeToFile(outbuffer, cons);
		}
		Logging.logTimeMessage(" ");
		Logging.log("\n======== Running the main analysis");
		runOnOneInput(extraTrees, toRemoveExtraTrees,outbuffer, mainTrees, bootstraps, 
				outgroup, options, constraintTree);

		outbuffer.close();
	}

	protected Tree runOnOneInput(List<Tree> extraTrees,
			List<Tree> toRemoveExtraTrees, BufferedWriter outbuffer, List<Tree> input, 
            Iterable<Tree> bootstraps, String outgroup, Options options, List<Tree> constraintTree) {
		long startTime;
		startTime = System.currentTimeMillis();
		AbstractInference inference =
				Factory.instance.newInference(options, input, extraTrees, toRemoveExtraTrees);
		if (constraintTree != null)
			inference.setConstraintTree(constraintTree);

		inference.setup(); 

		List<Solution> solutions = inference.inferSpeciesTree();

		Logging.logTimeMessage(" CommandLine 667: ");
		Logging.log("Optimal tree inferred in "
				+ (System.currentTimeMillis() - startTime) / 1000.0D + " secs.");
		//Logging.log("Weight calculation using polytrees cumulatively took " + Polytree.time / 1000000000.0D + " secs");

		Tree st = processSolution(outbuffer, bootstraps, outgroup, inference, solutions);

		return st;
	}

	protected boolean isGeneResamplign(JSAPResult config) {
		return config.getBoolean("gene-sampling") || config.getBoolean("gene-only") ;
	}

	protected Tree processSolution(BufferedWriter outbuffer,
			Iterable<Tree> bootstraps, String outgroup,
			AbstractInference inference, List<Solution> solutions) {
		Logging.logTimeMessage(" CommandLine 684: ");
		Tree st = solutions.get(0).getTree();
		Logging.logTimeMessage(" CommandLine 690: ");
		Logging.log(st.toNewick());

		st.rerootTreeAtNode(st.getNode(outgroup));
		Trees.removeBinaryNodes((MutableTree) st);

		// TODO: MULTIND. 
		GlobalMaps.taxonNameMap.getSpeciesIdMapper().stToGt((MutableTree) st);
		inference.scoreSpeciesTreeWithGTLabels(st, false);
		GlobalMaps.taxonNameMap.getSpeciesIdMapper().gtToSt((MutableTree) st);
		Iterator<TNode> ci = (Iterator<TNode>) st.getRoot().getChildren().iterator();
		TNode c = ci.next();
		while (c.isLeaf()) c=ci.next();
		c.setParentDistance(0);

		if ((bootstraps != null) && (bootstraps.iterator().hasNext())) {
			for (Solution solution : solutions) {
				Utils.computeEdgeSupports((STITree<Double>) solution.getTree(), bootstraps);
			}
		}
		writeTreeToFile(outbuffer, solutions.get(0).getTree());

		return st;
	}

	protected List<String> readTreeFileAsString(File file)
			throws FileNotFoundException, IOException {
		String line;		
		List<String> trees = new ArrayList<String>();
		BufferedReader treeBufferReader = new BufferedReader(new FileReader(file));
		while ((line = treeBufferReader .readLine()) != null) {
			if (line.length() > 0) {
				line = line.replaceAll("\\)[^,);:]*", ")");
				trees.add(line);
			}
		}
		treeBufferReader.close();
		return trees;

	}

	protected void readInputTrees(List<Tree> trees, List<String> lines, 
			boolean rooted, boolean checkCompleteness, boolean stLablel,
			Integer minleaves, int annotation, String outgroup)
					throws FileNotFoundException, IOException {

		List<Integer> skipped = new Stack<Integer>();
		int l = 0;			
		try {
			TreeSet<String> allleaves = new TreeSet<String>();
			for (String line : lines) {
				l++;
				Set<String> previousTreeTaxa = new HashSet<String>();
				if (line.length()  == 0) {
					continue;
				}
				NewickReader nr = new NewickReader(new StringReader(line));
				if (rooted) {
					STITree<Double> gt = new STITree<Double>(true);
					nr.readTree(gt);
					if (checkCompleteness) {
						if (previousTreeTaxa.isEmpty()) {
							previousTreeTaxa.addAll(Arrays.asList(gt.getLeaves()));
						} else {
							if (!previousTreeTaxa.containsAll(Arrays.asList(gt.getLeaves()))) {
								throw new RuntimeException( "Not all trees are on the same set of taxa: "
										+ gt.getLeaves() + "\n" + previousTreeTaxa);
							}
						}
					}
					if (minleaves == null || gt.getLeafCount() >= minleaves) {
						trees.add(gt);
					} else {
						skipped.add(l);
					}
				} else {	
					MutableTree tr = nr.readTree();
					if (minleaves == null || tr.getLeafCount() >= minleaves) {
						trees.add(tr);
					} else {
						skipped.add(l);
					}
					if (outgroup != null) {
						tr.rerootTreeAtNode(tr.getNode(outgroup));
					}
					Trees.removeBinaryNodes(tr);
					if (stLablel) {
						GlobalMaps.taxonNameMap.getSpeciesIdMapper().stToGt(tr);
					}
					String[] leaves = tr.getLeaves().clone();
					if (annotation != 6) {
						for (int i = 0; i < leaves.length; i++) {
							//if (!stLablel) {
							GlobalMaps.taxonIdentifier.taxonId(leaves[i]);
							//} else {
							//   GlobalMaps.taxonNameMap.getSpeciesIdMapper().speciesId(leaves[i]);
							//}
						}
					} else{
						allleaves.addAll(Arrays.asList(leaves));
					}
				}
				if (annotation == 6) {
					for (String leaf: allleaves) {
						GlobalMaps.taxonIdentifier.taxonId(leaf);
					}
				}
			}
		} catch (ParseException e) {
			throw new RuntimeException("Failed to Parse Tree number: " + l ,e);
		}
		if (skipped.size() > 0) {
			Logging.log("Skipping the following tree(s) because they had less than " + minleaves + " leaves: \n" + skipped);
		}
	}


	protected void writeTreeToFile(BufferedWriter outbuffer, Tree t) {
		try {
			outbuffer.write(t.toStringWD()+ " \n");
			outbuffer.flush();
		} catch (IOException e) {
			Logging.log("Error when writing the species tree");
			Logging.log(e.getMessage());
			e.printStackTrace();
		}
	}

	protected void renmaeFromGTtoST(List<Tree> mainTrees, BufferedWriter outbuffer) {

		Map<String, Set<String>> newNameMap = new HashMap<String, Set<String>>();
		SpeciesMapper spm = GlobalMaps.taxonNameMap.getSpeciesIdMapper();
		if (spm.isSingleIndividual()) {
			throw new RuntimeException("You seem to already have a single-individual input; make sure you provided the mapping file using the -a option.");
		}
		for (Tree t: mainTrees) {
			MutableTree gt = (MutableTree) t;
			Map <String,Integer> used = new HashMap<String,Integer>();
			for (String leave : gt.getLeaves()) {
				TMutableNode node = gt.getNode(leave);
				String speciesName = spm.getSpeciesNameForTaxonName(leave);
				Integer i = used.getOrDefault(speciesName, 0);
				String newName = speciesName + "_" + i;
				i++;
				used.put(speciesName, i);
				if (!newNameMap.containsKey(speciesName)) {
					newNameMap.put(speciesName,  new TreeSet<String>());
				}
				newNameMap.get(speciesName).add(newName);
				node.setName(newName);
			}
			writeTreeToFile(outbuffer, gt);
		}

		for (Entry<String, Set<String>> e : newNameMap.entrySet()) {
			System.out.print(e.getKey()+" "+e.getValue().size()+" ");
			for (String g : e.getValue()) {
				System.out.print(g+" ");
			}
			System.out.println();
		}

	}


	public static void main(String[] args) throws Exception{
		CommandLine cm = new CommandLine();
		cm.process(args);
	}

}
