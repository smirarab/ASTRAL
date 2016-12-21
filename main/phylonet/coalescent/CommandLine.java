package phylonet.coalescent;

import static org.jocl.CL.CL_CONTEXT_PLATFORM;
import static org.jocl.CL.CL_DEVICE_LOCAL_MEM_SIZE;
import static org.jocl.CL.CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE;
import static org.jocl.CL.CL_DEVICE_NAME;
import static org.jocl.CL.CL_DEVICE_TYPE_ALL;
import static org.jocl.CL.CL_MEM_COPY_HOST_PTR;
import static org.jocl.CL.CL_MEM_OBJECT_IMAGE2D;
import static org.jocl.CL.CL_MEM_READ_ONLY;
import static org.jocl.CL.CL_MEM_READ_WRITE;
import static org.jocl.CL.CL_MEM_WRITE_ONLY;
import static org.jocl.CL.CL_RGBA;
import static org.jocl.CL.CL_SIGNED_INT16;
import static org.jocl.CL.CL_TRUE;
import static org.jocl.CL.clBuildProgram;
import static org.jocl.CL.clCreateBuffer;
import static org.jocl.CL.clCreateCommandQueue;
import static org.jocl.CL.clCreateContext;
import static org.jocl.CL.clCreateKernel;
import static org.jocl.CL.clCreateProgramWithSource;
import static org.jocl.CL.clEnqueueNDRangeKernel;
import static org.jocl.CL.clEnqueueReadBuffer;
import static org.jocl.CL.clEnqueueWriteBuffer;
import static org.jocl.CL.clGetDeviceIDs;
import static org.jocl.CL.clGetDeviceInfo;
import static org.jocl.CL.clGetPlatformIDs;
import static org.jocl.CL.clSetKernelArg;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentLinkedQueue;

import org.jocl.CL;
import org.jocl.Pointer;
import org.jocl.Sizeof;
import org.jocl.cl_command_queue;
import org.jocl.cl_context;
import org.jocl.cl_context_properties;
import org.jocl.cl_device_id;
import org.jocl.cl_image_desc;
import org.jocl.cl_image_format;
import org.jocl.cl_kernel;
import org.jocl.cl_mem;
import org.jocl.cl_platform_id;
import org.jocl.cl_program;

import phylonet.tree.io.NewickReader;
import phylonet.tree.io.ParseException;
import phylonet.tree.model.MutableTree;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.util.Trees;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import com.martiansoftware.jsap.SimpleJSAP;
import com.martiansoftware.jsap.Switch;
import com.martiansoftware.jsap.stringparsers.FileStringParser;

public class CommandLine{
	
	protected static String _version = "4.11.0";

    public static ConcurrentLinkedQueue<Tripartition> queue1 = new ConcurrentLinkedQueue<Tripartition>();
    public static ConcurrentLinkedQueue<Long> queue2 = new ConcurrentLinkedQueue<Long>();
    
    public static long workGroupSize = 1L<<13;
    
    //public static int SPECIES_WORD_LENGTH;
    
    private static void exitWithErr(String extraMessage, SimpleJSAP jsap) {
        System.err.println();
        System.err.println(extraMessage);
        System.err.println();
        System.err.println("Usage: java -jar astral."+_version+".jar "+ jsap.getUsage());
        System.err.println();
        System.err.println(jsap.getHelp());
        System.exit( 1 );
    }

    
    private static SimpleJSAP getJSAP() throws JSAPException {
        return new SimpleJSAP(
                "ASTRAL (version" + _version + ")",
                "species tree inference from unrooted gene trees. "
                + "The ASTRAL algorithm maximizes the number of shared quartet trees with"
                + " the collection of all gene trees. The result of this optimization problem"
                + " is statistically consistent under the multi-species coalescent model."
                + " This software can also solve MGD and MGDL problems (see options) instead of ASTRAL.",
                    
                new Parameter[] {
                    
                    new FlaggedOption("input file", 
                            FileStringParser.getParser().setMustExist(true), null, JSAP.REQUIRED, 
                            'i', "input",
                            "a file containing input gene trees in newick format. (required)"),
                            
                    new FlaggedOption( "output file",
                            FileStringParser.getParser(), null, JSAP.NOT_REQUIRED,
                            'o', "output",
                            "a filename for storing the output species tree. Defaults to outputting to stdout."),

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
                            + "4: only three alternative posterior probabilities."),
                            
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
	                        + " -k bootstraps: outputs individual bootstrap replicates to a file called [output file name].[i].bs\n"
	                        + " -k bootstraps_norun: just like -k bootstraps, but exits after outputting bootstraps.\n"
	                        + " -k searchspace_norun: outputs the search space and exits; use -k searchspace to continue the run after outputting the search space."
	                        + "When -k option is used, -o option needs to be given. "
	                        + "The file name specified using -o is used as the prefix for the name of the extra output files.").setAllowMultipleDeclarations(true),

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
	
	                new FlaggedOption("minleaves", 
	                        JSAP.INTEGER_PARSER, null, JSAP.NOT_REQUIRED, 
	                        'm', "minleaves",
	                        "Remove genes with less than specified number of leaves "),
	
	                new Switch( "duplication",
	                        'd', "dup",
	                        "Solves MGD problem. Minimizes the number duplications required to explain "
	                        + "gene trees using DynaDup algorithm (Bayzid, 2011). Note that with this option, "
	                        + "DynaDyp would be used *instead of* ASTRAL."),
	                                            
                    new Switch("exact",
                            'x', "exact",
                            "find the exact solution by looking at all clusters - recommended only for small (<18) number of taxa."),

/*                    new Switch("scoreall",
                            'y', "scoreall",
                            "score all possible species trees."),*/

                    new FlaggedOption("extraLevel",
                    		JSAP.INTEGER_PARSER, "1", JSAP.NOT_REQUIRED,
                            'p', "extraLevel",
                            "How much extra bipartitions should be added: 0, 1, or 2. "
                            + "0: adds nothing extra. "
                            + "1 (default): adds to X but not excessively (greedy resolutions). "
                            + "2: adds a potentially large number and therefore can be slow (quadratic distance-based)."),
   
                    new FlaggedOption("extra trees", 
                            FileStringParser.getParser().setMustExist(true), null, JSAP.NOT_REQUIRED, 
                            'e', "extra",
                            "provide extra trees (with gene labels) used to enrich the set of clusters searched"),

                    new FlaggedOption("extra species trees", 
                            FileStringParser.getParser().setMustExist(true), null, JSAP.NOT_REQUIRED, 
                            'f', "extra-species",
                            "provide extra trees (with species labels) used to enrich the set of clusters searched"),

                    new FlaggedOption( "duploss weight",
                            JSAP.STRING_PARSER, null, JSAP.NOT_REQUIRED,
                            'l', "duploss",
                            "Solves MGDL problem. Minimizes the number duplication and losses required"
                            + " to explain gene trees using DynaDup algorithm. Note that with this option, "
                            + "DynaDyp would be used *instead of* ASTRAL. "
                            + "Use -l 0 for standard (homomorphic) definition, and -l 1 for our new bd definition. "
                            + "Any value in between weights the impact of missing taxa somewhere between these two extremes. "
                            + "-l auto will automatically pick this weight. "), 
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
		List<Tree> mainTrees;
		List<List<String>> bstrees = new ArrayList<List<String>>();
		List<List<String>> bootstrapInputSets = new ArrayList<List<String>>();
		List<Tree> extraTrees = new ArrayList<Tree>();
		double wh = 1.0D;
		//int addExtra;
		int k = 0;
		//int annotate = 1;
		Integer minleaves = null;
        BufferedWriter outbuffer;
        Set<String> keepOptions = new HashSet<String>();
        String outfileName = null;
        
		
        jsap = getJSAP();     
        config = jsap.parse(args);  
        if ( jsap.messagePrinted() ) {
            exitWithErr("",jsap);
        }
        
        if (config.getBoolean("gene-only") && config.getFile("bootstraps") != null) {
        	exitWithErr("--gene-only and -b cannot be used together",jsap);
        }
        
        File outfile = config.getFile("output file");  
        if (outfile == null) {
            outbuffer = new BufferedWriter(new OutputStreamWriter(System.out));
        } else {
            outbuffer = new BufferedWriter(new FileWriter(outfile));
			outfileName = config.getFile("output file") == null? 
					null: config.getFile("output file").getCanonicalPath();
        }
        
        
        
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
        System.err.println("This is ASTRAL version " + _version);

        System.err.println("Gene trees are treated as " + (rooted ? "rooted" : "unrooted"));

        GlobalMaps.random = new Random(config.getLong("seed"));
        
        
        if (config.getFile("mapping file") != null) {

        	System.err.println("****** WARNNING ******\n"
        			+ "		For multi-individual data, please use the code\n"
        			+ "		available at the multiind branch of the ASTRAL github.\n"
        			+ "		This branch does not yet include our improvements for the\n"
        			+ "		multi-individual inputs.\n"
        			+ "****** END OF WARNNING ******\n");
            BufferedReader br = new BufferedReader(new FileReader(
                    config.getFile("mapping file")));

            taxonMap = new HashMap<String, String>();
            String s;
            try {
            while ((s = br.readLine()) != null) {
                s = s.trim();
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
                        System.err
                        .println("The input file is not in correct format");
                        System.err
                        .println("Any gene name can only map to one species");
                        System.exit(-1);
                    } else if (alleles.length > 1 && allele.equals(species)) {
                        System.err
                        .println("Error: The species name cannot be identical to gene names when"
                        		+ "multiple alleles exist for the same gene"+ allele);
                        System.exit(-1);
                	}
                    //System.err.println("Mapping '"+allele+"' to '"+species+"'");
                    taxonMap.put(allele, species);
                }
            }

            } catch (Exception e) {
            	throw new RuntimeException("\n** Error **: Your name mapping file looks incorrect.\n   Carefully check its format. ", e);
            }
            br.close();
        }
        
        minleaves = config.contains("minleaves")? config.getInt("minleaves"):null;
        
        try {
        	
        	//GlobalMaps.taxonIdentifier.taxonId("0");

        	mainTrees = readInputTrees(
        			readTreeFileAsString(config.getFile("input file")),
        					rooted, true, false, minleaves, 
        					config.getInt("branch annotation level"), null);			
            k = mainTrees.size();
            System.err.println(k+" trees read from " + config.getFile("input file"));
            
            GlobalMaps.taxonIdentifier.lock();        
            
            //johng23
            System.err.println("global work group size is : " + workGroupSize);
            final int platformIndex = 0;
    		final long deviceType = CL_DEVICE_TYPE_ALL;
    		final int deviceIndex = 0;

    		// Enable exceptions and subsequently omit error checks in this sample
    		CL.setExceptionsEnabled(true);

    		// Obtain the number of platforms
    		int numPlatformsArray[] = new int[1];
    		clGetPlatformIDs(0, null, numPlatformsArray);
    		int numPlatforms = numPlatformsArray[0];

    		// Obtain a platform ID
    		cl_platform_id platforms[] = new cl_platform_id[numPlatforms];
    		clGetPlatformIDs(platforms.length, platforms, null);
    		cl_platform_id platform = platforms[platformIndex];

    		// Initialize the context properties
    		cl_context_properties contextProperties = new cl_context_properties();
    		contextProperties.addProperty(CL_CONTEXT_PLATFORM, platform);

    		// Obtain the number of devices for the platform
    		int numDevicesArray[] = new int[1];
    		clGetDeviceIDs(platform, deviceType, 0, null, numDevicesArray);
    		int numDevices = numDevicesArray[0];

    		// Obtain a device ID
    		cl_device_id devices[] = new cl_device_id[numDevices];
    		clGetDeviceIDs(platform, deviceType, numDevices, devices, null);
    		for (int i=0; i<numDevices; i++)
            {
                String deviceName = getString(devices[i], CL_DEVICE_NAME);
                System.out.println("Device "+i+" of "+numDevices+": "+deviceName);
            }
    		cl_device_id device = devices[deviceIndex];


        } catch (IOException e) {
            System.err.println("Error when reading trees.");
            System.err.println(e.getMessage());
            e.printStackTrace();
            return;
        } 
        
        if (mainTrees.size() == 0) {
            System.err.println("Empty list of trees. The function exits.");
            return;
        }
                
        if (taxonMap != null) {
            GlobalMaps.taxonNameMap = new TaxonNameMap(taxonMap);
        } else if (replace != null) {   
            GlobalMaps.taxonNameMap = new TaxonNameMap (pattern, replace);
        } else {
            GlobalMaps.taxonNameMap = new TaxonNameMap();
        }
        
        Options options = newOptions(criterion, rooted, extrarooted, 
        		1.0D, 1.0D, wh, keepOptions, config, outfileName);
        
        /*Options bsoptions = newOptions(criterion, rooted, extrarooted, 
        		1.0D, 1.0D, wh, keepOptions, config);
        bsoptions.setBranchannotation(0);*/
       
        String outgroup = GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSpeciesName(0);
        
        List<String> toScore = null;
        //List<StringBuffer> toScoreStrings = null;
        if (config.getFile("score species trees") != null) {
        	toScore = readTreeFileAsString(config.getFile("score species trees"));
        } else if (false) { //config.getBoolean("scoreall")) {
            toScore = Utils.generateAllBinaryTreeStrings(
            		GlobalMaps.taxonNameMap.getSpeciesIdMapper().getAllSpeciesNames());
            /*for (String trs : toScore) {
                System.out.println(trs);
                System.out.flush();
            }*/
        }
        
        if (toScore != null ) {	
            runScore(criterion, rooted, mainTrees, extraTrees, outbuffer,
					options, outgroup, toScore);
        } else {
        
	        runInference(jsap, config, criterion, rooted, extrarooted,
					mainTrees, bstrees, bootstrapInputSets, extraTrees, k,
					minleaves, outbuffer, keepOptions, outfile, options,
					outgroup);
        }
		
	    System.err.println("ASTRAL finished in "  + 
	            (System.currentTimeMillis() - startTime) / 1000.0D + " secs");
	}


	private static void runScore(int criterion, boolean rooted,
			List<Tree> mainTrees, List<Tree> extraTrees,
			BufferedWriter outbuffer, Options options, String outgroup,
			List<String> toScore) throws FileNotFoundException, IOException {
            System.err.println("Scoring: " + toScore.size() +" trees");
            
            AbstractInference inference =
                    initializeInference(criterion, mainTrees, extraTrees, options, true);           
       		double score = Double.NEGATIVE_INFINITY;
       		List<Tree> bestTree = new ArrayList<Tree>(); 
            for (String trs : toScore) {     
            	Tree tr = readInputTrees(Arrays.asList(new String[]{trs}),
                         rooted, true, true, null, 1, false? //config.getBoolean("scoreall")? 
                        		 outgroup: null).get(0);

			double nscore = inference.scoreSpeciesTreeWithGTLabels(tr, true);
			
				if (nscore > score) {
					score = nscore;
					bestTree.clear();
					bestTree.add(tr);
				} else if (nscore == score) {
					bestTree.add(tr);
				}
				
			if (!GlobalMaps.taxonNameMap.getSpeciesIdMapper().isSingleIndividual()) {
				System.err.println("Scored tree with gene names:\n"+tr.toNewickWD());
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
        

	private static void runInference(SimpleJSAP jsap, JSAPResult config,
			int criterion, boolean rooted, boolean extrarooted,
			List<Tree> mainTrees, List<List<String>> bstrees,
			List<List<String>> bootstrapInputSets, List<Tree> extraTrees,
			int k, Integer minleaves, BufferedWriter outbuffer,
			Set<String> keepOptions, File outfile, Options options,
			String outgroup) throws JSAPException, IOException,
			FileNotFoundException {
		System.err.println("All output trees will be *arbitrarily* rooted at "+outgroup);
        //SPECIES_WORD_LENGTH = GlobalMaps.taxonIdentifier.taxonCount()/64 + 1;

        if (config.getStringArray("keep") != null && config.getStringArray("keep").length != 0) {
			if (options.getOutputFile() == null) {
        		throw new JSAPException("When -k option is used, -o is also needed.");
        	}
        	for (String koption : config.getStringArray("keep")) {
        		if ("completed".equals(koption) ||
        			"bootstraps".equals(koption) ||
        			"bootstraps_norun".equals(koption)||
        			"searchspace_norun".equals(koption)||
        			"searchspace".equals(koption)) {
        			keepOptions.add(koption);
        		} else {
        			throw new JSAPException("-k "+koption+" not recognized.");
        		}
        	}
        }
        
        try {

            if (config.getFile("extra trees") != null) {
            	extraTrees = readInputTrees(
                	readTreeFileAsString(config.getFile("extra trees")), 
                        extrarooted, true, false, null, 1, null);
                System.err.println(extraTrees.size() + " extra trees read from "
                        + config.getFile("extra trees"));
            }
            
            if (config.getFile("extra species trees") != null) {
            	extraTrees = readInputTrees(
                	readTreeFileAsString(config.getFile("extra species trees")), 
                        extrarooted, true, true, null, 1, null);
                System.err.println(extraTrees.size() + " extra trees read from "
                        + config.getFile("extra trees"));
            }
            
        } catch (IOException e) {
            System.err.println("Error when reading extra trees.");
            System.err.println(e.getMessage());
            e.printStackTrace();
		    System.exit(1);
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
            System.err.println("Error when reading bootstrap trees.");
            System.err.println(e.getMessage());
            e.printStackTrace();
		    System.exit(1);
        } 

	
		if (config.getFile("bootstraps") != null || config.getBoolean("gene-only")) {
	        System.err.println("Bootstrapping with seed "+config.getLong("seed"));
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
			        }
			        else {   		        
    		            for (List<String> gene : bstrees) {
    		                input.add(gene.get(i));
    		            }
    		        }
		        } catch (IndexOutOfBoundsException e) {
		            exitWithErr("Error: You seem to have asked for "+config.getInt("replicates")+
		                    " but only "+ i +" replicates could be created.\n" + 
		                    " Note that for gene resampling, you need more input bootstrap" +
		                    " replicates than the number of species tree replicates.", jsap);
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
		    	System.err.println("bootstrap files written to files "+ outfile + ( "." + 0 + ".bs" ) + 
		    			" to "+outfile + ( "." + config.getInt("replicates") + ".bs" ));
		    }
		    if (keepOptions.contains("bootstraps_norun")) {
		    	System.err.println("Exiting after outputting the bootstrap files");
		    	System.exit(0);
		    }
		}

        int j = 0;
        List<Tree> bootstraps = new ArrayList<Tree>();
		for ( List<String> input : bootstrapInputSets) {  
	        System.err.println("\n======== Running bootstrap replicate " + j++);
	        bootstraps.add(runOnOneInput(criterion, 
	                 extraTrees, outbuffer, 
                    readInputTrees(input, rooted, false, false, minleaves,
                    		config.getInt("branch annotation level"), null),
                    null, outgroup, options));
		}
	    
		if (bootstraps != null && bootstraps.size() != 0) {
            STITree<Double> cons = (STITree<Double>) Utils.greedyConsensus(bootstraps,false, GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTTaxonIdentifier());
            cons.rerootTreeAtNode(cons.getNode(outgroup));
			Trees.removeBinaryNodes(cons);
            Utils.computeEdgeSupports(cons, bootstraps);
            writeTreeToFile(outbuffer, cons);
		}
		
        System.err.println("\n======== Running the main analysis");
        runOnOneInput(criterion, extraTrees, outbuffer, mainTrees, bootstraps, 
                outgroup, options);
           
		outbuffer.close();
	}


	private static int getSpeciesWordLength(){
		return (GlobalMaps.taxonIdentifier.taxonCount()/64 + 1);
	}
    private static Tree runOnOneInput(int criterion, List<Tree> extraTrees,
    		BufferedWriter outbuffer, List<Tree> input, 
            Iterable<Tree> bootstraps, String outgroup, Options options) {
        long startTime;
        startTime = System.currentTimeMillis();
        
        AbstractInference inference = initializeInference(criterion, input, extraTrees, options, true);
        inference.queue2 = queue2;
        
        inference.setup(); 
        
        AbstractInferenceNoCalculations inferenceNoCalc = new WQInferenceNoCalculations((AbstractInference) inference.semiDeepCopy());

        inferenceNoCalc.queue1 = queue1;
        int counter = 0;
        long[] allArray = new long[ ((WQDataCollection)inference.dataCollection).treeAllClusters.size()*getSpeciesWordLength()];
        for(int i = 0; i < ((WQDataCollection)inference.dataCollection).treeAllClusters.size(); i++) {
        	for(int j = getSpeciesWordLength() - 1 ; j >= 0; j--)
        		allArray[counter++] = ((WQDataCollection)inference.dataCollection).treeAllClusters.get(i).getBitSet().words[j];
        }
        /*int[] geneTreesAsInts = new int[((WQDataCollection)inference.dataCollection).geneTreesAsInts.length];
        for(int i = 0; i < geneTreesAsInts.length; i++) {
        	geneTreesAsInts[i] = ((WQDataCollection)inference.dataCollection).geneTreesAsInts[i];
        }*/
        TurnTaskToScores threadgpu = new TurnTaskToScores(inference, queue1, queue2, ((WQWeightCalculator)inference.weightCalculator).geneTreesAsInts(), allArray);
		WriteTaskToQueue thread1 = new WriteTaskToQueue(inferenceNoCalc, threadgpu);

		(new Thread(thread1)).start();
		(new Thread(threadgpu)).start();

        List<Solution> solutions = inference.inferSpeciesTree();
        
        System.err.println("Optimal tree inferred in "
        		+ (System.currentTimeMillis() - startTime) / 1000.0D + " secs.");
        
        Tree st = processSolution(outbuffer, bootstraps, outgroup, inference, solutions);
        
        return st;
    }

    private static boolean isGeneResamplign(JSAPResult config) {
    	return config.getBoolean("gene-sampling") || config.getBoolean("gene-only") ;
    }

	private static Tree processSolution(BufferedWriter outbuffer,
			Iterable<Tree> bootstraps, String outgroup,
			AbstractInference inference, List<Solution> solutions) {
        Tree st = solutions.get(0)._st;
        
        System.err.println(st.toNewick());
        
        st.rerootTreeAtNode(st.getNode(outgroup));
        
		Trees.removeBinaryNodes((MutableTree) st);

		GlobalMaps.taxonNameMap.getSpeciesIdMapper().stToGt((MutableTree) st);
		
		inference.scoreSpeciesTreeWithGTLabels(st, false);
		
		GlobalMaps.taxonNameMap.getSpeciesIdMapper().gtToSt((MutableTree) st);
		
        if ((bootstraps != null) && (bootstraps.iterator().hasNext())) {
            for (Solution solution : solutions) {
                Utils.computeEdgeSupports((STITree<Double>) solution._st, bootstraps);
            }
        }
        writeTreeToFile(outbuffer, solutions.get(0)._st);
        
        return st;
    }
    public static class TurnTaskToScores implements Runnable {
		public ConcurrentLinkedQueue<Long> queue2;
		public ConcurrentLinkedQueue<Tripartition> queue1;
		public AbstractInference inference;
		public long[] tripartition1;
		public long[] tripartition2;
		public long[] tripartition3;
		
		public long[] all;
		public int tripCounter = 0;
		public boolean done = false;
		public GPUCall gpu;
		
		public final boolean pjohng23 = false;
		public TurnTaskToScores(AbstractInference inf, ConcurrentLinkedQueue<Tripartition> queue1, ConcurrentLinkedQueue<Long> queue2, int[] geneTreeAsInts, long[] all) {
			this.inference = inf;
			this.queue1 = queue1;
			this.queue2 = queue2;
			this.all = all;
			tripartition1 = new long[(int)(getSpeciesWordLength() * workGroupSize)];
			tripartition2 = new long[(int)(getSpeciesWordLength() * workGroupSize)];
			tripartition3 = new long[(int)(getSpeciesWordLength() * workGroupSize)];
			
			gpu = new GPUCall(geneTreeAsInts, all, tripartition1, tripartition2, tripartition3, inference, pjohng23);
		}

		public void run() {
			
			Tripartition task = null;
			((WQWeightCalculator)inference.weightCalculator).lastTime = System.currentTimeMillis();
			while(!done || !queue1.isEmpty()) {
				if(!queue1.isEmpty()) {

					task = queue1.remove();
					
					for(int i = getSpeciesWordLength() - 1; i >= 0; i--)
						//tripartition1[tripCounter * SPECIES_WORD_LENGTH + SPECIES_WORD_LENGTH - i - 1] = task.trip.cluster1.getBitSet().words[i];
						tripartition1[(getSpeciesWordLength() - i - 1) * (int)workGroupSize + tripCounter]=task.cluster1.getBitSet().words[i];
					for(int i = getSpeciesWordLength() - 1; i >= 0; i--)
						//tripartition2[tripCounter * SPECIES_WORD_LENGTH + SPECIES_WORD_LENGTH - i - 1] = task.trip.cluster2.getBitSet().words[i];
						tripartition2[(getSpeciesWordLength() - i - 1) * (int)workGroupSize + tripCounter]=task.cluster2.getBitSet().words[i];
					for(int i = getSpeciesWordLength() - 1; i >= 0; i--)
						//tripartition3[tripCounter * SPECIES_WORD_LENGTH + SPECIES_WORD_LENGTH - i - 1] = task.trip.cluster3.getBitSet().words[i];
						tripartition3[(getSpeciesWordLength() - i - 1) * (int)workGroupSize + tripCounter]=task.cluster3.getBitSet().words[i];
					tripCounter++;
					if(tripCounter == workGroupSize) {
						gpu.compute(workGroupSize);
						tripCounter = 0;
						for(int i = 0; i < gpu.weightArray.length; i++) {
							queue2.add(gpu.weightArray[i]);
						}
					}

			
				}
			}
			
			gpu.compute(tripCounter);

			for(int i = 0; i < tripCounter; i++) {
				queue2.add(gpu.weightArray[i]);
			}
			
			if(pjohng23) {
    				clEnqueueReadBuffer(gpu.commandQueue, gpu.d_profile, CL_TRUE, 0L, Sizeof.cl_long * gpu.profile.length, Pointer.to(gpu.profile), 0, null, null);
				long total = 0;	
				for(int i = 0; i < 4; i++) {
					total += gpu.profile[i];
				}
				System.out.println("intersecting with the all arrays takes: " + 4*(double)gpu.profile[0]/1000000000);
				System.out.println("adding numbers to the stack takes: " + 4*(double)gpu.profile[1]/1000000000);
				System.out.println("calculating weight of a tripartition takes: " + 4*(double)gpu.profile[2]/1000000000);
				System.out.println("calculating weight of a polytomy takes: " + 4*(double)gpu.profile[3]/1000000000);
			}
			inference.weightCalculator.done = true;
			
		}
		
		
	}
    public static class GPUCall {
    	public long[] tripartitions1;
    	public long[] tripartitions2;
    	public long[] tripartitions3;
    	
    	public int[] geneTreesAsInts;
    	public short[] geneTreesAsShorts;
    	public long[] allArray;
    	public long[] weightArray;
    	public short[] stack;
    	public long[] profile;
    	public AbstractInference inference;
    	
	public boolean p;
    	private cl_context context;
    	private cl_command_queue commandQueue;
    	private cl_kernel kernel;
    	private cl_mem d_geneTreesAsInts;
    	private cl_mem[] d_geneTreesAsIntsConst;
    	private cl_mem d_tripartitions1;
    	private cl_mem d_tripartitions2;
    	private cl_mem d_tripartitions3;
    	private cl_mem d_allArray;
    	private cl_mem d_weightArray;
    	private cl_mem d_stack;
    	private cl_mem d_profile;
    	public GPUCall (int[] geneTreesAsInts, long[] all, long[] trip1, long[] trip2, long[] trip3, AbstractInference inference, boolean p) {
		this.p = p;
    		this.geneTreesAsInts = geneTreesAsInts;
    		geneTreesAsShorts = new short[(geneTreesAsInts.length/4 + 1) * 4];
    		for(int i = 0; i < geneTreesAsInts.length; i++) {
    			geneTreesAsShorts[i] = (short)geneTreesAsInts[i];
    			if(geneTreesAsInts[i] == Integer.MIN_VALUE)
    				geneTreesAsShorts[i] = Short.MIN_VALUE;
    		}
    		allArray = all;
    		tripartitions1 = trip1;
    		tripartitions2 = trip2;
    		tripartitions3 = trip3;
    		weightArray = new long[trip1.length / getSpeciesWordLength()];
    		this.inference = inference;
    		initCL();

    		prepare();

    	}
    	public void initCL() {
    		final int platformIndex = 0;
    		final long deviceType = CL_DEVICE_TYPE_ALL;
    		final int deviceIndex = 0;

    		// Enable exceptions and subsequently omit error checks in this sample
    		CL.setExceptionsEnabled(true);

    		// Obtain the number of platforms
    		int numPlatformsArray[] = new int[1];
    		clGetPlatformIDs(0, null, numPlatformsArray);
    		int numPlatforms = numPlatformsArray[0];

    		// Obtain a platform ID
    		cl_platform_id platforms[] = new cl_platform_id[numPlatforms];
    		clGetPlatformIDs(platforms.length, platforms, null);
    		cl_platform_id platform = platforms[platformIndex];

    		// Initialize the context properties
    		cl_context_properties contextProperties = new cl_context_properties();
    		contextProperties.addProperty(CL_CONTEXT_PLATFORM, platform);

    		// Obtain the number of devices for the platform
    		int numDevicesArray[] = new int[1];
    		clGetDeviceIDs(platform, deviceType, 0, null, numDevicesArray);
    		int numDevices = numDevicesArray[0];

    		// Obtain a device ID
    		cl_device_id devices[] = new cl_device_id[numDevices];
    		clGetDeviceIDs(platform, deviceType, numDevices, devices, null);
    		for (int i=0; i<numDevices; i++)
            {
                String deviceName = getString(devices[i], CL_DEVICE_NAME);
                System.out.println("Device "+i+" of "+numDevices+": "+deviceName);
            }
    		cl_device_id device = devices[deviceIndex];

    		// Create a context for the selected device
    		context = clCreateContext(contextProperties, 1, new cl_device_id[] { device }, null, null, null);

    		// Create a command-queue for the selected device
    		commandQueue = clCreateCommandQueue(context, device, 0, null);

    		//moved to here to edit source
    		cl_image_format geneTreesImageFormat = new cl_image_format();
    		geneTreesImageFormat.image_channel_data_type = CL_SIGNED_INT16;
    		geneTreesImageFormat.image_channel_order = CL_RGBA;
    		
    		cl_image_desc geneTreesImageDesc = new cl_image_desc();
    		geneTreesImageDesc.image_type = CL_MEM_OBJECT_IMAGE2D;
    		geneTreesImageDesc.image_width = 4*((int)Math.sqrt((double)geneTreesAsShorts.length)/16);
    		geneTreesImageDesc.image_height = geneTreesAsShorts.length/geneTreesImageDesc.image_width+1;
    		geneTreesImageDesc.image_depth = 0;
    		geneTreesImageDesc.image_array_size = 0;
    		geneTreesImageDesc.image_row_pitch = 0;
    		geneTreesImageDesc.image_slice_pitch = 0;
    		geneTreesImageDesc.num_mip_levels = 0;
    		geneTreesImageDesc.num_samples = 0;
    		geneTreesImageDesc.buffer = null;
    		
    		// getting the tree height
    		int treeheight = ((WQWeightCalculator)inference.weightCalculator).maxHeight();
    		System.out.println("TREE HEIGHT IS: " + treeheight);
    		// Program Setup
    		String source = readFile(getClass().getResourceAsStream("/phylonet/coalescent/calculateWeight.cl"));
    		source = source.replaceAll("SPECIES_WORD_LENGTH - 1", Long.toString(getSpeciesWordLength() - 1));
    		source = source.replaceAll("SPECIES_WORD_LENGTH", Long.toString(getSpeciesWordLength()));
    		source = source.replaceAll("LONG_BIT_LENGTH", "64");
//    		source = source.replaceAll("(STACK_SIZE + 1) * 3", Integer.toString((treeheight + 1) * 3));
//    		source = source.replaceAll("(STACK_SIZE + 2) * 3", Integer.toString((treeheight + 2) * 3));
//    		source = source.replaceAll("(STACK_SIZE + 1) * 2", Integer.toString((treeheight + 1) * 2));
//    		source = source.replaceAll("(STACK_SIZE + 2) * 2", Integer.toString((treeheight + 2) * 2));
//    		source = source.replaceAll("(STACK_SIZE + 1)", Integer.toString(treeheight + 1));
//    		source = source.replaceAll("(STACK_SIZE + 2)", Integer.toString(treeheight + 2));
    		source = source.replaceAll("STACK_SIZE", Integer.toString(treeheight));
    		source = source.replaceAll("TAXON_SIZE", Integer.toString(GlobalMaps.taxonIdentifier.taxonCount()));
    		source = source.replaceAll("INT_MIN", "SHRT_MIN");
    		source = source.replaceAll("WORK_GROUP_SIZE", Long.toString(workGroupSize));
    		source = source.replaceAll("IMAGE_WIDTH", Long.toString(geneTreesImageDesc.image_width));
    		
    		long[] localmemsize = new long[1];
    		long[] constmemsize = new long[1];
    		clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, 8, Pointer.to(localmemsize), null);
    		clGetDeviceInfo(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, 8, Pointer.to(constmemsize), null);
    		System.out.println("LOCAL MEMORY SIZE IS: " + localmemsize[0]);
    		System.out.println("CONSTANT MEMORY SIZE IS: " + constmemsize[0]);
    		// Create the program
    		cl_program cpProgram = clCreateProgramWithSource(context, 1, new String[] { source }, null, null);

    		// Build the program
    		if(p)
    			clBuildProgram(cpProgram, 0, null, "-cl-opt-disable", null, null);
		else
    			clBuildProgram(cpProgram, 0, null, "-cl-mad-enable -cl-strict-aliasing", null, null);

    		// Create the kernel
    		kernel = clCreateKernel(cpProgram, "calcWeight", null);

    		// Create the memory object which will be filled with the
    		// pixel data
    		
    		d_geneTreesAsInts = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Sizeof.cl_short * geneTreesAsInts.length, Pointer.to(geneTreesAsShorts), null);
    		d_allArray = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, Sizeof.cl_long * allArray.length,
    				Pointer.to(allArray), null);
    		d_tripartitions1 = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, Sizeof.cl_long * tripartitions1.length,
    				Pointer.to(tripartitions1), null);
    		d_tripartitions2 = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, Sizeof.cl_long * tripartitions2.length,
    				Pointer.to(tripartitions2), null);
    		d_tripartitions3 = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, Sizeof.cl_long * tripartitions3.length,
    				Pointer.to(tripartitions3), null);
    		
    		//stack = new short[(int) (Sizeof.cl_ushort * 3 * (2 + treeheight) * workGroupSize)];
    		d_stack = clCreateBuffer(context, CL_MEM_READ_WRITE , Sizeof.cl_ushort * 3 * (2 + treeheight) * workGroupSize,
    				null, null);
    		d_weightArray = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, Sizeof.cl_long * weightArray.length,
    				Pointer.to(weightArray), null);
		if(p){
    			profile = new long[20];
    			d_profile = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR, Sizeof.cl_long * profile.length,
    				Pointer.to(profile), null);
    		}
		//d_c = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, Sizeof.cl_int,
    		//		Pointer.to(c), null);

    	}


    	private String readFile(InputStream inputStream) {
    		try {
    			BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
    			StringBuffer sb = new StringBuffer();
    			String line = null;
    			while (true) {
    				line = br.readLine();
    				if (line == null) {
    					break;
    				}
    				sb.append(line).append("\n");
    			}
    			return sb.toString();
    		} catch (IOException e) {
    			e.printStackTrace();
    			System.exit(1);
    			return null;
    		}
    	}

    	public void prepare() {

    		
    		clSetKernelArg(kernel, 0, Sizeof.cl_mem, Pointer.to(d_geneTreesAsInts));
    		clSetKernelArg(kernel, 1, Sizeof.cl_int, Pointer.to(new int[] {geneTreesAsInts.length}));
    		clSetKernelArg(kernel, 2, Sizeof.cl_mem, Pointer.to(d_allArray));
    		clSetKernelArg(kernel, 3, Sizeof.cl_mem, Pointer.to(d_tripartitions1));
    		clSetKernelArg(kernel, 4, Sizeof.cl_mem, Pointer.to(d_tripartitions2));
    		clSetKernelArg(kernel, 5, Sizeof.cl_mem, Pointer.to(d_tripartitions3));
       		clSetKernelArg(kernel, 6, Sizeof.cl_mem, Pointer.to(d_weightArray));
       		clSetKernelArg(kernel, 7, Sizeof.cl_mem, Pointer.to(d_stack));
		if(p)
       			clSetKernelArg(kernel, 8, Sizeof.cl_mem, Pointer.to(d_profile));

//     		clSetKernelArg(kernel, 8, Sizeof.cl_mem, Pointer.to(d_c));
       		
//    		clSetKernelArg(kernel, 7, Sizeof.cl_long * workGroupSize * SPECIES_WORD_LENGTH, null);
//    		clSetKernelArg(kernel, 8, Sizeof.cl_long * workGroupSize * SPECIES_WORD_LENGTH, null);
//    		clSetKernelArg(kernel, 9, Sizeof.cl_long * workGroupSize * SPECIES_WORD_LENGTH, null);

    	}
    	
    	public void compute(long workSize) {
    		
    		// Set work size and execute the kernel
    		
    		clEnqueueWriteBuffer(commandQueue, d_tripartitions1, CL_TRUE, 0L, Sizeof.cl_long * tripartitions1.length, Pointer.to(tripartitions1), 0,
    				null, null);
    		clEnqueueWriteBuffer(commandQueue, d_tripartitions2, CL_TRUE, 0L, Sizeof.cl_long * tripartitions2.length, Pointer.to(tripartitions2), 0,
    				null, null);
    		clEnqueueWriteBuffer(commandQueue, d_tripartitions3, CL_TRUE, 0L, Sizeof.cl_long * tripartitions3.length, Pointer.to(tripartitions3), 0,
    				null, null);
//    		if(workSize >= 32 && workSize % 32 == 0) {
//    			clEnqueueNDRangeKernel(commandQueue, kernel, 1, null, new long[]{workSize}, new long[]{32L}, 0, null, null);
//    		}
//    		else	
    			clEnqueueNDRangeKernel(commandQueue, kernel, 1, null, new long[]{workSize}, null, 0, null, null);	

    		clEnqueueReadBuffer(commandQueue, d_weightArray, CL_TRUE, 0L, Sizeof.cl_long * weightArray.length, Pointer.to(weightArray), 0, null, null);
    		
    	}
    }
    public static class WriteTaskToQueue implements Runnable {
    	AbstractInferenceNoCalculations inf;
    	TurnTaskToScores threadgpu;
		public WriteTaskToQueue(AbstractInferenceNoCalculations inf, TurnTaskToScores threadgpu) {
			this.inf = inf;
			this.threadgpu = threadgpu;
		}
		public void run() {
			// TODO Auto-generated method stub
			inf.inferSpeciesTree();
			threadgpu.done = true;
		}
    	
    }
	private static String getString(cl_device_id device, int paramName)
    {
        long size[] = new long[1];
        clGetDeviceInfo(device, paramName, 0, null, size);
        byte buffer[] = new byte[(int)size[0]];
        clGetDeviceInfo(device, paramName, 
            buffer.length, Pointer.to(buffer), null);
        return new String(buffer, 0, buffer.length-1);
    }
  
	/*private static Tree processSolution(BufferedWriter outbuffer,
			Iterable<Tree> bootstraps, String outgroup,
			AbstractInference inference, List<Solution> solutions) {
		Tree st = solutions.get(0)._st;
        
        System.err.println(st.toNewick());
        
        st.rerootTreeAtNode(st.getNode(outgroup));
		Trees.removeBinaryNodes((MutableTree) st);
   
		GlobalMaps.taxonNameMap.getSpeciesIdMapper().stToGt((MutableTree) st);
		
		inference.scoreSpeciesTree(st, false);
		
		GlobalMaps.taxonNameMap.getSpeciesIdMapper().gtToSt((MutableTree) st);
		
        if ((bootstraps != null) && (bootstraps.iterator().hasNext())) {
            for (Solution solution : solutions) {
                Utils.computeEdgeSupports((STITree<Double>) solution._st, bootstraps);
            }
        }
        writeTreeToFile(outbuffer, solutions.get(0)._st);
		return st;
	}*/
    static private Options newOptions(int criterion, boolean rooted,
            boolean extrarooted, 
            double cs, double cd, double wh, 
            Set<String> keepOptions, JSAPResult config, String outfileName) {
    	Options options = new Options(
    			rooted, extrarooted, 
    			config.getBoolean("exact"), 
    			criterion > 0, 1, 
    			config.getInt("extraLevel"),
    			keepOptions.contains("completed"), 
    			keepOptions.contains("searchspace_norun") || keepOptions.contains("searchspace"), 
    			!keepOptions.contains("searchspace_norun"),
    			config.getInt("branch annotation level"), 
    			config.getDouble("lambda"),
    			outfileName);
    	options.setDLbdWeigth(wh); 
    	options.setCS(cs);
    	options.setCD(cd);
    	
    	
    	return options;
    }

    private static AbstractInference initializeInference(int criterion, 
            List<Tree> trees, List<Tree> extraTrees,
            Options options, boolean calculations) {
        AbstractInference inference;		
		if (criterion == 1 || criterion == 0) {
			inference = new DLInference(options, 
					trees, extraTrees);			
		} else if (criterion == 2) {
			if(calculations)
				inference = new WQInference(options, trees, 
					extraTrees );
			else {
				inference = new WQInferenceNoCalculations(options, trees, 
						extraTrees );
			}
		} else {
			throw new RuntimeException("criterion not set?");
		}		
        return inference;
    }

    private static List<String> readTreeFileAsString(File file)
    				throws FileNotFoundException, IOException {
    	String line;		
    	List<String> trees = new ArrayList<String>();
    	BufferedReader treeBufferReader = new BufferedReader(new FileReader(file));
		while ((line = treeBufferReader .readLine()) != null) {
    		if (line.length() > 0) {
    			line = line.replaceAll("\\)[^,);]*", ")");
    			trees.add(line);
    		}
    	}
    	treeBufferReader.close();
    	return trees;

    }

    private static List<Tree> readInputTrees(List<String> lines, 
    		boolean rooted, boolean checkCompleteness, boolean stLablel,
    		Integer minleaves, int annotation, String outgroup)
    				throws FileNotFoundException, IOException {
    	List<Tree> trees = new ArrayList<Tree>();
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
    						previousTreeTaxa.addAll(Arrays.asList(gt
    								.getLeaves()));
    					} else {
    						if (!previousTreeTaxa.containsAll(Arrays.asList(gt
    								.getLeaves()))) {
    							throw new RuntimeException(
    									"Not all trees are on the same set of taxa: "
    											+ gt.getLeaves() + "\n"
    											+ previousTreeTaxa);
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
    	        		Trees.removeBinaryNodes(tr);
                	}
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
    		System.err.println("Skipping the following tree(s) because they had less than " + minleaves+" leaves: \n" + skipped);
    	}
    	return trees;
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
    
}
