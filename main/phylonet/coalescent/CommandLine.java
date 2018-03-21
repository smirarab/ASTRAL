package phylonet.coalescent;

import static org.jocl.CL.CL_CONTEXT_PLATFORM;
import static org.jocl.CL.CL_DEVICE_NAME;
import static org.jocl.CL.CL_DEVICE_TYPE_ALL;
import static org.jocl.CL.CL_DEVICE_VENDOR;

import static org.jocl.CL.clCreateContext;
import static org.jocl.CL.clGetDeviceIDs;
import static org.jocl.CL.clGetDeviceInfo;
import static org.jocl.CL.clGetPlatformIDs;

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
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;
import java.util.concurrent.Executors;

import org.jocl.CL;
import org.jocl.Pointer;
import org.jocl.cl_context_properties;
import org.jocl.cl_device_id;
import org.jocl.cl_platform_id;

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

public class CommandLine {
	protected static String _version = "5.12.0";

	protected static SimpleJSAP jsap;

	private static void exitWithErr(String extraMessage) {
		System.err.println();
		System.err.println(extraMessage);
		System.err.println();
		System.err.println("Usage: java -jar astral." + _version + ".jar "
				+ jsap.getUsage());
		System.err.println();
		System.err.println(jsap.getHelp());
		System.exit(1);
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
						new Switch("cpu only", 'C', "cpu-only"),

						new FlaggedOption("cpu threads", JSAP.INTEGER_PARSER,
								"-1", JSAP.NOT_REQUIRED, 'T', "cpu-threads"),

						new FlaggedOption("input file", FileStringParser
								.getParser().setMustExist(true), null,
								JSAP.REQUIRED, 'i', "input",
								"a file containing input gene trees in newick format. (required)"),

						new FlaggedOption(
								"output file",
								FileStringParser.getParser(),
								null,
								JSAP.NOT_REQUIRED,
								'o',
								"output",
								"a filename for storing the output species tree. Defaults to outputting to stdout."),

						new FlaggedOption(
								"score species trees",
								FileStringParser.getParser().setMustExist(true),
								null, JSAP.NOT_REQUIRED, 'q', "score-tree",
								"score the provided species tree and exit"),

						new FlaggedOption(
								"branch annotation level",
								JSAP.INTEGER_PARSER,
								"3",
								JSAP.NOT_REQUIRED,
								't',
								"branch-annotate",
								"How much annotations should be added to each branch: 0, 1, or 2. \n"
										+ "0: no annotations. \n"
										+ "1: only the quartet support for the main resolution. \n"
										+ "2: full annotation (quartet support, quartet frequency, and posterior probability for all three alternatives, "
										+ "plus total number of quartets around the branch and effective number of genes).\n"
										+ "3 (default): only the posterior probability for the main resolution.\n"
										+ "4: three alternative posterior probabilities.\n"
										+ "8: three alternative quartet scores.\n"
										+ "10: p-values of a polytomy null hypothesis test (arxiv: 1708.08916)."),

						new FlaggedOption(
								"bootstraps",
								FileStringParser.getParser().setMustExist(true),
								null,
								JSAP.NOT_REQUIRED,
								'b',
								"bootstraps",
								"perform multi-locus bootstrapping using input bootstrap replicate files (use --rep to change the number of replications). "
										+ "The file given with this option should have a list of the gene tree bootstrap files, one per line, and each line corresponding to one gene. "
										+ "By default performs site-only resampling, but gene/site resampling can also be used. "),

						new FlaggedOption("replicates", JSAP.INTEGER_PARSER,
								"100", JSAP.NOT_REQUIRED, 'r', "reps",
								"Set the number of bootstrap replicates done in multi-locus bootstrapping. "),

						new FlaggedOption("seed", JSAP.LONG_PARSER, "692",
								JSAP.NOT_REQUIRED, 's', "seed",
								"Set the seed number used in multi-locus bootstrapping. "),

						new Switch(
								"gene-sampling",
								'g',
								"gene-resampling",
								"perform gene tree resampling in addition to site resampling. Useful only with the -b option."),

						new Switch(
								"gene-only",
								JSAP.NO_SHORTFLAG,
								"gene-only",
								"perform bootstrapping but only with gene tree resampling. Should not be used with the -b option."),

						new FlaggedOption(
								"keep",
								JSAP.STRING_PARSER,
								null,
								JSAP.NOT_REQUIRED,
								'k',
								"keep",
								" -k completed: outputs completed gene trees (i.e. after adding missing taxa) to a file called [output file name].completed_gene_trees.\n"
										+ " -k bootstraps: outputs individual bootstrap replicates to a file called [output file name].[i].bs\n"
										+ " -k bootstraps_norun: just like -k bootstraps, but exits after outputting bootstraps.\n"
										+ " -k searchspace_norun: outputs the search space and exits; use -k searchspace to continue the run after outputting the search space."
										+ "When -k option is used, -o option needs to be given. "
										+ "The file name specified using -o is used as the prefix for the name of the extra output files.")
								.setAllowMultipleDeclarations(true),

						new FlaggedOption(
								"lambda",
								JSAP.DOUBLE_PARSER,
								"0.5",
								JSAP.NOT_REQUIRED,
								'c',
								"lambda",
								"Set the lambda parameter for the Yule prior used in the calculations"
										+ " of branch lengths and posterior probabilities. Set to zero to get ML branch "
										+ "lengths instead of MAP."
										+ " Higher values tend to shorten estimated branch lengths and very"
										+ " high values can give inaccurate results (or even result in underflow)."),

						new FlaggedOption(
								"mapping file",
								FileStringParser.getParser().setMustExist(true),
								null,
								JSAP.NOT_REQUIRED,
								'a',
								"namemapfile",
								"a file containing the mapping between names in gene tree and names in the species tree. "
										+ "The mapping file has one line per species, with one of two formats:\n"
										+ " species: gene1,gene2,gene3,gene4\n"
										+ " species 4 gene1 gene2 gene3 gene4\n"),

						new FlaggedOption("minleaves", JSAP.INTEGER_PARSER,
								null, JSAP.NOT_REQUIRED, 'm', "minleaves",
								"Remove genes with less than specified number of leaves "),

						new FlaggedOption(
								"samplingrounds",
								JSAP.INTEGER_PARSER,
								null,
								JSAP.NOT_REQUIRED,
								JSAP.NO_SHORTFLAG,
								"samplingrounds",
								"For multi-individual datasets, perform these many rounds of individual sampling for"
										+ " building the set X. The program"
										+ " automatically picks this parameter if not provided or if below one."),

						new FlaggedOption("gene repetition",
								JSAP.INTEGER_PARSER, "1", JSAP.NOT_REQUIRED,
								'w', "generepeat",
								"the number of trees sampled for each locus. "),

						new FlaggedOption(
								"polylimit",
								JSAP.INTEGER_PARSER,
								null,
								JSAP.NOT_REQUIRED,
								JSAP.NO_SHORTFLAG,
								"polylimit",
								"Sets a limit for size of polytomies in greedy consensus trees where O(n) number"
										+ " of new  resolutions are added. ASTRAL-III sets automatic limits to guarantee polynomial"
										+ " time running time."),

						new Switch(
								"duplication",
								JSAP.NO_SHORTFLAG,
								"dup",
								"Solves MGD problem. Minimizes the number duplications required to explain "
										+ "gene trees using DynaDup algorithm (Bayzid, 2011). Note that with this option, "
										+ "DynaDyp would be used *instead of* ASTRAL."),

						new Switch(
								"exact",
								'x',
								"exact",
								"find the exact solution by looking at all clusters - recommended only for small (<18) number of taxa."),

						/*
						 * new Switch("scoreall", 'y', "scoreall",
						 * "score all possible species trees."),
						 */

						new FlaggedOption(
								"extraLevel",
								JSAP.INTEGER_PARSER,
								"1",
								JSAP.NOT_REQUIRED,
								'p',
								"extraLevel",
								"How much extra bipartitions should be added: 0, 1, or 2. "
										+ "0: adds nothing extra. "
										+ "1 (default): adds to X but not excessively (greedy resolutions). "
										+ "2: adds a potentially large number and therefore can be slow (quadratic distance-based)."),

						new FlaggedOption(
								"extra trees",
								FileStringParser.getParser().setMustExist(true),
								null,
								JSAP.NOT_REQUIRED,
								'e',
								"extra",
								"provide extra trees (with gene labels) used to enrich the set of clusters searched"),

						new FlaggedOption(
								"extra species trees",
								FileStringParser.getParser().setMustExist(true),
								null,
								JSAP.NOT_REQUIRED,
								'f',
								"extra-species",
								"provide extra trees (with species labels) used to enrich the set of clusters searched"),

						new FlaggedOption(
								"duploss weight",
								JSAP.STRING_PARSER,
								null,
								JSAP.NOT_REQUIRED,
								'l',
								"duploss",
								"Solves MGDL problem. Minimizes the number duplication and losses required"
										+ " to explain gene trees using DynaDup algorithm. Note that with this option, "
										+ "DynaDyp would be used *instead of* ASTRAL. "
										+ "Use -l 0 for standard (homomorphic) definition, and -l 1 for our new bd definition. "
										+ "Any value in between weights the impact of missing taxa somewhere between these two extremes. "
										+ "-l auto will automatically pick this weight. "), });
	}

	static Options readOptions(int criterion, boolean rooted,
			boolean extrarooted, double wh, JSAPResult config,
			List<Tree> mainTrees, List<List<String>> bootstrapInputSets)
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

		if (config.getBoolean("gene-only")
				&& config.getFile("bootstraps") != null) {
			exitWithErr("--gene-only and -b cannot be used together");
		}

		if (outfile == null) {
			if (config.getInt("branch annotation level") == 16) {
				File extraTreeFile = config.getFile("score species trees");
				freqPath = extraTreeFile.getAbsoluteFile().getParentFile()
						.getAbsolutePath();
			}
		} else {
			if (config.getInt("branch annotation level") == 16) {
				freqPath = outfile.getAbsoluteFile().getParentFile()
						.getAbsolutePath();
			}
			outfileName = config.getFile("output file") == null ? null : config
					.getFile("output file").getCanonicalPath();
		}

		if (config.getBoolean("duplication")
				&& config.contains("duploss weight")) {
			exitWithErr("dup and duploss options cannot be used together. Choose only one. ");
		}
		// johng23
		if (!config.getBoolean("cpu only")) {
			final int platformIndex = 0;
			final long deviceType = CL_DEVICE_TYPE_ALL;
			final int deviceIndex = 0;

			// Enable exceptions and subsequently omit error checks in this
			// sample
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
			Threading.contextProperties = new cl_context_properties();
			Threading.contextProperties.addProperty(CL_CONTEXT_PLATFORM, platform);

			// Obtain the number of devices for the platform
			int numDevicesArray[] = new int[1];
			clGetDeviceIDs(platform, deviceType, 0, null, numDevicesArray);
			int numDevices = numDevicesArray[0];

			// Obtain a device ID
			cl_device_id devices[] = new cl_device_id[numDevices];
			clGetDeviceIDs(platform, deviceType, numDevices, devices, null);
			for (int i = 0; i < numDevices; i++) {
				String deviceName = getString(devices[i], CL_DEVICE_NAME);
				String deviceVendor = getString(devices[i], CL_DEVICE_VENDOR);
				System.out.println("Device " + (i + 1) + " of " + numDevices
						+ ": " + deviceName + " " + devices[i] + " Vendor: " + deviceVendor);
			}
			System.out
					.println("Please enter the devices you'd like this program to use separated by spaces: ");
			Scanner in = new Scanner(System.in);
			ArrayList<cl_device_id> usedDevicesAL = new ArrayList<cl_device_id>();
			ArrayList<String> deviceVendorsAL = new ArrayList<String>();
			// while(in.hasNext()) {
			// usedDevicesAL.add(devices[in.nextInt()-1]);
			// }
			// testing only
			deviceVendorsAL.add(getString(devices[0], CL_DEVICE_VENDOR));
			usedDevicesAL.add(devices[0]);
			// usedDevicesAL.add(devices[1]);
			// usedDevicesAL.add(devices[2]);
			// usedDevicesAL.add(devices[3]);
			Threading.usedDevices = new cl_device_id[usedDevicesAL.size()];
			Threading.usedDevices = usedDevicesAL.toArray(Threading.usedDevices);
			Threading.deviceVendors = new String[deviceVendorsAL.size()];
			Threading.deviceVendors = deviceVendorsAL.toArray(Threading.deviceVendors);
			// cl_device_id device = devices[deviceIndex];
			// context = clCreateContext(contextProperties, 1, new
			// cl_device_id[]{device}, null, null, null);
			// System.out.println(usedDevices.length + " " +
			// usedDevices[0].toString());
			Threading.context = clCreateContext(Threading.contextProperties, Threading.usedDevices.length,
					Threading.usedDevices, null, null, null);
			// johng23 end
		}
		int numThreads = config.getInt("cpu threads");
		if (numThreads == -1) {
			numThreads = Runtime.getRuntime().availableProcessors();
		}
		
		Threading.startThreading(numThreads);
		

		if (Logging.timerOn) {
			System.err.println("Timer starts here");
			Logging.timer = System.currentTimeMillis();
		}

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
						alleles = s.split(" ", 3);
						species = alleles[0];
						alleles = alleles[2].split(" ");
					}
					for (String allele : alleles) {
						allele = allele.trim();
						if (taxonMap.containsKey(allele)) {
							System.err
									.println("The name mapping file is not in the correct format");
							System.err
									.println("A gene name can map to one only species name; check: "
											+ allele
											+ " which seems to appear at least twice: "
											+ taxonMap.get(allele)
											+ " & "
											+ species);
							System.exit(-1);
						} else if (alleles.length > 1 && allele.equals(species)) {
							System.err
									.println("Error: The species name cannot be identical to gene names when"
											+ "multiple alleles exist for the same gene: "
											+ allele);
							System.exit(-1);
						}
						// System.err.println("Mapping '"+allele+"' to
						// '"+species+"'");
						taxonMap.put(allele, species);
					}
				}

			} catch (Exception e) {
				br.close();
				throw new RuntimeException(
						"\n** Error **: Your name mapping file looks incorrect.\n   Carefully check its format. ",
						e);
			}
			br.close();
		}

		minleaves = config.contains("minleaves") ? config.getInt("minleaves")
				: null;
		samplingrounds = config.contains("samplingrounds") ? config
				.getInt("samplingrounds") : null;
		polylimit = config.contains("polylimit") ? config.getInt("polylimit")
				: null;

		try {

			// GlobalMaps.taxonIdentifier.taxonId("0");

			// System.err.println("Main input file: "+config.getFile("input file"));
			readInputTrees(mainTrees,
					readTreeFileAsString(config.getFile("input file")), rooted,
					true, false, minleaves,
					config.getInt("branch annotation level"), null);
			System.err.println(mainTrees.size() + " trees read from "
					+ config.getFile("input file"));

			GlobalMaps.taxonIdentifier.lock();

			Logging.logTimeMessage("");

		} catch (IOException e) {
			System.err.println("Error when reading trees.");
			System.err.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}

		if (mainTrees == null || mainTrees.size() == 0) {
			System.err.println("Empty list of trees. The function exits.");
			System.exit(1);
		} else {
			k = mainTrees.size();
		}

		if (taxonMap != null) {
			GlobalMaps.taxonNameMap = new TaxonNameMap(taxonMap);
		} else if (replace != null) {
			GlobalMaps.taxonNameMap = new TaxonNameMap(pattern, replace);
		} else {
			GlobalMaps.taxonNameMap = new TaxonNameMap();
		}

		if (config.getStringArray("keep") != null
				&& config.getStringArray("keep").length != 0) {
			if (outfileName == null) {
				throw new JSAPException(
						"When -k option is used, -o is also needed.");
			}
			for (String koption : config.getStringArray("keep")) {
				if ("completed".equals(koption) || "bootstraps".equals(koption)
						|| "bootstraps_norun".equals(koption)
						|| "searchspace_norun".equals(koption)
						|| "searchspace".equals(koption)) {
					keepOptions.add(koption);
				} else {
					throw new JSAPException("-k " + koption
							+ " not recognized.");
				}
			}
		}

		try {
			if (config.getFile("bootstraps") != null) {
				String line;
				BufferedReader rebuff = new BufferedReader(new FileReader(
						config.getFile("bootstraps")));
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

		if (config.getFile("bootstraps") != null
				|| config.getBoolean("gene-only")) {
			System.err.println("Bootstrapping with seed "
					+ config.getLong("seed"));
			for (int i = 0; i < config.getInt("replicates"); i++) {
				List<String> input = new ArrayList<String>();
				bootstrapInputSets.add(input);
				try {
					if (config.getBoolean("gene-sampling")) {
						for (int j = 0; j < k; j++) {
							input.add(bstrees.get(GlobalMaps.random.nextInt(k))
									.remove(0));
						}
					} else if (config.getBoolean("gene-only")) {
						for (int j = 0; j < k; j++) {
							input.add(mainTrees.get(
									GlobalMaps.random.nextInt(k)).toString());
						}
					} else {
						for (List<String> gene : bstrees) {
							input.add(gene.get(i));
						}
					}
				} catch (IndexOutOfBoundsException e) {
					exitWithErr("Error: You seem to have asked for "
							+ config.getInt("replicates")
							+ " but only "
							+ i
							+ " replicates could be created.\n"
							+ " Note that for gene resampling, you need more input bootstrap"
							+ " replicates than the number of species tree replicates.");
				}
				if (keepOptions.contains("bootstraps_norun")
						|| keepOptions.contains("bootstraps")) {
					String bsfn = outfile + ("." + i + ".bs");
					BufferedWriter bsoutbuffer = new BufferedWriter(
							new FileWriter(bsfn));
					for (String tree : input) {
						bsoutbuffer.write(tree + " \n");
					}
					bsoutbuffer.close();
				}
			}
			if (keepOptions.contains("bootstraps_norun")
					|| keepOptions.contains("bootstraps")) {
				System.err.println("bootstrap files written to files "
						+ outfile + ("." + 0 + ".bs") + " to " + outfile
						+ ("." + config.getInt("replicates") + ".bs"));
			}
			if (keepOptions.contains("bootstraps_norun")) {
				System.err
						.println("Exiting after outputting the bootstrap files");
				System.exit(0);
			}
		}

		Options options = new Options(rooted, extrarooted,
				config.getBoolean("exact"), criterion > 0, 1,
				config.getInt("extraLevel"), keepOptions.contains("completed"),
				keepOptions.contains("searchspace_norun")
						|| keepOptions.contains("searchspace"),
				!keepOptions.contains("searchspace_norun"),
				config.getInt("branch annotation level"),
				config.getDouble("lambda"), outfileName,
				samplingrounds == null ? -1 : samplingrounds,
				polylimit == null ? -1 : polylimit, freqPath, minleaves,
				config.getInt("gene repetition"));
		options.setDLbdWeigth(wh);
		options.setCS(1d);
		options.setCD(1d);

		return options;
	}

	public static void main(String[] args) throws Exception {

		long startTime = System.currentTimeMillis();

		JSAPResult config;
		int criterion = 2; // 2 for ASTRAL, 0 for dup, 1 for duploss
		boolean rooted = false;
		boolean extrarooted = false;
		double wh = 1.0D;

		List<Tree> mainTrees = new ArrayList<Tree>();
		List<List<String>> bootstrapInputSets = new ArrayList<List<String>>();
		BufferedWriter outbuffer;

		System.err
				.println("\n================== ASTRAL ===================== \n");
		System.err.println("This is ASTRAL version " + _version);

		jsap = getJSAP();
		config = jsap.parse(args);
		if (jsap.messagePrinted()) {
			exitWithErr("");
		}

		if (config.getBoolean("duplication")) {
			criterion = 0;
			rooted = true;
			extrarooted = true;
			System.err
					.println("Using DynaDup application, minimizing MGD (not ASTRAL).");
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
					exitWithErr("duploss weight has to be between 0 and 1");
				}
				;
			}
			System.err
					.println("Using DynaDup application, minimizing MGDL (not ASTRAL).");
		}

		System.err.println("Gene trees are treated as "
				+ (rooted ? "rooted" : "unrooted"));

		GlobalMaps.random = new Random(config.getLong("seed"));

		Options options = readOptions(criterion, rooted, extrarooted, wh,
				config, mainTrees, bootstrapInputSets);

		File outfile = config.getFile("output file");
		if (outfile == null) {
			outbuffer = new BufferedWriter(new OutputStreamWriter(System.out));

		} else {

			outbuffer = new BufferedWriter(new FileWriter(outfile));
		}

		String outgroup = GlobalMaps.taxonNameMap.getSpeciesIdMapper()
				.getSpeciesName(0);

		List<String> toScore = null;
		if (config.getFile("score species trees") != null) {
			System.err.println("Scoring "
					+ config.getFile("score species trees"));
			toScore = readTreeFileAsString(config
					.getFile("score species trees"));
			runScore(criterion, rooted, mainTrees, outbuffer, options,
					outgroup, toScore);
		} else {

			runInference(config, criterion, rooted, extrarooted, mainTrees,
					outbuffer, bootstrapInputSets, options, outgroup);
		}
		
		Threading.shutdown();
		
		
		System.err.println("ASTRAL finished in "
				+ (System.currentTimeMillis() - startTime) / 1000.0D + " secs");
	}

	private static void runScore(int criterion, boolean rooted,
			List<Tree> mainTrees, BufferedWriter outbuffer, Options options,
			String outgroup, List<String> toScore)
			throws FileNotFoundException, IOException {
		System.err.println("Scoring: " + toScore.size() + " trees");

		AbstractInference inference = initializeInference(criterion, mainTrees,
				new ArrayList<Tree>(), options);
		double score = Double.NEGATIVE_INFINITY;
		List<Tree> bestTree = new ArrayList<Tree>();
		for (String trs : toScore) {
			List<Tree> trees = new ArrayList<Tree>();
			readInputTrees(trees, Arrays.asList(new String[] { trs }), rooted,
					true, true, null, 1, false ? // config.getBoolean("scoreall")?
					outgroup
							: null);
			Tree tr = trees.get(0);

			double nscore = inference.scoreSpeciesTreeWithGTLabels(tr, true);

			if (nscore > score) {
				score = nscore;
				bestTree.clear();
				bestTree.add(tr);
			} else if (nscore == score) {
				bestTree.add(tr);
			}

			if (!GlobalMaps.taxonNameMap.getSpeciesIdMapper()
					.isSingleIndividual()) {
				System.err.println("Scored tree with gene names:\n"
						+ tr.toNewickWD());
			}

			GlobalMaps.taxonNameMap.getSpeciesIdMapper().gtToSt(
					(MutableTree) tr);

			if (options.getBranchannotation() != 12) {
				writeTreeToFile(outbuffer, tr);
			}
		}
		if (options.getBranchannotation() == 12) {
			for (Tree bt : bestTree)
				writeTreeToFile(outbuffer, bt);
		}

		outbuffer.close();
	}

	private static void runInference(JSAPResult config, int criterion,
			boolean rooted, boolean extrarooted, List<Tree> mainTrees,
			BufferedWriter outbuffer, List<List<String>> bootstrapInputSets,
			Options options, String outgroup) throws JSAPException,
			IOException, FileNotFoundException {

		System.err.println("All output trees will be *arbitrarily* rooted at "
				+ outgroup);
		List<Tree> extraTrees = new ArrayList<Tree>();

		try {

			if (config.getFile("extra trees") != null) {
				readInputTrees(extraTrees,
						readTreeFileAsString(config.getFile("extra trees")),
						extrarooted, true, false, null, 1, null);
				System.err.println(extraTrees.size()
						+ " extra trees read from "
						+ config.getFile("extra trees"));
			}

			if (config.getFile("extra species trees") != null) {
				readInputTrees(extraTrees,
						readTreeFileAsString(config
								.getFile("extra species trees")), extrarooted,
						true, true, null, 1, null);
				System.err.println(extraTrees.size()
						+ " extra trees read from "
						+ config.getFile("extra trees"));
			}

		} catch (IOException e) {
			System.err.println("Error when reading extra trees.");
			System.err.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}

		int j = 0;
		List<Tree> bootstraps = new ArrayList<Tree>();
		for (List<String> input : bootstrapInputSets) {
			System.err.println("\n======== Running bootstrap replicate " + j++);
			List<Tree> trees = new ArrayList<Tree>();
			readInputTrees(trees, input, rooted, false, false,
					options.getMinLeaves(),
					config.getInt("branch annotation level"), null);
			bootstraps.add(runOnOneInput(criterion, extraTrees, outbuffer,
					trees, null, outgroup, options));
		}

		if (bootstraps != null && bootstraps.size() != 0) {
			STITree<Double> cons = (STITree<Double>) Utils
					.greedyConsensus(bootstraps, false, GlobalMaps.taxonNameMap
							.getSpeciesIdMapper().getSTTaxonIdentifier(), false);
			cons.rerootTreeAtNode(cons.getNode(outgroup));
			Trees.removeBinaryNodes(cons);
			Utils.computeEdgeSupports(cons, bootstraps);
			writeTreeToFile(outbuffer, cons);
		}

		Logging.logTimeMessage(" ");

		System.err.println("\n======== Running the main analysis");
		runOnOneInput(criterion, extraTrees, outbuffer, mainTrees, bootstraps,
				outgroup, options);

		outbuffer.close();
	}

	private static Tree runOnOneInput(int criterion, List<Tree> extraTrees,
			BufferedWriter outbuffer, List<Tree> input,
			Iterable<Tree> bootstraps, String outgroup, Options options) {
		long startTime;
		startTime = System.currentTimeMillis();

		AbstractInference inferenceConsumer = initializeInference(criterion, input, extraTrees, options);
		inferenceConsumer.setup();
		

		List<Solution> solutions = inferenceConsumer.inferSpeciesTree();
		Logging.logTimeMessage(" CommandLine 667: ");

		System.err.println("Optimal tree inferred in "
						+ (System.currentTimeMillis() - startTime) / 1000.0D
						+ " secs.");
		System.err.println("Weight calculation cumulatively took "
				+ Polytree.time / 1000000000.0D + " secs");

		Tree st = processSolution(outbuffer, bootstraps, outgroup, inferenceConsumer, solutions);

		return st;
	}


	private static boolean isGeneResamplign(JSAPResult config) {
		return config.getBoolean("gene-sampling")
				|| config.getBoolean("gene-only");
	}

	private static Tree processSolution(BufferedWriter outbuffer,
			Iterable<Tree> bootstraps, String outgroup,
			AbstractInference inference, List<Solution> solutions) {
		Logging
				.logTimeMessage(" CommandLine 684: ");

		Tree st = solutions.get(0)._st;
		Logging
				.logTimeMessage(" CommandLine 690: ");

		System.err.println(st.toNewick());

		st.rerootTreeAtNode(st.getNode(outgroup));

		Trees.removeBinaryNodes((MutableTree) st);

		// TODO: MULTIND.
		GlobalMaps.taxonNameMap.getSpeciesIdMapper().stToGt((MutableTree) st);

		inference.scoreSpeciesTreeWithGTLabels(st, false);

		GlobalMaps.taxonNameMap.getSpeciesIdMapper().gtToSt((MutableTree) st);

		if ((bootstraps != null) && (bootstraps.iterator().hasNext())) {
			for (Solution solution : solutions) {
				Utils.computeEdgeSupports((STITree<Double>) solution._st,
						bootstraps);
			}
		}
		writeTreeToFile(outbuffer, solutions.get(0)._st);

		return st;
	}

	private static AbstractInference initializeInference(int criterion,
			List<Tree> trees, List<Tree> extraTrees, Options options) {
		AbstractInference inference;
		if (criterion == 1 || criterion == 0) {
			inference = new DLInference(options, trees, extraTrees);
		} else if (criterion == 2) {
			inference = new WQInference(options, trees, extraTrees);
		} else {
			throw new RuntimeException("criterion not set?");
		}
		return inference;
	}

	private static List<String> readTreeFileAsString(File file)
			throws FileNotFoundException, IOException {
		String line;
		List<String> trees = new ArrayList<String>();
		BufferedReader treeBufferReader = new BufferedReader(new FileReader(
				file));
		while ((line = treeBufferReader.readLine()) != null) {
			if (line.length() > 0) {
				line = line.replaceAll("\\)[^,);]*", ")");
				trees.add(line);
			}
		}
		treeBufferReader.close();
		return trees;

	}

	private static void readInputTrees(List<Tree> trees, List<String> lines,
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
				if (line.length() == 0) {
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
					// System.err.println(".");
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
							// if (!stLablel) {
							GlobalMaps.taxonIdentifier.taxonId(leaves[i]);
							// } else {
							// GlobalMaps.taxonNameMap.getSpeciesIdMapper().speciesId(leaves[i]);
							// }
						}
					} else {
						allleaves.addAll(Arrays.asList(leaves));
					}
				}
				if (annotation == 6) {
					for (String leaf : allleaves) {
						GlobalMaps.taxonIdentifier.taxonId(leaf);
					}
				}
			}
		} catch (ParseException e) {
			throw new RuntimeException("Failed to Parse Tree number: " + l, e);
		}
		if (skipped.size() > 0) {
			System.err
					.println("Skipping the following tree(s) because they had less than "
							+ minleaves + " leaves: \n" + skipped);
		}
	}

	private static void writeTreeToFile(BufferedWriter outbuffer, Tree t) {
		try {
			outbuffer.write(t.toStringWD() + " \n");
			outbuffer.flush();
		} catch (IOException e) {
			System.err.println("Error when writing the species tree");
			System.err.println(e.getMessage());
			e.printStackTrace();
		}
	}

	private static String getString(cl_device_id device, int paramName) {
		long size[] = new long[1];
		clGetDeviceInfo(device, paramName, 0, null, size);
		byte buffer[] = new byte[(int) size[0]];
		clGetDeviceInfo(device, paramName, buffer.length, Pointer.to(buffer),
				null);
		return new String(buffer, 0, buffer.length - 1);
	}
}
