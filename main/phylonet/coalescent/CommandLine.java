package phylonet.coalescent;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import phylonet.coalescent.GlobalMaps.TaxonNameMap;
import phylonet.tree.io.NewickReader;
import phylonet.tree.io.ParseException;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITree;

public class CommandLine {
	
	protected static boolean _print = true;
	
	protected static String _versinon = "4.2.1";


	public static void main(String[] args) {
		if ((args == null) || args.length == 0 || (args[0].equals("-h"))
				|| (args.length < 1)) {
			printUsage();
			return;
		}
		
		int criterion = 2;
		boolean rooted = false;
		//boolean fast = false;
		boolean extrarooted = false;
		boolean exactSolution = false;
		
		Map<String, String> taxonMap = null;
		String rep = null;
		String pattern = null;
		List<Tree> trees = null;
		List<Tree> extraTrees = null;
		String output = null;
		STITree scorest = null;
		// boolean explore = false;
		// double proportion = 0.0D;
		// boolean exhaust = false;
		double bootstrap = 1.0D;
		double cs = 1.0D;
		double cd = 1.0D;
		double time = -1.0D;
		double wd = 1.0D;
		double wh = 1.0D;
		boolean unresolved = false;		
		long startTime = System.currentTimeMillis();
		int alg = -1;
		String line;
		BufferedReader treeBufferReader = null;
		BufferedReader extraTreebuffer = null;
		try {
			List<String[]> options = getOptions(args);
			for (String[] option : options) {
				if (option[0].equals("-i")) {
					if (option.length != 2) {
						printUsage();
						return;
					}
					treeBufferReader = new BufferedReader(new FileReader(
							option[1]));

					trees = new ArrayList();
					// String line;
				} else if (option[0].equals("-ex")) {
					if (option.length != 2) {
						printUsage();
						return;
					}					
					extraTreebuffer = new BufferedReader(new FileReader(
							option[1]));

					extraTrees = new ArrayList();					
				} /*else if (option[0].equals("-st")) {
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
				}*//* else if (option[0].equals("-a")) {
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
				}*/ else if (option[0].equals("-o")) {
					if (option.length != 2) {
						printUsage();
						return;
					}
					output = option[1];
					setPrinting (false);
				} else if (option[0].equals("-x")) {
					if (option.length != 1) {
						printUsage();
						return;
					}
					// exhaust = true;
				} else if (option[0].equals("-e")) {
					if (option.length > 2) {
						printUsage();
						return;
					}
					// explore = true;
					if (option.length != 2)
						continue;
					try {
						// proportion = Double.parseDouble(option[1]);
					} catch (NumberFormatException e) {
						System.err.println("Error in reading parameter");
						printUsage();
						return;
					}

				} else if (option[0].equals("-b")) {
					if (option.length > 2) {
						printUsage();
						return;
					}
					if (option.length != 2)
						continue;
					try {
						bootstrap = Double.parseDouble(option[1]);
						if ((bootstrap <= 1.0D) && (bootstrap > 0.0D))
							continue;
						printUsage();
					} catch (NumberFormatException e) {
						System.err.println("Error in reading parameter");
						printUsage();
						return;
					}

				} else if (option[0].equals("-fa")) {
					alg = Integer.parseInt(option[1]);

				}else if (option[0].equals("-cs")) {
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

				} else if (option[0].equals("-ur")) {
					if ((option.length != 1) || (time != -1.0D)) {
						printUsage();
						return;
					}
					unresolved = true;
//				} else if (option[0].equals("-f")) {
//					if ((option.length != 1)) {
//						printUsage();
//						return;
//					}
//					fast = true;
				} else if (option[0].equals("-u")) {
					if ((option.length != 1) || (time != -1.0D)) {
						printUsage();
						return;
					}
					rooted = false;
				} else if (option[0].equals("-xu")) {
					if ((option.length != 1) || (time != -1.0D)) {
						printUsage();
						return;
					}
					extrarooted = false;
				} else if (option[0].equals("-t")) {
					if ((option.length != 2) || (unresolved)) {
						printUsage();
						return;
					}
					if (option.length != 2)
						continue;
					try {
						time = Double.parseDouble(option[1]);
						if (time > 0.0D)
							continue;
						printUsage();
					} catch (NumberFormatException e) {
						System.err.println("Error in reading parameter");
						printUsage();
						return;
					}
				} else if (option[0].equals("-xt")) {
					if (option.length != 1) {
						printUsage();
						return;
					}
					exactSolution = true;
				} else if (option[0].equals("-dl")) {
					if (criterion != 2) {
						System.err.println("You should choose only of the following options: -d (duplications) -dl (duploss)");
						printUsage();
						return;
					}
					criterion = 1;
					rooted = true;
					extrarooted = true;
					if (option.length != 2) {
						printUsage();
						return;
					}					
					try {
						if (option[1].equals("auto")) {
							wh = -1;
							continue;
						} else {
							wh = Double.parseDouble(option[1]);
							if (wh >= 0.0D)
								continue;
						}
						printUsage();
						return;
					} catch (NumberFormatException e) {
						System.err.println("Error in reading parameter wd");
						printUsage();
						return;
					}
				} else if (option[0].equals("-d")) {
					if (criterion != 2) {
						System.err.println("You should choose only of the following options: -d (duplications) -dl (duploss))");
						printUsage();
						return;
					}
					criterion = 0;
					rooted = true;
					extrarooted = true;
					if (option.length != 1) {
						printUsage();
						return;
					}
				}  /*else if (option[0].equals("-wq")) {
					if (criterion != -1) {
						System.err.println("You should choose only of the following options: -d (duplications) -dl (duploss) -wq (weighted qurtets)");
						printUsage();
						return;
					}
					System.err.println("Criterion set to weighted quartets. Gene trees will be treated as unrooted.");
					criterion = 2;
					rooted = false;
					extrarooted = false;
					if (option.length != 1) {
						printUsage();
						return;
					}
				}*/ else {
					printUsage();
					return;
				}
			}

			if (treeBufferReader == null) {
				System.err.println("The input file has not been specified.");
				printUsage();
				return;
			}

			System.err.println("Gene trees are treated as "
					+ (rooted ? "rooted" : "unrooted"));
			int l = 0;
			
			try {
				while ((line = treeBufferReader.readLine()) != null) {
					l++;
					Set<String> previousTreeTaxa = new HashSet<String>();
					if (line.length() > 0) {
						line = line.replaceAll("\\)[^,);]*", ")");
						NewickReader nr = new NewickReader(new StringReader(line));
						if (rooted) {
							STITree gt = new STITree(true);
							nr.readTree(gt);
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
							trees.add(gt);
						} else {						
							Tree tr = nr.readTree();
							trees.add(tr);
						}
					}
				}
				treeBufferReader.close();
			} catch (ParseException e) {
				treeBufferReader.close();
				throw new RuntimeException("Failed to Parse Tree number: " + l ,e);
			}			

			if (extraTreebuffer != null) {
				while ((line = extraTreebuffer.readLine()) != null) {
					if (line.length() > 0) {	
						line = line.replaceAll("\\)[^,);]*", ")");
						NewickReader nr = new NewickReader(
								new StringReader(line));
						if (extrarooted) {
							STITree gt = new STITree(true);
							nr.readTree(gt);
							extraTrees.add(gt);
						} else {
							Tree tr = nr.readTree();
							extraTrees.add(tr);
						}
					}
				}				
				extraTreebuffer.close();
			}
			
		} catch (IOException e) {
			System.err.println("Error when reading trees. The function exits.");
			System.err.println(e.getMessage());
			e.printStackTrace();
			return;
		} catch (ParseException e) {
			System.err
					.println("Error when parsing the Newick representation from input file.");
			System.err.println(e.getMessage());
			e.printStackTrace();
			return;
		}

		if (trees.size() == 0) {
			System.err.println("Empty list of trees. The function exits.");
			return;
		}

		if (_print) {
			System.err.println("Reading trees in "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		}

		startTime = System.currentTimeMillis();
				
		if (taxonMap != null) {
			GlobalMaps.taxonNameMap = new TaxonNameMap(taxonMap);
		} else if (rep != null) {	
			GlobalMaps.taxonNameMap = new TaxonNameMap (pattern, rep);
		}

		
		Inference inference;
		
		if (criterion == 1 || criterion == 0) {
			inference = new DLInference(rooted, extrarooted, 
					trees, extraTrees, exactSolution,criterion > 0);			
		} else if (criterion == 2) {
			inference = new WQInference(rooted, extrarooted, 
					trees, extraTrees, exactSolution,criterion > 0, alg);
		} else {
			throw new RuntimeException("You should choose one of the following options: -d (duplications) -dl (duploss) -wq (weighted qurtets)");
		}
		
		inference.setDLbdWeigth(wh); 
		inference.setCS(cs);
		inference.setCD(cd);

		if (scorest != null) {
			inference.scoreGeneTree(scorest);
			System.exit(0);
		}
		
		List<Solution> solutions = inference.inferSpeciesTree();

		//if (_print) {
			System.err.println("Optimal tree inferred in "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		//}

		if ((_print)) {
			String metric;
			if (criterion == 0) {
				metric = "duplications";
			} else if (criterion == 1) {
				metric = "duplication+loss (homomorphic)";
			} else {
				metric = "weighted quartet sum";
			}
			for (Solution s : solutions)
				System.out.println(
					        s._st.toStringWD()
						+ " \n"
						+ s._totalCoals
						+ " " + metric + " in total");
		} else
			try {
				FileWriter fw = new FileWriter(output);
				for (Solution s : solutions) {
					fw.write(s._st.toString()+ " \n");
				}
				fw.close();
			} catch (IOException e) {
				System.err.println("Error when writing the species tree");
				System.err.println(e.getMessage());
				e.printStackTrace();
			}
	}


	protected static void printUsage() {
		System.out
				.println("ASTRAL (version "+_versinon + " ): species tree inference from unrooted gene trees.\n"
						+ "The ASTRAL algorithm maximizes the number of shared quartet trees with the collection of all gene trees. \n"
						+ "The result of this optimization problem is statistically consistent under the multi-species coalesent model. \n"
						+ "Use -d to instead minimize duplications or -dl to minimize dupication and loss (using the DynaDup algorithm).");
		System.out.println("\nUsage is:");
		System.out
				.println("\t-i input [-o output] [-ex extra_trees] [-xt] [-d | -dl N] ");// [-s species tree] [-cs number] [-cd number] [-a mapping] [-wd duplication weight]");
		System.out
				.println("\t-i input file: A file containing gene trees. (required)");
/*		System.out
				.println("\t-st species tree file: The file containing a species tree to be scored.\n" +
						 "\t                       If this option is provided the software only scores the species tree.");
		System.out
				.println("\t-a mapping file: The file containing the mapping from alleles to speceis if multiple alleles sampled.\n" +
						 "\t                 Alternatively, two reqular expressions for automatic name conversion (optional)");
*/		System.out
				.println("\t-o output file: A filename for storing the output species tree. Defaults to stdout. (optional)");
		System.out.println("\t-xt: find the exact solution by looking at all clusters (recommended only for small (<18) numer of taxa.");
		//System.out.println("\t-u treat input gene trees as unrooted (Not implemented!)");
		System.out.println("\t-ex: provide extra trees used to enrich the set of clusters searched");
		//System.out.println("\t-xu treat extra trees input gene trees as unrooted (Not implemented!)");
		System.out.println("\t-d: optimizes duplications (DynaDup). This option invokes DynaDyp instead of ASTRAL.");
		System.out.println("\t-dl N: optimize duplications and losses (DynaDup). This option invokes DynaDyp instead of ASTRAL.\n"
				+ 		   "\t    Use -dl 0 for standard (homomorphic) definition, and -dl 1 for ``bd'' definition.\n"+
						   "\t    Any value in between weights the impact of missing taxa on the tree.");
		//System.out.println("\t-wq optimizes weighted quartet score (useful for ILS).");
		System.out.println("\t-cs and -cd: these two EXPERIMENTAL options set two parameters (cs and cd) to a value between 0 and 1. \n" +
						   "\t    For any cluster C if |C| >= cs*|taxa|, we add complementary clusters (with respect to C) of all subclusters of C\n" +
						   "\t    if size of the subcluster is >= cd*|C|.\n" +
						   "\t    By default cs = cd = 1; so no extra clusters are added. Lower cs and cd values could result in better scores\n" +
						   "\t    (especially when gene trees have missing data) but can also increase the running time quite substantially.");
		
		//System.out.println("\t-f perform fast and less-accurate subtree-bipartition based search (Not implemented!).");
		System.out.println();
	}


	protected static List<String[]> getOptions(String[] args) {
		LinkedList<String []> opts = new LinkedList<String[]>();
		LinkedList<String> argList = new LinkedList<String>();
	
		int i = 0;
		while (i < args.length) {
			if (args[i].charAt(0) != '-') {
				printUsage();
				System.exit(-1);
			}
	
			argList.clear();
			argList.addFirst(args[i]);
			i++;
	
			while ((i < args.length) && (args[i].charAt(0) != '-')) {
				argList.addLast(args[i]);
				i++;
			}
			String[] argArray = new String[argList.size()];
			argList.toArray(argArray);
	
			opts.addLast(argArray);
		}
		return opts;
	}


	public static void setPrinting(boolean print) {
		CommandLine._print = print;
	}
}
