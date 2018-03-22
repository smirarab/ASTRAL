package phylonet.coalescent;

/*
######################
######################
PhyDstar program
written in JAVA 1.5
by Alexis Criscuolo

criscuol@lirmm.fr
criscuol@dpt-info.u-strasbg.fr
criscuol@pasteur.fr

MPG group, ISEM
Universite de Montpellier II - CC 065
34095 MONTPELLIER Cedex 05
FRANCE

MAB team, LIRMM
161, rue Ada
34000 MONTPELLIER
FRANCE

BT team, LSIIT
Pole API, Boulevard Sebastien Brant - BP 10413
67412 ILLKIRCH Cedex
FRANCE

Phylogenomics team
Institut Pasteur, 28 rue du Dr Roux
75724 PARIS cedex 15
FRANCE

Reference to cite:
     Fast NJ-like algorithms to deal with incomplete distance matrices
     Criscuolo, Alexis; Gascuel, Olivier
     BMC Bioinformatics 2008, 9:166
######################
######################
*/

import java.io.*;
import java.util.*;

public class PhyDstar
{
    //##### menu #####
     String[] option;
     BufferedReader cin;
     String ligne;
     int multifile;                     // number of dataset
     int p;                             // threshold parameter for agglomeration
     String choix;                      // for the start menu
     String outgroup;                   // name of the outgroup taxon
     boolean ok, nj, unj, bionj, mvr;   // to make some tests

    //##### data structures #####
     int n;                // number of taxa
     String[] taxa;        // list of taxa name
     Matrix dM, dMinit;    // distance matrix
     int[][] tT;           // tree topology
     double[][] tD;        // tree branch lengths
     int center;           // the coordinate of the center inside the tree structure
     boolean[] activeTaxa; // list of taxa present in the distance matrix
     int[] tree2dist;      // indicate the rank of a node in the distance matrix
     int hole;             // number of missing distances
     int[] missing;        // number of missing distances per taxa
     Matrix vM;            // variance matrix
    
    //##### reading data #####
     String dinfile, vinfile;
     BufferedReader ind;   // to read the distance input file 
     BufferedReader inv;   // to read the variance input file 
     int endOfLine;        // for the lower triangular format
     int b;                // temp var
     int b2;               // temp var
     String tempString;    // temp String var
     double tempDouble;    // temp double var
     double tempDouble2;   // temp double var
     int tempInt;          // temp int var
     String line;
    
    //##### building tree #####
     int r;                // size of the distance matrix during the agglomerative algorithm
     double[] sumR;        // sum of the known distances for Q criterion
     double[] sumS;        // sum of the weights
     double _Q;            // current value of the Q criterion
     double _Qmax;         // max value of the Q criterion
     Matrix _Qstar;        // values of the Q* criterion
     Matrix _Sstar;        // normalizers of the Q* criterion
     double _Qcur;         // current value of the Q* criterion
     double _Scur;         // current normalizer of the Q* criterion
     Comparator<Object> compDesc;
     TreeSet<CriterionTuple> tabQstar;
     CriterionTuple ct;
     Iterator<CriterionTuple> it;
     double tabQstar_min;  // min value of the tabQ table
     double dNcur;         // current value of the discrete N* criterion
     double cNcur;         // current value of the continuous N* criterion
     double dNmax;         // max value of the discrete N* criterion
     double cNmax;         // max value of the continuous N* criterion
     double _Smax;         // normalizer of the max value of the discrete N* criterion
     int _j;               // ij is the current neighborable pair to test
     int _i;               // ij is the current neighborable pair to test
     int x;                // first taxon of the agglomerated pair
     int y;                // second taxon of the agglomerated pair
     double sum1;
     double sum2;
     int u;                // rank of the created node
     int rankmin;          // rank of the minimum rank active taxon
     double txu;           // length of branch T_xu
     double tyu;           // length of branch T_yu
     double un[];          // weights for computing reduction distance matrice (UNJ*)
     double mu;            // constant for the computing of w's (MVR*)
     double[] w;           // weights for computing branch lengths
     double[] lambda;      // weights for computing reduction distance matrices 

    //##### newick output tree #####
     Writer out;        // to write the trees inside an output file 
     String newick;             // outputed tree in newick format
     BitSet nodeFixed;          // to set the postfixed nodes
     int root;                  // the leaf root of the outputed tree
    
    public  String main(String[] args) throws IOException {

    compDesc = new Comparator<Object>() {     // a comparator dedicated to the DistanceTuple class
        CriterionTuple ct1; CriterionTuple ct2;
        public int compare (Object o1 , Object o2) {
        if ( (! ( o1 instanceof CriterionTuple )) || (! ( o2 instanceof CriterionTuple )) ) throw new ClassCastException();
        ct1 = (CriterionTuple) o1; ct2 = (CriterionTuple) o2;
        if ( ct1.Q < ct2.Q ) return 1;
        else {
            if ( ct1.Q > ct2.Q ) return -1;
            else { // ct1.Q == ct2.Q 
            if ( ct1.I < ct2.I ) return -1;
            else {
                if ( ct1.I > ct2.I ) return 1;
                else { // ct1.I == ct2.I
                if ( ct1.J < ct2.J ) return -1;
                else {
                    if ( ct1.J > ct2.J ) return 1;
                    else // ct1.J == ct2.J
                    return 0;
                }
                }
            }
            }
        }
        }
        public boolean equals (Object o1 , Object o2) {
        if ( (! ( o1 instanceof CriterionTuple )) || (! ( o2 instanceof CriterionTuple )) ) throw new ClassCastException();
        ct1 = (CriterionTuple) o1; ct2 = (CriterionTuple) o2;
        return ( (ct1.Q == ct2.Q) && (ct1.I == ct2.I) && (ct1.J == ct2.J) );
        }
    };
    


    option = new String[11];
    option[0] = "No";   // run the algorithms
    option[1] = "No";   // lower-traingular distance matrix?
    option[2] = "No";   // subreplicates matrix?
    option[3] = "No";   // multiple datasets?
    option[4] = "No";   // outgroup root?
    option[5] = "Yes";  // negative branch lengths? 
    option[6] = "Yes";  // binary tree?
    option[7] = "BioNJ*";  // method?
    option[8] = "No";   // input.v variance matrix?
    option[9] = "15";   // number of neighborable pairs
    option[10] = "No";  // confidence value at branches
    dinfile = ""; vinfile = "";
    multifile = 1;
    p = (new Integer(option[9])).intValue();
    outgroup = "";
    nj = false;
    unj = false;
    bionj = true;
    mvr = false;

    
    if ( args.length > 0 ) {  // command line options
        b = -1;
        while ( ++b < args.length ) {
        if ( args[b].toLowerCase().equals("-i") ) {
           
            dinfile = args[++b];
            ind = new BufferedReader(new StringReader(dinfile));
            
            continue;
        }
        if ( args[b].toLowerCase().equals("-v") ) {
            try { 
            vinfile = args[++b];
            inv = new BufferedReader(new FileReader(new File(vinfile)));
            } catch( FileNotFoundException e ) {
            System.err.println("Incorrect input file name (option -v)");
            System.exit(0);
            }
            option[8] = "Yes";
            continue;
        }
        if ( args[b].toLowerCase().equals("-d") ) {
            line = args[++b].toUpperCase();
            if ( line.equals("NJ") ) { nj = true; unj = false; bionj = false; mvr = false; option[7] = "NJ*"; }
            if ( line.equals("UNJ") ) { nj = false; unj = true; bionj = false; mvr = false; option[7] = "UNJ*"; }
            if ( line.equals("BIONJ") ) { nj = false; unj = false; bionj = true; mvr = false; option[7] = "BioNJ*"; }
            if ( line.equals("MVR") ) { nj = false; unj = false; bionj = false; mvr = true; option[7] = "MVR*"; }
            continue;
        }
        if ( args[b].toLowerCase().equals("-p") ) {
            line = args[++b];
            try { 
            p = Integer.parseInt(line); 
            if ( p < 1 ) p = 15;
            }
            catch ( NumberFormatException e ) {}
            if ( p > 1 ) option[9] = "" + p;
            continue;
        }
        if ( args[b].toLowerCase().equals("-n") ) {
            line = args[++b].toUpperCase();
            if ( line.equals("Y") ) option[5] = "Yes";
            if ( line.equals("N") ) option[5] = "No";
            continue;
        }
        if ( args[b].toLowerCase().equals("-b") ) {
            line = args[++b].toUpperCase();
            if ( line.equals("Y") ) option[6] = "Yes";
            if ( line.equals("N") ) option[6] = "No";
            continue;
        }
        if ( args[b].toLowerCase().equals("-o") ) {
            outgroup = args[++b];
            option[4] = "Yes (" + outgroup + ")";
            continue;
        }
        if ( args[b].toLowerCase().equals("-l") ) {
            option[1] = "Yes";
            continue;
        }
        if ( args[b].toLowerCase().equals("-s") ) {
            option[2] = "Yes";
            continue;
        }
        if ( args[b].toLowerCase().equals("-m") ) {
            line = args[++b];
            try { 
            multifile = Integer.parseInt(line); 
            if ( multifile < 1 ) multifile = 1;
            }
            catch ( NumberFormatException e ) {}
            if ( multifile > 1 ) option[3] = "Yes, " + multifile + " sets";
            continue;
        }
        if ( args[b].toLowerCase().equals("-c") ) {
            option[10] = "Yes";
            continue;
        }
        }
        if ( dinfile.equals("") ) {
        System.err.println("  No distance matrix file name");
        System.exit(0);
        }    
    }
    else {        // menu
        System.err.println("\nPhyD* version 1.1");
        cin = new BufferedReader((new InputStreamReader(System.in)));
        while ( option[0].equals("No") ) {
        System.err.println("\n\n\n\nSettings for this run:");
        System.err.println("   D             Method (BioNJ*, MVR*, NJ*, UNJ*)?  " + option[7]);
        if ( option[7].equals("MVR*") )
            System.err.println("   V                     Use variance matrix file?  " + option[8]);
        System.err.println("   P    Taxon pairs selected by NJ-like filtering?  " + option[9]);
        System.err.println("   N              Negative branch lengths allowed?  " + option[5]);
        System.err.println("   B                                  Binary tree?  " + option[6]);
        System.err.println("   O                                Outgroup root?  " + option[4]);
        System.err.println("   C                Confidence values at branches?  " + option[10]);
        System.err.println("   L                 Lower-triangular data matrix?  " + option[1]);
        System.err.println("   S                                Subreplicates?  " + option[2]);
        System.err.println("   M                   Analyse multiple data sets?  " + option[3]);
        System.err.print("\n   Y to accept these or type the letter for one to change   ");
        
        choix = cin.readLine();
        
        if ( (choix.equals("Y")) || (choix.equals("y")) )
            option[0] = "Yes";
        
        if ( (choix.equals("L")) || (choix.equals("l")) )
            if ( option[1].equals("No") ) option[1] = "Yes"; else option[1] = "No";
        
        if ( (choix.equals("S")) || (choix.equals("s")) )
            if ( option[2].equals("No") ) option[2] = "Yes"; else option[2] = "No";
        
        if ( (choix.equals("M")) || (choix.equals("m")) )
            if ( multifile == 1 ) {
            System.err.print("   How many data sets? ");
            choix = cin.readLine();
            multifile = (new Integer(choix)).intValue();
            option[3] = "Yes, " + multifile + " sets";
            }
            else {
            multifile = 1;
            option[3] = "No";
            }
        
        if ( (choix.equals("O")) || (choix.equals("o")) )
            if ( outgroup.equals("") ) {
            System.err.print("   Please enter the outgroup taxon name> ");
            choix = cin.readLine();
            outgroup = choix;
            option[4] = "Yes (" + outgroup + ")";
            }
            else {
            outgroup = "";
            option[4] = "No";
            }
        
        if ( (choix.equals("N")) || (choix.equals("n")) )
            if ( option[5].equals("No") ) option[5] = "Yes"; else option[5] = "No";
        
        if ( (choix.equals("B")) || (choix.equals("b")) )
            if ( option[6].equals("No") ) option[6] = "Yes"; else option[6] = "No";
        
        if ( (choix.equals("D")) || (choix.equals("d")) ) {
            if ( option[7].equals("NJ*") ) {
            option[7] = "UNJ*"; nj = false; unj = true;
            }
            else {
            if ( option[7].equals("UNJ*") ) {
                option[7] = "BioNJ*"; unj = false; bionj = true;
            }
            else {
                if ( option[7].equals("BioNJ*") ) {
                option[7] = "MVR*"; bionj = false; mvr = true;
                }
                else {
                option[7] = "NJ*"; mvr = false; nj = true;
                }
            }
            }
        }
        
        if ( ((choix.equals("V")) || (choix.equals("v"))) && option[7].equals("MVR*") ) 
            if ( option[8].equals("No") ) option[8] = "Yes"; else option[8] = "No";
        
        if ( (choix.equals("P")) || (choix.equals("p")) ) {
            System.err.print("   How many neighborable pairs to test? ");
            choix = cin.readLine();
            p = (new Integer(choix)).intValue();
            option[9] = "" + p;
        }

        if ( (choix.equals("C")) || (choix.equals("c")) ) 
            if ( option[10].equals("No") ) option[10] = "Yes"; else option[10] = "No";
        
        }
    
        choix = "input.d";
        ok = true;
        ind = null; 
        while (ok) {
        try {
            dinfile = choix;
            ind = new BufferedReader (new FileReader (dinfile));
            ok = false;
        }
        catch(FileNotFoundException e) {
            System.err.println("\nCan't find input file \"" + dinfile + "\"");
            System.err.print("Please enter a new file name> ");
            choix = cin.readLine();
        }
        }
        
        inv = null; 
        if ( option[7].equals("MVR*") 
         && option[8].equals("Yes") ) {
        choix = "input.v";
        ok = true;
        
        while (ok) {
            try {
            vinfile = choix;
            inv = new BufferedReader (new FileReader (vinfile));
            ok = false;
            }
            catch(FileNotFoundException e) {
            System.err.println("\nCan't find input file \"i" + vinfile + "\"");
            System.err.print("Please enter a new file name> ");
            choix = cin.readLine();
            }
        }
        }
    
        //System.err.print("\n");
    }

    if ( bionj ) out = new StringWriter();
    if ( nj )  out = new BufferedWriter(new FileWriter(new File (dinfile + "_nj.t")));
    if ( unj ) out = new BufferedWriter(new FileWriter(new File (dinfile + "_unj.t")));
    if ( mvr ) out = new BufferedWriter(new FileWriter(new File (dinfile + "_mvr.t")));


    
    while (multifile > 0) {

        multifile--;
        if (multifile > 1)
        System.err.print("\r                                   \r   Computing tree " + multifile);
        //if ( multifile == 1)
        //System.err.print("\r                                   ");

        // #################################################
        // ##### reading inpud.d and initializing data #####
        // #################################################
        ligne = ind.readLine();
        while ( ((ligne.trim()).hashCode() == 0)
            || (ligne.charAt(0) == '%')
            || (ligne.charAt(0) == '#') )
        ligne = ind.readLine();
        
        n = Integer.parseInt(ligne.trim());
        taxa = new String[n];
        dM = new Matrix(n);
        dMinit = new Matrix(n);
        hole = 0;

        for (int lig = 0 ; lig < n ; lig++) { // reading distance matrix
        ligne = ind.readLine().trim() + "  ";

        b = ligne.indexOf(' ');
        taxa[lig] = ligne.substring(0 , b);
        dMinit.setTaxon(taxa[lig] , lig);
            
        ligne = (ligne.substring(b)).trim();

        for (int col = 0 ; col < lig ; col++) {
            b = ligne.indexOf(' ');
            if (b == -1) dM.setValue(lig , col , (new Double(ligne.substring(0))).doubleValue());
            else dM.setValue(lig , col , (new Double(ligne.substring(0 , b))).doubleValue());
            if ( dM.getValue(lig , col) == 0.0) dM.setValue(lig , col , 0.00000000001);

            dMinit.setValue(lig , col , dM.getValue(lig , col));

            if ( dM.getValue(lig , col) == -99.0 ) {
            dM.setAbsence(lig , col);
            dMinit.setAbsence(lig , col);
            hole++;
            }
                    
            if (b != -1) ligne = (ligne.substring(b)).trim() + "  ";

            if (option[2].equals("Yes")) { // subreplicate format
            b = ligne.indexOf(' ');
            tempString = ligne.substring(0 , b);
            if (! tempString.equals("1")) {
                dM.setAbsence(lig , col);
                dMinit.setAbsence(lig , col);
                hole++;
            }
            ligne = (ligne.substring(b)).trim();
            }           
        }
        }
        if ( option[10].equals("No") ) dMinit = null;


        // #################################################
        // ##### reading inpud.v and initializing data #####
        // #################################################
        vM = null;
        if ( bionj 
         || mvr ) {
        vM = new Matrix(n);
        
        if ( bionj 
             || option[8].equals("No") ) { // the variance is computed from the distance

            for (int lig = 0 ; lig < n ; lig++)
            for (int col = 0 ; col < lig ; col++) {
                if ( option[7].equals("BioNJ*") )
                vM.setValue(lig , col , dM.getValue(lig , col));
                else
                vM.setValue(lig , col , Math.pow(dM.getValue(lig , col) , 2.0)); //1.823));
                
                if ( vM.getValue(lig , col) == 0.0 ) // no null variance
                vM.setValue(lig , col , 0.00000000001);
            }
        }
        else { // input.v for MVR*

            if ( mvr ) {

            ligne = inv.readLine();
            while ( ((ligne.trim()).hashCode() == 0)
                || (ligne.charAt(0) == '%') )
                ligne = inv.readLine();
            
            for (int lig = 0 ; lig < n ; lig++) { // reading variance matrix
                ligne = (inv.readLine()).trim() + "  ";
                b = ligne.indexOf(' ');
                ligne = (ligne.substring(b)).trim();
                
                for (int col = 0 ; col < lig ; col++) {
                b = ligne.indexOf(' ');
                if (b == -1)
                    vM.setValue(lig , col , (new Double(ligne.substring(0))).doubleValue());
                else
                    vM.setValue(lig , col , (new Double(ligne.substring(0 , b))).doubleValue());
                
                if ( vM.getValue(lig , col) == -99.0 )
                    vM.setAbsence(lig , col);
                
                if ( vM.getValue(lig , col) == 0.0 ) // no null variance
                    vM.setValue(lig , col , 0.00000000001);
                
                if (b != -1)
                    ligne = (ligne.substring(b)).trim() + "  ";
                
                if (option[2].equals("Yes")) { // subreplicate format
                    b = ligne.indexOf(' ');
                    tempString = ligne.substring(0 , b);
                    if (! tempString.equals("1"))
                    vM.setAbsence(lig , col);
                    ligne = (ligne.substring(b)).trim();
                }           
                } // for ( col
            } // for ( lig
            }
        } // else
        }

        // #############################
        // ##### initializing tree #####
        // #############################
        center = 2*n;                             // the center of the star tree is on the last line of the tree structure
        tT = new int[center + 1][];
        tD = new double[center + 1][];
        activeTaxa = new boolean[center + 1];
        tree2dist = new int[center + 1];
        
        for (int lig = 0 ; lig < center ; lig++) { // initializing the tree
        tT[lig] = new int[3];
        tD[lig] = new double[3];
        
        activeTaxa[lig] = false;
        tree2dist[lig] = -1;
        
        tT[lig][0] = -1;
        tT[lig][1] = -1;
        tT[lig][2] = -1;
        tD[lig][0] = 0; //Math.PI;
        tD[lig][1] = 0; //Math.PI;
        tD[lig][2] = 0; //Math.PI;
        }
        
        tT[center] = new int[n];                    // the center of the star tree
        tD[center] = new double[n];  
        
        for (int lig = 0 ; lig < n ; lig++) {       // the tree is initialized as a star
        tT[lig][0] = center;
        tD[lig][0] = Math.E;
        activeTaxa[lig] = true;
        tree2dist[lig] = lig;
        }
        for (int col = 0 ; col < n ; col++) {
        tT[center][col] = col;
        tD[center][col] = Math.E;
        }
        
        

        // ##################################
        // ##### initializing variables #####
        // ##################################
        sumR = new double[n];
        sumS = new double[n];

        _Qstar = null;
        _Sstar = null;
        if ( hole > 0 ) {
        _Qstar = new Matrix(n);
        _Sstar = new Matrix(n);
        }
        missing = new int[n];
        
        if ( p < 2 )
        p = 1;
        if ( p > (n*(n-1))/2 )
        p = (n*(n-1))/2;
        
        if ( p != 1 )
        tabQstar = new TreeSet<CriterionTuple>(compDesc);
        
        if ( unj )
        un = new double[n];
        mu = 0;
        w = new double[n];
        lambda = new double[n];
        for (int i = 0 ; i < n ; i++) {
        lambda[i] = 0.5;
        if ( unj )
            un[i] = 1.0;
            missing[i] = 0;
        }

        r = n;
        u = n;        // rank of the created nodes
        rankmin = 0;
        

        // ###################################################################################
        // ##### precomputing _Qstar and _Sstar tables if there exists missing distances #####
        // ###################################################################################
        if ( hole > 0 ) {

        for (int i = n ; --i >= 0 ; ) {
            //System.err.print("\r        \r" + i);
            for (int j = i ; --j >= 0 ; ) {
            if ( dM.exists(tree2dist[i] , tree2dist[j]) ) {
                _Qstar.setValue(tree2dist[i] , tree2dist[j] , 0);
                _Sstar.setValue(tree2dist[i] , tree2dist[j] , 0);
                
                for (int ii = n ; --ii >= 0 ; ) {
                if ( dM.exists(tree2dist[i] , tree2dist[ii])
                     && dM.exists(tree2dist[j] , tree2dist[ii]) ) {
                    _Qstar.add(tree2dist[i] , tree2dist[j] , dM.getValue(tree2dist[i] , tree2dist[ii])
                           + dM.getValue(tree2dist[j] , tree2dist[ii]));
                    _Sstar.add(tree2dist[i] , tree2dist[j] , 1.0);
                }
                } // ii
            }
            else {
                missing[tree2dist[i]]++;
                missing[tree2dist[j]]++;
            }
            } // j
        } // i
        }
        else {
        for (int z = u ; --z >= rankmin ; ) {
            if ( activeTaxa[z] ) {
            sumR[tree2dist[z]] = 0.0;
            
            for (int i = u ; --i >= rankmin ; ) {
                if ( activeTaxa[i] 
                 && (z != i) ) 
                sumR[tree2dist[z]] += dM.getValue(tree2dist[i] , tree2dist[z]);
            }
            }
        }
        }

        //System.err.print("##########\n");

        // ################################
        // ##### agglomerative scheme #####
        // ################################
        while (r > 3) {

        if ( p != 1 ) tabQstar.clear();
        else _Qmax = - Double.MAX_VALUE;
                    
        x = 0;
        y = 0;
    
        if ( hole == 0 ) {

            for (int i = u ; --i >= rankmin ; ) {
            if ( activeTaxa[i] ) {
                for (int j = i ; --j >= rankmin ; ) {
                if ( activeTaxa[j] ) {
                     //&& dM.exists(tree2dist[i] , tree2dist[j]) ) {

                    _Q = ((double) (2 - r)) * dM.getValue(tree2dist[i] , tree2dist[j])
                    + sumR[tree2dist[i]] 
                    + sumR[tree2dist[j]];

                    if ( p == 1 ) {
                    if ( _Q > _Qmax ) {
                        _Qmax = _Q;
                        x = j; y = i;
                    }
                    }
                    else { // p > 1

                    if ( (tabQstar.size() < p)
                         || ((tabQstar.size() == p)
                         && ((CriterionTuple) tabQstar.last()).Q <= _Q) )
                        tabQstar.add( new CriterionTuple( _Q , i , j) );
                    
                    if ( tabQstar.size() > p )
                        tabQstar.remove( tabQstar.last() );

                    }
                }
                }
            }
            }
        }

        else { // hole > 0

            for (int i = u ; --i >= rankmin ; ) {
            if ( activeTaxa[i] ) {
                for (int j = i ; --j >= rankmin ; ) {
                if ( activeTaxa[j] && dM.exists(tree2dist[i] , tree2dist[j]) ) {

                    _Qcur = _Qstar.getValue(tree2dist[i] , tree2dist[j]);
                    _Scur = _Sstar.getValue(tree2dist[i] , tree2dist[j]);

                    _Q = (_Qcur / ( _Scur - 2.0 )) - dM.getValue(tree2dist[i] , tree2dist[j]);

                    //_Q = (_Qcur / ( _Scur - 2.0 )) 
                    //    - ( ((double)(r-4))/((double)(r-2)) + 2.0/(_Scur-2.0) ) * dM.getValue(tree2dist[i] , tree2dist[j]);

                    if ( p == 1 ) {
                    if ( _Q > _Qmax ) { _Qmax = _Q; x = j; y = i; }
                    }
                    else { // p > 1
                    if ( (tabQstar.size() < p) || ((tabQstar.size() == p) && (tabQstar.last().Q <= _Q)) )
                        tabQstar.add( new CriterionTuple( _Q , i , j) );
                    if ( tabQstar.size() > p ) tabQstar.remove( tabQstar.last() );
                    }
                }
                }
            }
            }
        }

            
        if ( p > 1 ) { // testing p agglomerable pairs of taxa
                        
            dNmax =  0.0;
            cNmax = -999999.0;
            _Smax = - Double.MAX_VALUE;  // it is important that _Smax >> dNmax, cNmax
            _Qmax = - Double.MAX_VALUE;

            it = tabQstar.iterator();
            while ( it.hasNext() ) {
            ct = it.next();
            _j = ct.J;  // ij is the current neighborable pair to test
            _i = ct.I; 
                    
            if ( (_i != _j) && activeTaxa[_i] && activeTaxa[_j] ) {
                dNcur = 0.0; cNcur = 0.0; _Scur = 0.0;
                        
                for (int ii = u ; --ii >= rankmin ; ) {
                if ( activeTaxa[ii]
                     && (ii != _i)
                     && (ii != _j)
                     && ( (hole == 0) 
                      || (dM.exists(tree2dist[_i] , tree2dist[ii])
                          || dM.exists(tree2dist[_j] , tree2dist[ii])) ) ) {
                    for (int jj = ii ; --jj >= rankmin ; ) {
                    if ( activeTaxa[jj]
                         && (jj != _i)
                         && (jj != _j)
                         && ( (hole == 0) 
                          || (dM.exists(tree2dist[ii] , tree2dist[jj])
                              && ( dM.exists(tree2dist[_i] , tree2dist[jj]) 
                               || dM.exists(tree2dist[_j] , tree2dist[jj]) )) ) ) {
                        tempDouble = dM.getValue(tree2dist[_i] , tree2dist[_j])
                        + dM.getValue(tree2dist[ii] , tree2dist[jj]);
                                        
                        if ( (hole == 0)
                         || (dM.exists(tree2dist[_i] , tree2dist[ii])
                             && dM.exists(tree2dist[_j] , tree2dist[jj])) )
                        sum1 = dM.getValue(tree2dist[_i] , tree2dist[ii])
                            + dM.getValue(tree2dist[_j] , tree2dist[jj])
                            - tempDouble;
                        else
                        sum1 = Math.E;
                                        
                        if ( (hole == 0)
                         || (dM.exists(tree2dist[_i] , tree2dist[jj])
                             && dM.exists(tree2dist[_j] , tree2dist[ii])) )
                        sum2 = dM.getValue(tree2dist[_i] , tree2dist[jj])
                            + dM.getValue(tree2dist[_j] , tree2dist[ii])
                            - tempDouble;
                        else
                        sum2 = Math.E;
                                        
                        if ( (sum1 != Math.E) || (sum2 != Math.E) ) {
                        if ( (sum1 != Math.E) && (sum2 == Math.E) ) {
                            _Scur++; cNcur += sum1;
                            if ( sum1 >= 0 ) dNcur++;
                        }
                        if ( (sum1 == Math.E) && (sum2 != Math.E) ) {
                            _Scur++; cNcur += sum2;
                            if ( sum2 >= 0 ) dNcur++;
                        }

                        if ( (sum1 != Math.E) && (sum2 != Math.E) ) {
                            _Scur += 2.0; cNcur += sum1 + sum2;
                            if ( (sum1 >= 0) && (sum2 >= 0) ) dNcur += 2.0;
                        }
                        }
                    }
                    } // jj
                }
                } // ii
                        

                if ( _Scur != 0.0 ) {
                if ( dNcur / _Scur > dNmax / _Smax ) {
                    dNmax = dNcur; cNmax = cNcur ; _Smax = _Scur; _Qmax = ct.Q; x = _j; y = _i;
                }
                else {
                    if ( dNcur / _Scur == dNmax / _Smax ) {
                    if ( _Scur > _Smax ) {
                        dNmax = dNcur; cNmax = cNcur ; _Smax = _Scur; _Qmax = ct.Q; x = _j; y = _i;
                    }
                    else {
                        if ( _Scur == _Smax ) {
                        if ( missing[tree2dist[_i]] + missing[tree2dist[_j]] 
                             > missing[tree2dist[x]] + missing[tree2dist[y]] ) {
                            dNmax = dNcur; cNmax = cNcur ; _Smax = _Scur; _Qmax = ct.Q; x = _j; y = _i;
                        }
                        else {
                            if ( (missing[tree2dist[_i]] + missing[tree2dist[_j]] 
                                  == missing[tree2dist[x]] + missing[tree2dist[y]]) && (cNcur > cNmax) ) {
                                dNmax = dNcur; cNmax = cNcur ; _Smax = _Scur; _Qmax = ct.Q; x = _j; y = _i;
                            }
                        }
                        }
                    }
                    }
                }
                }
                if ( (x == y) && (ct.Q > _Qmax) ) {
                dNmax = dNcur; cNmax = cNcur; _Smax = _Scur; _Qmax = ct.Q; x = _j; y = _i;
                }
            }
            
            } // z
        }

        //System.err.println(x + " " + y);
        if ( x == y ) { // no possible agglomeration => exiting the agglomerative scheme
            System.err.println( x+ " "+y+ (char) 7);
            break;
        }
        if ( x == rankmin ) rankmin++;
        if ( y == rankmin ) // note that always x < y
            rankmin++;
                



        // ##########
        // agglomeration of taxa x and y
        // ##########
        b = -1;
        for (int col = 0 ; col < n ; col++) { // replacing x and y by u in the star tree
            if ( (tT[center][col] == x)
             || (tT[center][col] == y) ) {
            tT[center][col] = u;
            b = col;
            break;
            }
        }
        for (int col = b ; col < n ; col++) {
            if ( (tT[center][col] == x)
             || (tT[center][col] == y) ) {
            b = col;
            break;
            }
        }
        for (int col = b ; col < n-1 ; col++) {
            tT[center][col] = tT[center][col+1];
            tD[center][col] = tD[center][col+1];
            if ( tT[center][col] == -1 )
            break;
        }
        tT[center][n-1] = -1;
        tD[center][n-1] = 0; //Math.PI;



        // ##########
        // first updating of pQ and nS tables if there exists missing distances:
        // removing values corresponding to the xy pair
        // ##########
        if ( hole > 0 ) {

            for (int i = u ; --i >= rankmin ; ) {
            if ( activeTaxa[i]
                 && (i != x)
                 && (i != y) ) {
                for (int j = i ; --j >= rankmin ; ) {
                if ( activeTaxa[j] 
                     && (j != x)
                     && (j != y)
                     && dM.exists(tree2dist[i] , tree2dist[j]) ) {
                    
                    if ( dM.exists(tree2dist[i] , tree2dist[x])
                     && dM.exists(tree2dist[j] , tree2dist[x]) ) {
                    _Qstar.add( tree2dist[i] , tree2dist[j] , - (dM.getValue(tree2dist[i] , tree2dist[x])
                                             + dM.getValue(tree2dist[j] , tree2dist[x])) );
                    _Sstar.add( tree2dist[i] , tree2dist[j] , -1.0 );
                    }

                    if ( dM.exists(tree2dist[i] , tree2dist[y])
                     && dM.exists(tree2dist[j] , tree2dist[y]) ) {
                    _Qstar.add( tree2dist[i] , tree2dist[j] , - (dM.getValue(tree2dist[i] , tree2dist[y])
                                             + dM.getValue(tree2dist[j] , tree2dist[y])) );
                    _Sstar.add(tree2dist[i] , tree2dist[j] , -1.0);
                    }
                }
                } // j
            }
            } // i
        }
        else {
            
            for (int z = u ; --z >= rankmin ; ) {
            if ( activeTaxa[z] ) 
                sumR[tree2dist[z]] -= dM.getValue(tree2dist[z] , tree2dist[x]) + dM.getValue(tree2dist[z] , tree2dist[y]);
            }
        }

        // ##########
        // computing w_is and lambda_is
        // ##########
        if ( nj ) {

            tempDouble = 0.0;
            for (int i = u ; --i >= rankmin ; ) {
            if ( activeTaxa[i] 
                 && ( i != x )
                 && ( i != y ) 
                 && ( (hole == 0)
                  || (dM.exists(tree2dist[x] , tree2dist[i])
                      && dM.exists(tree2dist[y] , tree2dist[i])) ) )
                tempDouble++;
            }
            tempDouble = Math.pow(2.0 * tempDouble , -1);

            for (int i = u ; --i >= rankmin ; ) { // computing w_is and lambda_is
            if ( activeTaxa[i] ) {
                w[tree2dist[i]] = 0.0;
                
                if ( hole == 0 ) {
                w[tree2dist[i]] = tempDouble;
                lambda[tree2dist[i]] = 0.5;
                }
                else {
                if ( dM.exists(tree2dist[x] , tree2dist[i])
                     && dM.exists(tree2dist[y] , tree2dist[i]) ) {
                    w[tree2dist[i]] = tempDouble;
                    lambda[tree2dist[i]] = 0.5;
                }
                else {
                    if ( (! dM.exists(tree2dist[x] , tree2dist[i]))
                     && dM.exists(tree2dist[y] , tree2dist[i]) )
                    lambda[tree2dist[i]] = 0.0; 
                    else {
                    if ( dM.exists(tree2dist[x] , tree2dist[i])
                         && (! dM.exists(tree2dist[y] , tree2dist[i])) )
                        lambda[tree2dist[i]] = 1.0;
                    }
                }
                }
            }
            }
        }
        else {
            if ( unj ) {

            tempDouble = 0.0;
            for (int i = u ; --i >= rankmin ; ) {
                if ( activeTaxa[i] 
                 && ( i != x )
                 && ( i != y ) 
                 && ( (hole == 0)
                      || (dM.exists(tree2dist[x] , tree2dist[i])
                      && dM.exists(tree2dist[y] , tree2dist[i])) ) )
                tempDouble += un[tree2dist[i]];
            }
            tempDouble = Math.pow(2.0 * tempDouble , -1);
            tempDouble2 = un[tree2dist[x]] 
                / ( un[tree2dist[x]] + un[tree2dist[y]] );

            for (int i = u ; --i >= rankmin ; ) { // computing w_is and lambda_is
                if ( activeTaxa[i] ) {
                w[tree2dist[i]] = 0.0;

                if ( hole == 0 ) {
                    w[tree2dist[i]] = un[tree2dist[i]] * tempDouble;
                    lambda[tree2dist[i]] = tempDouble2;
                }
                else {
                    if ( dM.exists(tree2dist[x] , tree2dist[i])
                     && dM.exists(tree2dist[y] , tree2dist[i]) ) {
                    w[tree2dist[i]] = un[tree2dist[i]] * tempDouble;
                    lambda[tree2dist[i]] = tempDouble2;
                    }
                    else {
                    if ( (! dM.exists(tree2dist[x] , tree2dist[i]))
                         && dM.exists(tree2dist[y] , tree2dist[i]) )
                        lambda[tree2dist[i]] = 0.0; 
                    else {
                        if ( dM.exists(tree2dist[x] , tree2dist[i])
                         && (! dM.exists(tree2dist[y] , tree2dist[i])) )
                        lambda[tree2dist[i]] = 1.0;
                    }
                    }
                }
                }
            }
            }
            else {

            if ( bionj ) {

                tempDouble = 0.0;
                tempDouble2 = 0.0;
                for (int  i = u ; --i >= rankmin ; ) {
                if ( activeTaxa[i] 
                     && ( i != x )
                     && ( i != y ) 
                     && dM.exists(tree2dist[x] , tree2dist[i])
                     && dM.exists(tree2dist[y] , tree2dist[i]) ) {
                    tempDouble++;
                    tempDouble2 += vM.getValue(tree2dist[y] , tree2dist[i])
                    - vM.getValue(tree2dist[x] , tree2dist[i]);
                }
                }

                for (int i = u ; --i >= rankmin ; ) { // computing w_is and lambda_is
                if ( activeTaxa[i] ) {
                    w[tree2dist[i]] = 0.0;

                    if ( dM.exists(tree2dist[x] , tree2dist[i])
                     && dM.exists(tree2dist[y] , tree2dist[i]) )
                    w[tree2dist[i]] = Math.pow(2.0 * tempDouble , -1);

                    if ( (! dM.exists(tree2dist[x] , tree2dist[i]))
                     && dM.exists(tree2dist[y] , tree2dist[i]) )
                    lambda[tree2dist[i]] = 0.0; 

                    if ( dM.exists(tree2dist[x] , tree2dist[i])
                     && (! dM.exists(tree2dist[y] , tree2dist[i])) )
                    lambda[tree2dist[i]] = 1.0;

                    if ( dM.exists(tree2dist[x] , tree2dist[i])
                     && dM.exists(tree2dist[y] , tree2dist[i]) )
                    lambda[tree2dist[i]] = 0.5 
                        + tempDouble2 
                        * Math.pow(2.0 * tempDouble * vM.getValue(tree2dist[x] , tree2dist[y]) , -1);

                    if ( ((vM.getValue(tree2dist[x] , tree2dist[y]) != 0.0)
                      && (lambda[tree2dist[i]] > 1.0))
                     || ((vM.getValue(tree2dist[x] , tree2dist[y]) == 0.0) // infinite lambda
                         && tempDouble2 > 0) )
                    lambda[tree2dist[i]] = 1.0;

                    if ( ((vM.getValue(tree2dist[x] , tree2dist[y]) != 0.0)
                      && (lambda[tree2dist[i]] < 0.0))
                     || ((vM.getValue(tree2dist[x] , tree2dist[y]) == 0.0) // infinite lambda
                         && tempDouble2 < 0) )
                    lambda[tree2dist[i]] = 0.0; 
                }
                }
            }
            else {  // mvr

                mu = 0; // computing mu
                for (int i = u  ; --i >= rankmin ; ) {
                if ( activeTaxa[i] 
                     && ( i != x )
                     && ( i != y) 
                     && dM.exists(tree2dist[i] , tree2dist[x])
                     && dM.exists(tree2dist[i] , tree2dist[y]) )
                    mu += Math.pow( vM.getValue(tree2dist[x] , tree2dist[i]) + vM.getValue(tree2dist[y] , tree2dist[i]) , -1 );
                }
                mu = Math.pow( 2.0 * mu , -1 );

                for (int i = u ; --i >= rankmin ; ) { // computing w_is and lambda_is
                if ( activeTaxa[i] ) {
                    w[tree2dist[i]] = 0.0;
                    
                    if ( dM.exists(tree2dist[x] , tree2dist[i])
                     && dM.exists(tree2dist[y] , tree2dist[i]) )
                    w[tree2dist[i]] = mu / (vM.getValue(tree2dist[x] , tree2dist[i]) + vM.getValue(tree2dist[y] , tree2dist[i]));
                            
                    if ( (! dM.exists(tree2dist[x] , tree2dist[i]))
                     && dM.exists(tree2dist[y] , tree2dist[i]) )
                    lambda[tree2dist[i]] = 0.0; 

                    if ( dM.exists(tree2dist[x] , tree2dist[i])
                     && (! dM.exists(tree2dist[y] , tree2dist[i])) )
                    lambda[tree2dist[i]] = 1.0;

                    if ( dM.exists(tree2dist[x] , tree2dist[i])
                     && dM.exists(tree2dist[y] , tree2dist[i]) )
                    lambda[tree2dist[i]] = vM.getValue(tree2dist[y] , tree2dist[i])
                        / ( vM.getValue(tree2dist[x] , tree2dist[i]) + vM.getValue(tree2dist[y] , tree2dist[i]) );
                }
                }
            } // else
            } // else
        } // else
        
  
            
        // ##########
        // branch length of T_xu and T_yu
        // ##########
        txu = 0;
        tempDouble = 0;
        for (int i = u ; --i >= rankmin ; ) {
            if ( activeTaxa[i] 
             && ( i != x )
             && ( i != y ) 
             && ( (hole ==0)
                  || (dM.exists(tree2dist[i] , tree2dist[x])
                  && dM.exists(tree2dist[i] , tree2dist[y])) ) )
            txu += w[tree2dist[i]]
                * ( dM.getValue(tree2dist[x] , tree2dist[i]) - dM.getValue(tree2dist[y] , tree2dist[i]) ); 
        }
        txu += 0.5 * dM.getValue(tree2dist[x] , tree2dist[y]);

        tyu = dM.getValue(tree2dist[x] , tree2dist[y]) - txu;

        if (x < n) { // if x is a leaf
            tT[x][0] = u; // creating edge T_xu
            tD[x][0] = txu;
        }
        else {
            tT[x][2] = u; // creating edge T_xu
            tD[x][2] = txu;
        }

        if (y < n) { // if y is a leaf
            tT[y][0] = u; // creating edge T_yu
            tD[y][0] = tyu;
        }
        else {
            tT[y][2] = u; // creating edge T_yu
            tD[y][2] = tyu;
        }

        tT[u][0] = x; // creating node u
        tD[u][0] = txu;
        tT[u][1] = y;
        tD[u][1] = tyu;
        tT[u][2] = center; // connecting node u in the star tree
        tD[u][2] = Math.E;

            
        // ##########
        // reduction of the distance matrix dM (and eventually the variance matrix vM)
        // ##########
        activeTaxa[u] = true;
        tree2dist[u] = tree2dist[x];
        missing[tree2dist[u]] = 0;

        if ( unj ) un[tree2dist[u]] = un[tree2dist[x]] + un[tree2dist[y]];

        //txu = 0;
        //tyu = 0;

        for (int i = u ; --i >= rankmin ; ) {
            if ( activeTaxa[i]
             && (x != i)
             && (y != i) ) { // && (u != i) )

            if ( ! dM.exists(tree2dist[x] , tree2dist[i]) ) {
                 missing[tree2dist[i]]--;
                 hole--;
            }
            if ( ! dM.exists(tree2dist[y] , tree2dist[i]) ) {
                 missing[tree2dist[i]]--;
                 hole--;
            }

            if ( (! dM.exists(tree2dist[x] , tree2dist[i]))
                 && (! dM.exists(tree2dist[y] , tree2dist[i])) ) {
                dM.setAbsence(tree2dist[u] , tree2dist[i]);
                missing[tree2dist[u]]++;
                missing[tree2dist[i]]++;
                hole++;

                if ( bionj
                 || mvr )
                vM.setAbsence(tree2dist[u] , tree2dist[i]);
            }
            else {
                dM.setValue(tree2dist[u] , tree2dist[i] , 
                    lambda[tree2dist[i]] * dM.getValue(tree2dist[x] , tree2dist[i])
                    + (1.0 - lambda[tree2dist[i]]) * dM.getValue(tree2dist[y] , tree2dist[i])
                    - lambda[tree2dist[i]] * txu
                    - (1.0 - lambda[tree2dist[i]]) * tyu);
                
                if ( bionj ) // reduction of vM
                vM.setValue(tree2dist[u] , tree2dist[i] , 
                        lambda[tree2dist[i]] * vM.getValue(tree2dist[x] , tree2dist[i])
                        + (1.0 - lambda[tree2dist[i]]) * vM.getValue(tree2dist[y] , tree2dist[i])
                        - lambda[tree2dist[i]] * (1.0 - lambda[tree2dist[i]]) * vM.getValue(tree2dist[x] , tree2dist[y]));
                
                if ( mvr ) { // reduction of vM
                if ( lambda[tree2dist[i]] == 1.0 )
                    vM.setValue(tree2dist[u] , tree2dist[i] , 
                        vM.getValue(tree2dist[x] , tree2dist[i]));
                            
                if ( lambda[tree2dist[i]] == 0.0 )
                    vM.setValue(tree2dist[u] , tree2dist[i] , 
                        vM.getValue(tree2dist[y] , tree2dist[i]));
                            
                if ( (lambda[tree2dist[i]] != 0.0)
                     && (lambda[tree2dist[i]] != 1.0) )
                    vM.setValue( tree2dist[u] , tree2dist[i] , 
                         ( vM.getValue(tree2dist[x] , tree2dist[i]) * vM.getValue(tree2dist[y] , tree2dist[i]) )
                         / ( vM.getValue(tree2dist[x] , tree2dist[i]) + vM.getValue(tree2dist[y] , tree2dist[i]) ) );
                }
            }
            }
        }
        activeTaxa[x] = false;
        tree2dist[x] = -1;
        activeTaxa[y] = false;
        tree2dist[y] = -1;

        

        // ##########
        // updating _Qstar and _Sstar tables if there exists missing distances :
        // ##########
        if ( hole > 0 ) {
            for (int i = u ; --i >= rankmin ; ) {
            if ( activeTaxa[i] ) {

                // ##########
                // second updating of _Qstar and _Sstar tables :
                // computing values corresponding to the ui pairs, for every i
                // ##########
                if ( dM.exists(tree2dist[i] , tree2dist[u]) ) {
                _Qstar.setValue( tree2dist[i] , tree2dist[u] , 0);
                _Sstar.setValue( tree2dist[i] , tree2dist[u] , 0);

                for (int ii = u+1 ; --ii >= rankmin ; ) { // u+1 car il faut compter cette valeur
                    if ( activeTaxa[ii] //&& (ii != i)
                     && ( dM.exists(tree2dist[i] , tree2dist[ii]) 
                          && dM.exists(tree2dist[u] , tree2dist[ii]) ) ) {
                    _Qstar.add(tree2dist[i] , tree2dist[u] , dM.getValue(tree2dist[i] , tree2dist[ii])
                           + dM.getValue(tree2dist[u] , tree2dist[ii]));
                    _Sstar.add(tree2dist[i] , tree2dist[u] , 1.0);
                    }
                } // ii
                }

                // ##########
                // third updating of _Qstar and _Sstar tables :
                // adding values corresponding to the ui and uj pairs, for every i j 
                // ##########
                for (int j = i ; --j >= rankmin ; ) {
                if ( activeTaxa[j] 
                     && dM.exists(tree2dist[i] , tree2dist[j])
                     && dM.exists(tree2dist[i] , tree2dist[u])
                     && dM.exists(tree2dist[j] , tree2dist[u]) ) {
                    _Qstar.add(tree2dist[i] , tree2dist[j] , dM.getValue(tree2dist[i] , tree2dist[u])
                           + dM.getValue(tree2dist[j] , tree2dist[u]));
                    _Sstar.add(tree2dist[i] , tree2dist[j] , 1.0);
                }
                } // j
            }
            } // i
        }
        else {
            
            for (int z = u ; --z >= rankmin ; ) {
            if ( activeTaxa[z] ) 
                sumR[tree2dist[z]] += dM.getValue(tree2dist[z] , tree2dist[u]);
            }
            sumR[tree2dist[u]] = 0;
            for (int i = u ; --i >= rankmin ; ) {
            if ( activeTaxa[i] 
                 && (u != i) ) {
                sumR[tree2dist[u]] += dM.getValue(tree2dist[i] , tree2dist[u]);
            }
            }
        }

        
        u++;
        r--;

        }  // end of the agglomerative scheme


        // ##########
        // last agglomeration(s)
        // ##########
        if ( r == 3 ) { // binary tree
        u--;

        x = 0;
        y = 0;
        b = 0;
        while ( tT[center][b] != u )
            b++;
        if (b == 0) { x = 1; y = 2; }
        else {
            if (b == 1) { x = 0; y = 2; }
            else { x = 0; y = 1; }
        }

        if ( (hole == 0)
             || (dM.exists(tree2dist[u] , tree2dist[tT[center][x]])
             && dM.exists(tree2dist[u] , tree2dist[tT[center][y]])
             && dM.exists(tree2dist[tT[center][x]] , tree2dist[tT[center][y]])) ) { // the three last distances exist
                
            if ( mvr ) {
            mu = 0; // computing mu
            for (int i = 0 ; i < center + 1 ; i++) {
                if ( activeTaxa[i] 
                 && ( i != tT[center][x] )
                 && ( i != tT[center][y] ) 
                 && dM.exists(tree2dist[i] , tree2dist[tT[center][x]])
                 && dM.exists(tree2dist[i] , tree2dist[tT[center][y]]) )
                mu += Math.pow( vM.getValue(tree2dist[tT[center][x]] , tree2dist[i]) 
                        + vM.getValue(tree2dist[tT[center][y]] , tree2dist[i]) , -1 );
            }
            mu = Math.pow( 2.0 * mu , -1 );
            
            for (int i = 0 ; i < center + 1 ; i++) { // computing w_is
                if ( activeTaxa[i] ) {
                w[tree2dist[i]] = 0.0;
                
                if ( dM.exists(tree2dist[tT[center][x]] , tree2dist[i])
                     && dM.exists(tree2dist[tT[center][y]] , tree2dist[i]) )
                    w[tree2dist[i]] = mu 
                    / ( vM.getValue(tree2dist[tT[center][x]] , tree2dist[i]) 
                        + vM.getValue(tree2dist[tT[center][y]] , tree2dist[i]) );
                }
            }
            }
            
            txu = 0;
            tempDouble = 0;
            for (int i = 0 ; i < center + 1 ; i++) {
            if ( activeTaxa[i] 
                 && ( i != tT[center][x] )
                 && ( i != tT[center][y] ) 
                 && dM.exists(tree2dist[i] , tree2dist[tT[center][x]])
                 && dM.exists(tree2dist[i] , tree2dist[tT[center][y]]) ) {
                if ( nj
                 || bionj ) {
                txu += 0.5 
                    * ( dM.getValue(tree2dist[tT[center][x]] , tree2dist[i]) 
                    - dM.getValue(tree2dist[tT[center][y]] , tree2dist[i]) );
                tempDouble++;
                }
                else {
                if ( unj )
                    txu += 0.5 * un[tree2dist[i]] 
                    * ( dM.getValue(tree2dist[tT[center][x]] , tree2dist[i]) 
                        - dM.getValue(tree2dist[tT[center][y]] , tree2dist[i]) );
                else // mvr
                    txu += w[tree2dist[i]]
                    * ( dM.getValue(tree2dist[tT[center][x]] , tree2dist[i]) 
                        - dM.getValue(tree2dist[tT[center][y]] , tree2dist[i]) );
                }
            }
            }
            
            if ( nj
             || bionj )
            txu = 0.5 * dM.getValue(tree2dist[tT[center][x]] , tree2dist[tT[center][y]]) 
                + txu / tempDouble;
            else {
            if ( unj )
                txu = 0.5 * dM.getValue(tree2dist[tT[center][x]] , tree2dist[tT[center][y]]) 
                + txu / ( ((double) n) - un[tree2dist[tT[center][x]]] - un[tree2dist[tT[center][y]]] );
            else // mvr
                txu += 0.5 * dM.getValue(tree2dist[tT[center][x]] , tree2dist[tT[center][y]]);
            }
            //System.err.println(txu);
            
            tD[center][x] = txu; //dM[tree2dist[tT[2*n][x]]][tree2dist[u]] - tD[u][2];
            if (tT[center][x] < n)
            tD[tT[center][x]][0] = tD[center][x];
            else
            tD[tT[center][x]][2] = tD[center][x];
            
            tD[center][y] = dM.getValue(tree2dist[tT[center][x]] , tree2dist[tT[center][y]]) - txu;
            if (tT[center][y] < n)
            tD[tT[center][y]][0] = tD[center][y];
            else
            tD[tT[center][y]][2] = tD[center][y];
            
            tD[u][2] = dM.getValue(tree2dist[tT[center][x]] , tree2dist[u]) - txu;
            tD[center][b] = tD[u][2];
            
        }
        }
        else { // non binary tree => artificial agglomeration
        
        while ( r > 3 ) {
            x = tT[center][0];
            y = tT[center][1];
            
            // ##########
            // agglomeration of taxa x and y
            // ##########
            b = -1;
            for (int col = 0 ; col < n ; col++) { // replacing x and y by u in the star tree
            if ( (tT[center][col] == x)
                 || (tT[center][col] == y) ) {
                tT[2*n][col] = u;
                b = col;
                break;
            }
            }
            for (int col = b ; col < n ; col++) {
            if ( (tT[center][col] == x)
                 || (tT[center][col] == y) ) {
                b = col;
                break;
            }
            }
            for (int col = b ; col < n-1 ; col++) {
            tT[center][col] = tT[center][col+1];
            tD[center][col] = tD[center][col+1];
            if ( tT[center][col] == -1 )
                break;
            }
            tT[center][n-1] = -1;
            tD[center][n-1] = 0; //Math.PI;

            txu = Math.E;
            tyu = Math.E;

            if (x < n) { // if x is a leaf
            tT[x][0] = u; // creating edge T_xu
            tD[x][0] = txu;
            }
            else {
            tT[x][2] = u; // creating edge T_xu
            tD[x][2] = txu;
            }

            if (y < n) { // if y is a leaf
            tT[y][0] = u; // creating edge T_yu
            tD[y][0] = tyu;
            }
            else {
            tT[y][2] = u; // creating edge T_yu
            tD[y][2] = tyu;
            }

            tT[u][0] = x; // creating node u
            tD[u][0] = txu;
            tT[u][1] = y;
            tD[u][1] = tyu;
            tT[u][2] = center; // connecting node u in the star tree
            tD[u][2] = Math.E;

            activeTaxa[u] = true;
            tree2dist[u] = tree2dist[x];
            activeTaxa[x] = false;
            tree2dist[x] = -1;
            activeTaxa[y] = false;
            tree2dist[y] = -1;

            u++;
            r--;
        }
        }

        // ##############################################
        // ##### writing  the tree in newick format #####
        // ##############################################
        root = tT[0][0];
        if ( ! outgroup.equals("") ) {
            for (int lig = 0 ; lig < n ; lig++) {
            if ( taxa[lig].equals(outgroup) ) {
                root = tT[lig][0];  break;
            }
            }
        }

        nodeFixed = new BitSet(center + 1);
        nodeFixed.set(root);

        newick = ( (new StringBuffer("(")).append(makeTree(tT[root][0])).append(',').append(makeTree(tT[root][1])).append(',').append(makeTree(tT[root][2])).append(");") ).toString();
        


        // ##############################################
        // ##### cleaning the tree in newick format #####
        // ##############################################
        b = 0;
        while ( newick.charAt(b) != ';' ) {
            if ( newick.charAt(b) == ':' ) {

            if ( option[6].equals("No")                                       // non binary tree
                 &&  ( ((newick.substring(b , b + 9)).equals(":2.718281"))    // unlengthed branch
                   || ( option[5].equals("No")
                    && (newick.charAt(b+1) == '-') ) ) ) {

                if ( newick.charAt(b-1) == ')' ) { // creating polytomy
                b2 = b - 2;
                r = 1;
                while ( r != 0 ) {
                    b2--;
                    if ( newick.charAt(b2) == ')' )
                    r++;
                    if ( newick.charAt(b2) == '(' )
                    r--;
                }
                newick = newick.substring(0 , b2)
                    + newick.substring(b2 + 1);
                b -= 2;
                b2 = b + 1;
                while ( (newick.charAt(b2) != ')')
                    && (newick.charAt(b2) != ',') )
                    b2++;
                newick = newick.substring(0 , b)
                    + newick.substring(b2);
                }
                else { // it's a leaf with no branch length
                newick = newick.substring(0 , b) 
                    + ":0.000000" 
                    +  newick.substring(b + 9);
                }
            }

            if ( b < newick.length() - 9 )
                if ( ((newick.substring(b , b + 9)).equals(":2.718281"))     // unlengthed branch
                 || ((newick.substring(b , b + 9)).equals(":-0.00000"))  // -0.0 length branch
                 || ( (option[5].equals("No"))                           // no negative branch length
                      && ((newick.charAt(b + 1)) == '-')) ) {
                
                if ( newick.indexOf(',' , b) < 0 )
                    b2 = newick.indexOf(')' , b);
                else 
                    b2 = Math.min(newick.indexOf(')' , b) , newick.indexOf(',' , b));
                    

                newick = newick.substring(0 , b) + ":0.000000" +  newick.substring(b2);
                }
            }
                
            b++;
        }

        if ( option[10].equals("Yes") ) newick = setConfidenceValues( dMinit , newick );
        //System.err.println(newick);
        out.write(newick + "\n");


        } // end of while(multifile


    String outS = out.toString();
    out.close();
//    System.err.println("");
    return outS;

    }




    //##### recursive function to compute the newick tree representation #####
     StringBuffer makeTree( int node ) {
    int x,y,u;
    nodeFixed.set(node);  // to mark the node
    if ( node < n )
        return new StringBuffer(taxa[node] + ":" + double2String(tD[node][0] , 6));
    else {
        if ( nodeFixed.get(tT[node][0]) ) {
        u = 0; // tT[node][0];
        x = tT[node][1]; y = tT[node][2];
        }
        else {
        if ( nodeFixed.get(tT[node][1]) ) {
            u = 1; // tT[node][1];
            x = tT[node][0]; y = tT[node][2];
        }
        else {
            u = 2; // tT[node][2];
            x = tT[node][0]; y = tT[node][1];
        }
        }
        return (new StringBuffer("(")).append(makeTree(x)).append(',').append(makeTree(y)).append("):").append(double2String(tD[node][u] , 6));
    }
    }



    //##### function to compute a clean String representation of a double #####
     StringBuffer double2String(double d , int limit) {
    return new StringBuffer(String.format(Locale.ENGLISH , "%." + limit + "f" , new Double(d)));
    }


    //##### function to compute distance-based confidence values at branches #####
     String setConfidenceValues(Matrix dm , String tr) {
    int _b1, _b2, _b3, _b4, _b5;
    ArrayList<String> _tax, _l1, _l2, _l3, _l4;
    StringBuffer _sb;
    int _x, _y, _u, _v, __x, __y, __u, __v;
    double _q_ok, _q_total;


    _tax = new ArrayList<String>(0);
    _b1 = -1; while ( ++_b1 < dm.getSize() ) _tax.add( dm.getTaxon(_b1) );
    
    _sb = new StringBuffer(tr); // rooting
    if ( _sb.charAt(1) == '(' ) {
        _b1 = parenthesisForward( _sb , 1 );
        _b1 = _sb.indexOf( "," , _b1 );
        _b1 = _sb.indexOf( "," , _b1+1 );
    }
    else {
        _b1 = _sb.indexOf("," , 1);
        _b1++;
        if ( _sb.charAt(_b1) == '(' ) {
        _b1 = parenthesisForward( _sb , _b1 );
        _b1 = _sb.indexOf( "," , _b1 );
        }
        else 
        _b1 = _sb.indexOf("," , _b1 );
    }
    _sb = _sb.insert( _b1 , "):-99" ).insert( 1 , "(" );
    //System.err.println( "\n" + _sb.toString() );
    
    _b1 = 0;
    while ( true ) {                           // exploration of all quartet taxon sets
        _b1 = _sb.indexOf("(" , ++_b1);
        if ( _b1 == -1 )
        break;
        _b2 = parenthesisForward( _sb , _b1 ); // the index of ')' corresponding to the '(' at _b1
        _l2 = taxonSet(_sb , _b1 , _b2);       // _l2 contains all taxa between _b1 and _b2
        _l4 = new ArrayList<String>(0);
        _b5 = -1;
        while ( ++_b5 < _tax.size() ) {        // _l4 = _tax - _l2
        if ( ! _l2.contains(_tax.get(_b5)) )
            _l4.add(_tax.get(_b5));
        }
        _l1 = new ArrayList<String>(0);        // _l2 => _l1 U _l2 thanks to the bipartition of taxa inside _l2
        if ( _l2.size() == 2 ) 
        _l1.add(_l2.remove(0));
        else {
        _b3 = _sb.indexOf("(" , _b1+1);
        _b4 = parenthesisForward( _sb , _b3 );
        _l1 = taxonSet( _sb , _b3 , _b4 );
        _b5 = -1;
        while ( ++_b5 < _l1.size() )
            _l2.remove( _l1.get(_b5) );
        }
        _b5 = 1;
        _b3 = _b1;                             // _b3 is the '(' index of the first clade containing _b1-_b2
        while ( _b5 != 0 ) {
        if ( _sb.charAt(--_b3) == ')' )
            _b5++;
        if ( _sb.charAt(_b3) == '(' )
            _b5--;
        }
        _b4 = parenthesisForward( _sb , _b3 );
        _l3 = taxonSet(_sb , _b3 , _b4);
        _b5 = -1;
        while ( ++_b5 < _l1.size() )
        _l3.remove(_l1.get(_b5));
        _b5 = -1;
        while ( ++_b5 < _l2.size() )
        _l3.remove(_l2.get(_b5));
        _b5 = -1;
        while ( ++_b5 < _l3.size() )
        _l4.remove(_l3.get(_b5));

        if ( _l4.isEmpty() && (_b3 == 0) && ( _b1 > 1) ) {  // special case for root trifurcation
        if ( _l3.size() == 2 ) 
            _l4.add(_l3.remove(0));
        else {
            _b3 = _sb.indexOf("(" , 2);
            _b4 = parenthesisForward( _sb , _b3 );
            _l4 = taxonSet( _sb , _b3 , _b4 );
            _b5 = -1;
            while ( ++_b5 < _l4.size() )
            _l3.remove( _l4.get(_b5) );
        }
        }

        if ( (! _l1.isEmpty()) && (! _l2.isEmpty()) && (! _l3.isEmpty()) && (! _l4.isEmpty()) ) {
        // here we have index _b2 that corresponds to quartet _l1_l2 | _l3_l4
        _q_ok = 0; _q_total = 0;
        _x = -1;
        while ( ++_x < _l1.size() ) {
            __x = _tax.indexOf(_l1.get(_x));
            _y = -1;
            while ( ++_y < _l2.size() ) {
            __y = _tax.indexOf(_l2.get(_y));
            if ( dm.exists(__x , __y) ) {
                _u = -1;
                while ( ++_u < _l3.size() ) {
                __u = _tax.indexOf(_l3.get(_u));
                if ( dm.exists(__u , __x) && dm.exists(__u , __y) ) {
                    _v = -1;
                    while ( ++_v < _l4.size() ) {
                    __v = _tax.indexOf(_l4.get(_v));
                    if ( dm.exists(__v , __x) && dm.exists(__v , __y) && dm.exists(__v , __u) ) {
                        if ( dm.getValue( __x , __y) + dm.getValue( __u , __v) 
                         <= Math.min( dm.getValue( __x , __u) + dm.getValue( __y , __v) ,
                                  dm.getValue( __x , __v) + dm.getValue( __y , __u) ) ) _q_ok++;
                        _q_total++;
                    }
                    }
                }
                }
            }
            }
        }
        //System.err.println(_q_ok + " " + _q_total);
        _sb = _sb.insert( _b2+1 , ((int) (100 * _q_ok / _q_total)) );
        }
        else {
        //System.err.println("\n" + _l1.toString());
        //System.err.println(_l2.toString());
        //System.err.println(_l3.toString());
        //System.err.println(_l4.toString());
        //System.err.println(_b1 + " " + _b3);
        }   
    }

    _b1 = parenthesisForward( _sb , 1 ); // unrooting
    _b2 = _sb.indexOf("," , _b1);
    _sb = _sb.delete(_b1 , _b2).deleteCharAt(1);

    return _sb.toString();
    }
    
    //##### returns the taxon set from newick that are from start to end #####
     ArrayList<String> taxonSet (StringBuffer newick , int start , int end) {
    ArrayList<String> set = new ArrayList<String>(0);
    int __b2;
    int __b1 = start - 1;
    while ( ++__b1 < end ) {
        try {
        if ( ((newick.charAt(__b1) == '(') || (newick.charAt(__b1) == ',')) && (newick.charAt(__b1+1) != '(') ) {
            __b2 = __b1;
            while ( newick.charAt(++__b2) != ':' ) {}
            set.add( newick.substring(__b1+1 , __b2) );
            __b1 = __b2;
        }
        }
        catch(NullPointerException e) {
        break;
        }
    }
    return set;
    }


    //##### returns the index of the closing parenthesis corresponding to the open parenthesis at the specified index #####
     int parenthesisForward(StringBuffer newick , int start) {
    int ___b, ___p;
    ___p = 1; ___b = start;
    while ( ___p != 0 ) {
        if ( newick.charAt(++___b) == '(' ) ___p++;
        if ( newick.charAt(___b) == ')' ) ___p--;
    }
    return ___b;
    }






    //##### class to store a pairwise distance and its two corresponding taxon numbers #####
     class CriterionTuple {
    double Q; int I, J;
    CriterionTuple(double q , int iii , int jjj) { this.Q = q; this.I = iii; this.J = jjj; }
    }



    //##### class for symetric matrices that shares only the lower triangular entries with a boolean table for missing values #####
    class Matrix {

    float[][] matrice;          // distance matrix
    boolean[][] presence;       // presence matrix
    int taille;                 // matrix size
    String[] taxa;              // taxa list
    boolean withTaxon;          // true if there are taxon names
    
    double d;
    boolean b;
    String s;
    
    Matrix (int n) {
        taille = n;
        matrice = new float[n][];
        presence = new boolean[n][];
        taxa = null; 
        withTaxon = false;
        for (int i = 0 ; i < n ; i++) { matrice[i] = new float[i]; presence[i] = new boolean[i]; }
        
    }
    
    public void setValue (int lig , int col , double val) {
        if ( lig != col ) {
        if ( col < lig ) { matrice[lig][col] = (float) val; presence[lig][col] = true; }
        else { matrice[col][lig] = (float) val; presence[col][lig] = true; }
        }
    }
    
    public double getValue (int lig , int col)
    {
        if ( lig != col ) {
        if ( col < lig )
            return (double) matrice[lig][col];
        else // if ( lig < col )
            return (double) matrice[col][lig];
        }
        return 0;
    }
    
    public void add (int lig , int col , double val) {
        if ( lig != col ) {
        if ( col < lig ) {
            if ( presence[lig][col] ) matrice[lig][col] += (float) val;
            else matrice[lig][col] = (float) val;
        } else {
            if ( presence[col][lig] ) matrice[col][lig] += (float) val;
            else matrice[col][lig] = (float) val;
        }
        }
    }
    
    public void setAbsence(int lig , int col) {
        if ( lig != col ) {
        if ( col < lig )
            presence[lig][col] = false;
        else // if ( lig < col )
            presence[col][lig] = false;
        }
    }   
    
    public boolean exists (int lig , int col) {
        if ( col != lig ) {
        if ( col < lig )
            return presence[lig][col];
        else // if ( lig < col )
            return presence[col][lig];
        }
        return true;
    }
    
    public int getSize () {
        return(taille);
    }
    
    public void setTaxon (String chaine , int lig) {
        if ( ! withTaxon ) { withTaxon = true; taxa = new String[taille]; }
        taxa[lig] = chaine;
    }
    
    public String getTaxon (int lig) {
        if ( withTaxon ) return taxa[lig];
        return "";
    }
    
    public void finalize() throws Throwable {
        matrice = null;
        taxa = null;
        super.finalize();
    }
    }


}

