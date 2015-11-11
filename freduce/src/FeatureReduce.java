/*
 * FeatureReduce.java - Todd Riley
 *
 */



import java.util.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;
import java.io.File;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.SwingConstants;
import java.awt.event.*;
import java.awt.*;
import java.text.NumberFormat;
import java.util.Hashtable;
import java.util.regex.Pattern;
import java.math.BigDecimal;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.net.URL;
import org.w3c.dom.*;
import org.xml.sax.*;

//import biojava.*;
//import biojava.BaumWelchTrainer;
//import biojava.TrainingAlgorithm;
import org.biojava.bio.*;
import org.biojava.bio.dist.*;
import org.biojava.bio.dp.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.db.*;
import org.biojava.bio.seq.io.*;
import org.biojava.bio.symbol.*;
import org.biojava.utils.*;
import org.biojava.bio.seq.impl.*;

import org.biojava.bio.mydp.*;

import org.rosuda.JRI.Rengine;
import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RList;
import org.rosuda.JRI.RVector;
import org.rosuda.JRI.RMainLoopCallbacks;

import org.apache.commons.lang.mutable.MutableDouble;
import org.apache.commons.lang.ArrayUtils;

public class FeatureReduce {

    ////////////////////////////////////////////////////////////////////////////////
    // Globals
    ////////////////////////////////////////////////////////////////////////////////

    private static String assigner = "=";
    private static String delimiter = ", ";
    private static String comma = ",";
    private static String commaWhiteSpace = "[,\\s]";
    private static String whiteSpace = "\\s+";
    private static String tabs = "\\t+";
    private static String tab = "\\t";
    private static Pattern assignerPattern = Pattern.compile(assigner);
    private static Pattern delimiterPattern = Pattern.compile(delimiter);
    private static Pattern commaPattern = Pattern.compile(comma);
    private static Pattern commaWhiteSpacePattern = Pattern.compile(commaWhiteSpace);
    private static Pattern whiteSpacePattern = Pattern.compile(whiteSpace);
    private static Pattern tabsPattern = Pattern.compile(tabs);
    private static Pattern tabPattern = Pattern.compile(tab);

    Output out = null;
    String outputLevelString = null;

    String initFile = null;
    KeyValueMap initFileCache = null;
    ReduceData reduceData = null;

    FiniteAlphabet alphabet = DNATools.getDNA();
    SymbolTokenization parser = null;
    ReversibleTranslationTable complementTable = DNATools.complementTable();

    Rengine rengine = null;
    WeightMatrix psam = null;
    double[][] initPositionalBias = null;
    int[] initStartPositions = null;
    double featureThresh = 0;

    //////////////////////////////////////////////////////
    // Globals for callbacks from R
    //////////////////////////////////////////////////////
    double[] globalY; // observed dependent vars
    double[][][] globalX; // [row][col] independent vars
    int globalNumWindows;
    boolean globalIsSelfRevComp;
    double globalLambda;
    double globalAlpha;
    double globalBeta;
    double globalNonSpecKa;
    double[] globalConcGammas;
    double[] globalLastParList;
    int globalNumStrands;

    //////////////////////////////////////////////////////
    // Globals for Poisson Regression
    //////////////////////////////////////////////////////
    double[] globalInitPoolFreqs = null;
    // double[] globalExpectedNoEnrichCounts = null;


    //////////////////////////////////////////////////////
    // Stuff that needs to be removed from GLOBALS
    //////////////////////////////////////////////////////
    ArrayList<FeaturedWeightMatrix> residualsFSAMs = new ArrayList<FeaturedWeightMatrix>();
    double[] residuals = null;
    boolean residualsSemaphore = false;
    double[] bgSubtractedIs = null;

    boolean[] removeProbesArray = null;
    int numRegressions = 1;

    double fsamUnivarCoeff = 0; // used to calculate the "offset(PSAM)"
    double fsamUnivarInter = 0;

    boolean psamOffset = false;
    LinkedHashMap<Object, FeaturedWeightMatrix> motifToFwmMap = null;

    WeightMatrix[] psamsArray = null; // used by -compareMotifs
    double[][][] psamsPosWeights = null; // used by -compareMotifs
    boolean addWeightsToLeft = true;
    boolean subtractOtherAffinities = false;

    String[] correspondenceStrings = null;

//     double intensitiesScaler = 1;
//     double intensitiesIntercept = 0;

    ////////////////////////////////////////////////////////////////////////////////
    // Constructors
    ////////////////////////////////////////////////////////////////////////////////

    public FeatureReduce() {
    }

    public FeatureReduce(String anInitFile) {
        this(
            anInitFile,
            null,
            null,
            null,
            null,
            null,
            null);
    }

    public FeatureReduce(
        String anInitFile,
        String seedFileString,
        String psamFileType,
        String psamPathName,
        String positionalBiasFileName) {
        this(
            anInitFile,
            seedFileString,
            psamFileType,
            psamPathName,
            positionalBiasFileName,
            null,
            null);
    }

    public FeatureReduce(
        String anInitFile,
        String seedFileString,
        String psamFileType,
        String psamPathName,
        String positionalBiasFileName,
        LinkedHashMap<Integer, LinkedHashMap<String, Double>> expNumToOrderedTrainIntensities,
        LinkedHashMap<Integer, LinkedHashMap<String, Double>> expNumToOrderedTestIntensities) {
        initialize(
            anInitFile,
            seedFileString,
            psamFileType,
            psamPathName,
            positionalBiasFileName,
            expNumToOrderedTrainIntensities,
            expNumToOrderedTestIntensities);
    }

    public void initialize(
        String anInitFile,
        String seedFileString,
        String psamFileType,
        String psamPathName,
        String positionalBiasFileName,
        LinkedHashMap<Integer, LinkedHashMap<String, Double>> expNumToOrderedTrainIntensities,
        LinkedHashMap<Integer, LinkedHashMap<String, Double>> expNumToOrderedTestIntensities) {

        try {
            // cache the anInitFile
            this.initFile = anInitFile;
	    if (initFileCache == null) {
		this.initFileCache = new KeyValueMap(anInitFile, "=");
	    }

            // initialize the output context
            this.outputLevelString = (String)initFileCache.get("OUTPUT", "Level");
            if (outputLevelString.equalsIgnoreCase("Debug")) {
                out = new Output(Output.Destination.STDOUT, (byte)0, false);
            }
            else if (outputLevelString.equalsIgnoreCase("Verbose")) {
                out = new Output(Output.Destination.STDOUT, (byte)25, false);
            }
            else {
                out = new Output(Output.Destination.STDOUT, (byte)100, false);
            }

            // create resultsDir
            String intensitiesDir = (String)initFileCache.get("DATA", "ProbeIntensitiesDir");
            String outputDir = (String)initFileCache.get("OUTPUT", "Directory");
            String resultsDir = intensitiesDir + File.separator + outputDir;
            FileTools.makeDir(resultsDir);

            // load initPositionalWeights
            if (positionalBiasFileName != null) {
                out.print("\nReading the positional-binding bias profile from "+positionalBiasFileName+"...");
                String[]  biasFileRows = FileTools.readStrings(positionalBiasFileName, 0);
                String[][] biasFileMatrix = StringTools.toStringMatrix(biasFileRows, null, true);
                this.initPositionalBias = StringTools.toDoubleMatrix(biasFileMatrix, true, true);
                out.println(" Done.\n");
                out.println("\nThe loaded positional-binding bias profile is..\n"+StringTools.toString(biasFileMatrix));
            }

            // init the reduceData context
            if (initFileCache.get("DREAM", "Label") == null) {
                reduceData = new ReduceData(initFileCache);
            }

            // create PSAM from seed seqs
            if (seedFileString != null) {
                // create the PSAM from the fasta file of aligned seqs
                out.print("\nCreating the PSAM from the same-length sequences found in the FASTA file "+seedFileString+"...");
                // Distribution[] dists = HMMTools.readDistsFromFastaFile(seedFileString, alphabet, "Yes", null, null, complementTable);
                Distribution[] dists = readDistsFromFastaFile(seedFileString, alphabet, "Yes", null, null, complementTable);
                psam = new SimpleWeightMatrix("seed PSAM", dists);
                out.println(" Done.\n");
            }
            // read PSAM from file
            else if (psamPathName != null) {
                // Load the PSAM from file
                out.print("\nReading the PSAM from "+psamFileType+" file "+psamPathName+"...");
                psam = loadPsam(psamFileType, psamPathName);
                out.println(" Done.\n");
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Methods
    ////////////////////////////////////////////////////////////////////////////////

    public WeightMatrix loadPsam(String format, String filePathName) {
        WeightMatrix pwm = null;
        if (format.equalsIgnoreCase("xml")) {
            pwm = WeightMatrixTools.readFromXML(filePathName);
        }
        else if (format.equalsIgnoreCase("table")) {
            pwm = WeightMatrixTools.readFromTable(filePathName);
        }
        else if (format.equalsIgnoreCase("mxr")) {
            pwm = WeightMatrixTools.readFromMatrixReduce(filePathName);
        }
        else if (format.equalsIgnoreCase("uniprobe")) {
            pwm = WeightMatrixTools.readFromUniprobe(filePathName);
        }
        else if (format.equalsIgnoreCase("jaspar")) {
            pwm = WeightMatrixTools.readFromJaspar(filePathName);
        }
        else {
            pwm = (WeightMatrix)FileTools.readSerializedFile(filePathName);
        }
        //pwm.setName(FileTools.stripPath(filePathName));
        pwm.setName("PSAM");
        return(pwm);
    }


    public static Distribution[] readDistsFromFastaFile(String fileName, Alphabet alphabet, String trainReverseCompString) {
        return(readDistsFromFastaFile(fileName, alphabet, trainReverseCompString, null, null, null));
    }
    public static Distribution[] readDistsFromFastaFile(String fileName,
        Alphabet alphabet,
        String trainReverseCompString,
        String[] aCorrStrings,
        String[] aFactorsStrings,
        ReversibleTranslationTable aComplementsTable) {
        Distribution[] dists = null;

        try {
            Map map = new HashMap();
            BufferedReader testBufferedReader = new BufferedReader(new FileReader(fileName));
            SymbolTokenization sTok = (SymbolTokenization) alphabet.getTokenization("token");
            SequenceIterator seqIter = SeqIOTools.readFasta(testBufferedReader, sTok);

            while (seqIter.hasNext()) {
                Sequence sequence = seqIter.nextSequence();
                map.put(sequence.getName(), sequence);
                if (trainReverseCompString.equalsIgnoreCase("Yes")) {
                    map.put(sequence.getName()+"ReverseComp", DNATools.reverseComplement(sequence));
                }
            }

            Alignment align = new SimpleAlignment(map);

            // Make a Distribution[] of the motif
            dists = DistributionTools.distOverAlignment(align);
            //dists = DistributionTools.distOverAlignment(align, false, .0000001);

            if (aCorrStrings != null) {
                DistributionCorrelations distCorrelations = new DistributionCorrelations(dists, aComplementsTable, aCorrStrings, aFactorsStrings);
                distCorrelations.correlate();
            }
        }

        catch (Exception ex) {
            ex.printStackTrace();
        }

        return(dists);
    }

    public void record(String aString) {
        out.println("\n"+aString);
    }

    private String[] getParams(String[] argsList, int index) {
        int counter = 0;
        while ((index+counter < argsList.length) && (!argsList[index+counter].startsWith("-"))) {
            counter++;
        }
        String[] params = new String[counter];
        for(int i=0; i < counter; i++) {
            params[i] = argsList[index+i];
        }
        //System.out.println("params_end="+params[params.length-1]);
        return(params);
    }

    public void parseCommandLine(String[] args) throws IOException {

        try {
            String title = "FeatureReduce";
            String version = "1.09";

            int numMandatoryArgs = 1;

            String genMatrixReduceInput = "No";
            String displayMotifsString = "Yes";
            boolean displayMotifs = true;
            boolean batch = false;

            String rnaBindingString = "Yes";
            boolean rnaBinding = true;

            boolean genWordCounts = false;

            boolean createFasta = false;
            //String motif = null;
            String motif = "1 2 3 4 5 6 7 8";
            String seed = null;
            SymbolList seedSymList = null;
            String[] referenceArray = null;
            String[] residualsArray = null;
            String[] residualsTypes = null;;
            double[] referenceKaArray = null;
            SymbolList[] refSymListArray = null;
            String seedFileString = null;
            String scoreType = null;
            String scoreArray[] = null;
            String scoreOutputFile = null;
            String psamFileType = null;
            String psamPathName = null;
            String savePosWeightsFile = null;
            String fsamPathName = null;
            String intensitiesPathName = null;
            String kmerLengthsArray[] = null;

            kmerLengthsArray = new String[7];
            kmerLengthsArray[0] = "4";
            kmerLengthsArray[1] = "5";
            kmerLengthsArray[2] = "6";
            kmerLengthsArray[3] = "7";
            kmerLengthsArray[4] = "8";
            kmerLengthsArray[5] = "9";
            kmerLengthsArray[6] = "10";

            String saveFormat = null;
            String saveFile = null;
            String savePsamFormat = null;
            String savePsamFile = null;

            //String saveLogoFormat = null;
            String saveLogoFile = null;

            String saveAllString = "Yes";
            boolean saveAll = true;

            String positionalBiasFileName = null;
            String compType = null;
            String compArray[] = null;
            String strandString = "both";
            String calcString = "NORMED_SUM";
            boolean equilibrate = false;

            String[] ids = null;
            String dreamLabel = null;
            String dreamInputFile = null;
            String dreamTemplateFile = null;
            String dreamOutputFile = null;
            int dreamLabelColumn = 0;
            int dreamSequenceColumn = 2;
            //int dreamIntensityColumn = 3; // signal mean
            int dreamIntensityColumn = 5; // singal median
            boolean dream = false;
            String kmerModelFlagString = "No";
            boolean kmerModelFlag = false;
            String noFeaturesFlagString = "No";
            boolean noFeaturesFlag = false;
            String noPosBiasFlagString = "No";
            boolean noPosBiasFlag = false;

            double affinitySphereThresh = 0;
            String affinitySphereFileName = null;

            String[] unifiedMotifs = null;
            String[] unifiedSeeds = null;
            String unifiedSeedFile = null;
            String seedType = null;

            String[] psamPathNames = null;
            String[] posWeightsPathNames = null;

            String proteinAlignmentFileString = null;
            String affinityModelListFileString = null;

            GridBagLayoutFrame psamLogo = null;
            GridBagLayoutFrame fsamLogo = null;
            GridBagLayoutFrame posBiasChart = null;

            System.out.println("\n"+title+" v"+version+" - Building biophysical protein-DNA and protein-RNA affinity models using robust regression routines in R and high-throughput binding assays.");
            System.out.println("\t\t    - written by Todd R. Riley\n");

            if ( (!(args.length==0 && numMandatoryArgs==0)) && (args.length < numMandatoryArgs || args[0].equalsIgnoreCase("-help") || args[0].equalsIgnoreCase("-?"))) {
                System.out.println("\n Usage: "+title+" <Init-File> -Option\n");
                System.out.println(" Options:");
                System.out.println();
                System.out.println("   -correspond [correspondence motif] ...");
                System.out.println();
                System.out.println("\t (Repeat Motif Example: -correspond \"1 2 3 4 1 2 3 4\")");
                System.out.println("\t (Palidromic Motif Example: -correspond \"1 2 3 4 4 3 2 1\")");
                System.out.println("\t (Reverse-Comp Palindromic Example: -correspond \"1 2 3 4 ~4 ~3 ~2 ~1\")");
                System.out.println("\t (Combined Example: -correspond \"1 2 3 ~3 ~2 ~1 1 2 3 ~3 ~2 ~1\")");
                System.out.println();
                System.out.println("\t If no motif is specified then all motifs of the given seed-length");
                System.out.println("\t in the motif library will be searched.");
                System.out.println();
                System.out.println("\n   -seed [iupac sequence] ...");
                System.out.println();
                System.out.println("\t Use -seed to specify the highest-affinity consensus to build a PSAM with.");
                System.out.println();
                System.out.println("\n   -seedFile [fasta file] ...");
                System.out.println();
                System.out.println("\t Use -seedFile to specify a list of  hight-affinity sites to build a PSAM with.");
                System.out.println("\t Note: All of the sites in the [fasta file] must be the same length.");
                System.out.println();
                System.out.println("\n   -seedType [all/palindrome/rcPalindrome] ...");
                System.out.println();
                System.out.println("\t Use -seedType to specify the type of seeds to test for highest correlation (to then build a PSAM with).");
                System.out.println();
                System.out.println("\n   -psam [xml/ser/table/mxr/uniprobe] File1 File2...");
                System.out.println();
                System.out.println("\t Use -psam to load in a starting PSAM for feature discovery.");
                System.out.println();
                System.out.println("\n   -fsam [FSAM1.ser] ...");
                System.out.println();
                System.out.println("\t Use -fsam to load an FSAM.");
                System.out.println();
                System.out.println("\n   -equilibrate");
                System.out.println();
                System.out.println("\t Use -equilibrate to re-calc the simpler models in an FSAM from the full model.");
                System.out.println();
                System.out.println("\n   -residuals [FSAM1.ser] ...");
                System.out.println();
                System.out.println("\t Use -residuals to load FSAMs and fit to their residuals.");
                System.out.println();
                System.out.println("\n   -positionalBias [file] ...");
                System.out.println();
                System.out.println("\t Use -positional bias to load a positional bias profile for binding to the probes.");
                System.out.println();
                System.out.println("\n   -displayMotifs <Yes/No> (default "+displayMotifsString+")");
                System.out.println();
                System.out.println("\t Use -displayMotifs to turn on/off the display of PSAM and FSAM logos.");
                System.out.println();
                System.out.println("\n   -score <seq/table/fasta> <*.fasta> ...");
                System.out.println();
                System.out.println("\n   -scoreOutputFile <File-Name> ...");
                System.out.println();
                System.out.println("\n   -strand <pos/neg/both> (default "+strandString+")");
                System.out.println();
                System.out.println("\n   -bothStrandsCalc <max/sum/union/normed_sum> (default "+calcString+")");
                System.out.println();
                System.out.println("\n   -rnaBinding <Yes/No> (default "+rnaBindingString+")");
                System.out.println();
                System.out.println("\n   -scoreOutputFile <File-Name> ...");
                System.out.println();
//                 System.out.println("\n   -dream <Input-File> <Template-File> <Output-File>");
//                 System.out.println();
                System.out.println("\n   -affinitySphere <FeatureThresh> <TotalAffinityThresh> <File-Name>");
                System.out.println();
                System.out.println("\t Saves the relative affinity sphere all sequences with relative affinity >= threshold.");
                System.out.println();
                System.out.println("\n   -reference <seq>");
                System.out.println();
                System.out.println("\n   -referenceKa <Ka>");
                System.out.println();
                System.out.println("\n   -comp table <file>");
                System.out.println();
                System.out.println("   -save <ser/xml> <File-Name>");
                System.out.println();
                System.out.println("   -savePsam <ser/xml/uniprobe> <File-Name>");
                System.out.println();
                System.out.println("   -savePosWeights <File-Name>");
                System.out.println();
                System.out.println("   -saveLogo <File-Name>");
                System.out.println();
                System.out.println("\n   -topN <count>");
                System.out.println();
                System.out.println("\t Use -topN to specify the number of top N intensities/counts for regression.");
                System.out.println();
                System.out.println("   -batch");
                System.out.println();
                System.out.println("\t Use -batch to force termination after execution and plot generation.");
                System.out.println();
                System.out.println("\n   -genWordCounts <length> <length> ...");
                System.out.println();
                System.out.println("\t Use -genWordCounts to generate the Word-Counts File for the");
                System.out.println("\t specified kmer-lengths. (Example: \"-genWordCounts 5 6 7\"");
                System.out.println("\t will generate 3 Word-Counts files for all 5mers, 6mers, and");
                System.out.println("\t 7mers, respectively.)");
                System.out.println();
                System.out.println("\n   -createFasta");
                System.out.println();
                System.out.println("\n   -genIntensities <outputPathName>");
                System.out.println();
                System.out.println("\t Use -genIntensities to generate an intensities file from a PSAM");
                System.out.println("\t and the sequence data.");
                System.out.println();
                System.out.println("\n   -clearIntensities");
                System.out.println();
                System.out.println("\t Use -clearIntensities if you have changed the probe intensity files or");
                System.out.println("\t normalization parameters. This option will force re-reads and");
                System.out.println("\t re-normalization of the probe-intensity data");
                System.out.println();
                System.out.println("\n   -clearCounts");
                System.out.println();
                System.out.println("\t Use -clearCounts if you have changed the probe-associated sequences");
                System.out.println("\t in any way. This option will force a re-read of the sequence data,");
                System.out.println("\t and re-calculation of all the occurrences of motifs that are");
                System.out.println("\t pertubations from the seed motif(s).");
                System.out.println();
                System.out.println("\n   -genInput <Yes/No> (default "+genMatrixReduceInput+")");
                System.out.println("\n   -help, -? \tDisplays this help message");
                System.out.println("\n   -version, \tDisplays version number");
                System.out.println();
                System.out.println("See reduce.init for the format of the <Init-File>.");
                System.out.println("");

                System.out.flush();

                if ((args.length > 0) && (args[0].equalsIgnoreCase("-help") || args[0].equalsIgnoreCase("-?"))) {
                    System.exit(0);
                }
                else {
                    System.exit(-1);
                }
            }

            String anInitFile=args[0];

            for (int a=0; a < args.length; a++) {
                if (args[a].equals("-version")) {
                    System.out.println("\n"+title+" v"+version+"\n");
                    System.exit(0);
                }
                else if (args[a].equals("-motif") || args[a].equals("-correspond")) { // "-motif" is for backawrd compatibility
                    motif = new String(args[a+1]);
                    this.correspondenceStrings = whiteSpacePattern.split(motif);
                    //motif = getParam(args, a+1);
                }
                else if (args[a].equals("-seed")) {
                    seed = new String(args[a+1]);
                    seedSymList = DNATools.createDNA(seed);
                    //seedsArray = getParams(args, a+1);
                }
                else if (args[a].equals("-seedFile")) {
                    seedFileString = new String(args[a+1]);
                    //seedsArray = getParams(args, a+1);
                }
                else if (args[a].equals("-seedType")) {
                    seedType = new String(args[a+1]);
                    //seedsArray = getParams(args, a+1);
                }
                else if (args[a].equals("-psam")) {
                    psamFileType = new String(args[a+1]);
                    psamPathName = new String(args[a+2]);
                    //psamsArray = getParams(args, a+1);
                }
                else if (args[a].equals("-fsam")) {
                    fsamPathName = new String(args[a+1]);
                    //fsamsArray = getParams(args, a+1);
                }
                else if (args[a].equals("-residuals")) {
                    residualsTypes = getParams(args, a+1);
                    residualsArray = getParams(args, a+2);

                    for (int i=0; i < residualsArray.length; i++) {

                        String aPsamFileType = residualsTypes[i];
                        String aPsamPathName = residualsArray[i];

                        System.out.print("\nReading the PSAM from "+aPsamFileType+" file "+aPsamPathName+"...");
                        WeightMatrix psam = loadPsam(aPsamFileType, aPsamPathName);
                        System.out.println(" Done.\n");

                        String name = ((SimpleWeightMatrix) psam).getName();
                        if ((name == null) || (name.startsWith("-"))) {
                            name = "FSAM";
                        }

                        FeaturedWeightMatrix anFsam = new FeaturedWeightMatrix(
                            //"FSAM",
                            //anInitFile,
                            name,
                            psam,
                            ReduceData.getStrand(strandString),
                            ReduceData.getBothStrandsCalc(calcString),
                            //                     WeightMatrixTools.BindingStrand.BOTH,
                            //                     WeightMatrixTools.BothStrandsCalc.NORMED_SUM,
                            complementTable);


                        // FeaturedWeightMatrix anFsam = (FeaturedWeightMatrix)FileTools.readSerializedFile(residualsArray[i]);

                        alphabet = (FiniteAlphabet)anFsam.getAlphabet();

                        //anFsam.reHashFeatures();
                        anFsam.removeFeatures();

                        residualsFSAMs.add(anFsam);
                    }
                }
                else if (args[a].equals("-positionalBias")) {
                    positionalBiasFileName = new String(args[a+1]);
                    //psamsArray = getParams(args, a+1);
                }
                else if (args[a].equals("-genWordCounts")) {
                    kmerLengthsArray = getParams(args, a+1);
                    genWordCounts = true;
                }
                else if (args[a].equals("-compareMotifs")) {
                    psamPathNames = getParams(args, a+1);
                }
                else if (args[a].equals("-comparePosWeights")) {
                    posWeightsPathNames = getParams(args, a+1);
                }
                else if (args[a].equals("-genIntensities")) {
                    intensitiesPathName = new String(args[a+2]);
                }
                else if (args[a].equals("-genInput")) {
                    genMatrixReduceInput = new String(args[a+1]);
                }
                else if (args[a].equals("-score")) {
                    scoreType = new String(args[a+1]);
                    scoreArray = getParams(args, a+2);
                }
                else if (args[a].equals("-kmer")) {
//                     kmerModelFlagString = new String(args[a+1]);
//                     kmerModelFlag = StringTools.parseBoolean(kmerModelFlagString);
                    kmerModelFlagString = new String("Yes");
                    kmerModelFlag = true;
                }
                else if (args[a].equals("-noFeatures")) {
                    noFeaturesFlagString = new String("Yes");
                    noFeaturesFlag = true;
                }
                else if (args[a].equals("-noPosBias")) {
                    noPosBiasFlagString = new String("Yes");
                    noPosBiasFlag = true;
                }


		////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////

                else if (args[a].equals("-s")) {
                    dream = true;
                    if (this.initFileCache == null) {
                        this.initFileCache = new KeyValueMap(anInitFile, "=");
                    }
                    int size = Integer.parseInt(args[a+1]);

                    switch (size) {
                    case 4:
                        motif = "1 2 3 4";
                        break;
                    case 5:
                        motif = "1 2 3 4 5";
                        break;
                    case 6:
                        motif = "1 2 3 4 5 6";
                        break;
                    case 7:
                        motif = "1 2 3 4 5 6 7";
                        break;
                    case 8: // default
                        motif = "1 2 3 4 5 6 7 8";
                        break;
                    case 9:
                        motif = "1 2 3 4 5 6 7";
                        initFileCache.put("AFFINITY MODELS", "SeedMotifMinExtensionLength", "1");
                        initFileCache.put("AFFINITY MODELS", "SeedMotifMaxExtensionLength", "1");
                        break;
                    case 10:
                        motif = "1 2 3 4 5 6 7 8";
                        initFileCache.put("AFFINITY MODELS", "SeedMotifMinExtensionLength", "1");
                        initFileCache.put("AFFINITY MODELS", "SeedMotifMaxExtensionLength", "1");
                        break;
                    case 11:
                        motif = "1 2 3 4 5 6 7";
                        initFileCache.put("AFFINITY MODELS", "SeedMotifMinExtensionLength", "1");
                        initFileCache.put("AFFINITY MODELS", "SeedMotifMaxExtensionLength", "2");
                        break;
                    case 12:
                        motif = "1 2 3 4 5 6 7 8";
                        initFileCache.put("AFFINITY MODELS", "SeedMotifMinExtensionLength", "1");
                        initFileCache.put("AFFINITY MODELS", "SeedMotifMaxExtensionLength", "2");
                        break;

                    }
                }
                else if (args[a].equals("-l")) {
                    dream = true;
                    if (this.initFileCache == null) {
                        this.initFileCache = new KeyValueMap(anInitFile, "=");
                    }
                    dreamLabel = new String(args[a+1]);
                }
                else if (args[a].equals("-c")) {
                    dream = true;
                    if (this.initFileCache == null) {
                        this.initFileCache = new KeyValueMap(anInitFile, "=");
                    }
                    dreamLabelColumn = Integer.parseInt(args[a+1]);
                    dreamSequenceColumn = Integer.parseInt(args[a+2]);
                    dreamIntensityColumn = Integer.parseInt(args[a+3]);
                }
                else if (args[a].equals("-ids")) {
                    ids = getParams(args, a+1);
                }
                else if (args[a].equals("-i")) {
                    dream = true;
                    if (this.initFileCache == null) {
                        this.initFileCache = new KeyValueMap(anInitFile, "=");
                    }
                    dreamInputFile = new String(args[a+1]);
                    initFileCache.put("SEQUENCE", "1", dreamInputFile);
                    initFileCache.put("INTENSITIES", "1", dreamInputFile);

                    batch = true;

                    // initial motif length
                    //motif = "1 2 3 4 5 6 7 8";
                    //motif = "1 2 3 4 5 6 7";
                    //motif = "1 2 3 4 5 6";

//                     kmerLengthsArray = new String[5];
//                     kmerLengthsArray[0] = "4";
//                     kmerLengthsArray[1] = "5";
//                     kmerLengthsArray[2] = "6";
//                     kmerLengthsArray[3] = "7";
//                     kmerLengthsArray[4] = "8";

//                     kmerLengthsArray = new String[3];
//                     kmerLengthsArray[0] = "4";
//                     kmerLengthsArray[1] = "6";
//                     kmerLengthsArray[2] = "8";


//                     kmerLengthsArray = new String[2];
//                     kmerLengthsArray[0] = "6";
//                     kmerLengthsArray[1] = "8";


                    // motif = "1 2 3 4 5 6 7";

                    // kmerLengthsArray = new String[5];
                    // kmerLengthsArray[0] = "4";
                    // kmerLengthsArray[1] = "5";
                    // kmerLengthsArray[2] = "6";
                    // kmerLengthsArray[3] = "7";
                    // kmerLengthsArray[4] = "8";

                    // kmerLengthsArray = new String[4];
                    // kmerLengthsArray[0] = "8";
                    // kmerLengthsArray[1] = "9";
                    // kmerLengthsArray[2] = "10";
                    // kmerLengthsArray[3] = "11";

                    // kmerLengthsArray = new String[1];
                    // kmerLengthsArray[0] = "8";

                }
                else if (args[a].equals("-a")) {
                    dream = true;
                    if (this.initFileCache == null) {
                        this.initFileCache = new KeyValueMap(anInitFile, "=");
                    }
                    dreamTemplateFile = new String(args[a+1]);
                    initFileCache.put("DREAM", "TemplateFile", dreamTemplateFile);
                }
                else if (args[a].equals("-o")) {
                    dream = true;
                    if (this.initFileCache == null) {
                        this.initFileCache = new KeyValueMap(anInitFile, "=");
                    }
                    dreamOutputFile = new String(args[a+1]);
                    initFileCache.put("DREAM", "OutputFile", dreamOutputFile);
                }

		////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////

                else if (args[a].equals("-comp")) {
                    compType = new String(args[a+1]);
                    compArray = getParams(args, a+2);
                }
                else if (args[a].equals("-scoreOutputFile")) {
                    scoreOutputFile = new String(args[a+1]);
                }
                else if (args[a].equals("-topN")) {
                    if (this.initFileCache == null) {
                        this.initFileCache = new KeyValueMap(anInitFile, "=");
                    }
                    String topN = new String(args[a+1]);
                    initFileCache.put("REGRESSION", "NumberOfProbes", topN);
                }
                else if (args[a].equals("-strand")) {
                    strandString = new String(args[a+1]);
                }
                else if (args[a].equals("-bothStrandsCalc")) {
                    calcString = new String(args[a+1]);
                }
                else if (args[a].equals("-reference")) {
                    //reference = new String(args[a+1]);
                    referenceArray = getParams(args, a+1);
                    refSymListArray = new SymbolList[referenceArray.length];
                    for (int i=0; i < referenceArray.length; i++) {
                        refSymListArray[i] = DNATools.createDNA(referenceArray[i]);
                    }
                }
                else if (args[a].equals("-referenceKa")) {
                    String[] referenceKaStrings = getParams(args, a+1);
                    referenceKaArray = new double[referenceKaStrings.length];
                    for (int i=0; i < referenceKaStrings.length; i++) {
                        referenceKaArray[i] = Double.parseDouble(args[a+1]);
                    }
                }
                else if (args[a].equals("-createFasta")) {
                    createFasta = true;
                }
                else if (args[a].equals("-displayMotifs")) {
                    displayMotifsString = new String(args[a+1]);
                    displayMotifs = StringTools.parseBoolean(displayMotifsString);
                }
                else if (args[a].equals("-rnaBinding")) {
                    rnaBindingString = new String(args[a+1]);
                    rnaBinding = StringTools.parseBoolean(rnaBindingString);
                }
                else if  (args[a].equalsIgnoreCase("-save")) {
                    saveFormat = new String(args[a+1]);
                    saveFile = new String(args[a+2]);
                }
                else if  (args[a].equalsIgnoreCase("-savePsam")) {
                    savePsamFormat = new String(args[a+1]);
                    savePsamFile = new String(args[a+2]);
                }
                else if  (args[a].equalsIgnoreCase("-savePosWeights")) {
                    savePosWeightsFile = new String(args[a+1]);
                }
                else if  (args[a].equalsIgnoreCase("-saveLogo")) {
                    // saveLogoFormat = new String(args[a+1]);
                    // saveLogoFile = new String(args[a+2]);

                    saveLogoFile = new String(args[a+1]);
                }
                else if (args[a].equals("-batch")) {
                    batch = true;
                }
                else if (args[a].equals("-equilibrate") || args[a].equals("-e")) {
                    equilibrate = true;
                }
                else if (args[a].equals("-unified")) {
                    unifiedMotifs = getParams(args, a+1);
                }
                else if (args[a].equals("-unifiedSeeds")) {
                    unifiedSeeds = getParams(args, a+1);
                }
                else if (args[a].equals("-unifiedSeedFile")) {
                    unifiedSeedFile = new String(args[a+1]);
                }
                else if (args[a].equals("-saveAll")) {
                    saveAllString = new String(args[a+1]);
                    saveAll = StringTools.parseBoolean(saveAllString);
                }
                else if (args[a].equals("-family")) {
                    proteinAlignmentFileString = new String(args[a+1]);
                    affinityModelListFileString = new String(args[a+2]);
                }
                else if (args[a].equalsIgnoreCase("-affinitySphere")) {
                    featureThresh = Double.parseDouble(args[a+1]);
                    affinitySphereThresh = Double.parseDouble(args[a+2]);
                    affinitySphereFileName = new String(args[a+3]);
                }
                else if (args[a].startsWith("-")) {
                    System.out.println("\nError: Unknown option \""+args[a]+"\"");
                    System.exit(-1);
                }
            }

            if (dream == true) {
                if (dreamLabel == null) {
                    dreamLabel = dreamInputFile;
                }

                initFileCache.put("DREAM", "Label", dreamLabel);
                initFileCache.put("DREAM", "LabelColumn", dreamLabelColumn);
                initFileCache.put("DREAM", "SequenceColumn", dreamSequenceColumn);
                initFileCache.put("DREAM", "IntensityColumn", dreamIntensityColumn);
            }

            if (kmerModelFlag == true) {

                // kmerLengthsArray = new String[5];
                // kmerLengthsArray[0] = "4";
                // kmerLengthsArray[1] = "5";
                // kmerLengthsArray[2] = "6";
                // kmerLengthsArray[3] = "7";
                // kmerLengthsArray[4] = "8";
            }

            if (!anInitFile.startsWith("-")) {
                initialize(anInitFile, seedFileString, psamFileType, psamPathName, positionalBiasFileName, null, null);
            }

            if (createFasta) {
                //out.print("\nCreating Fasta Files...");
                reduceData.createFastaFile("/users/triley/microarray/human/manley/manley.mat.bar.sorted.txt", ".", "tls.mat.sorted.hg18", 600, 25, "Center", "-", true);
                //out.println(" Done.\n");
                System.exit(0);
            }

            FeaturedWeightMatrix fsam = null;

            if (psamPathNames != null) {

                /////////////////////////////
                // Load the PSAMs
                /////////////////////////////
                ArrayList<WeightMatrix> psamsArrayList = new ArrayList<WeightMatrix>();

                for (String aPsamPathName : psamPathNames) {

                    String aPsamFileType = "xml";

                    // Load the PSAM from file
                    System.out.print("\nReading the PSAM from "+aPsamFileType+" file "+aPsamPathName+"...");
                    WeightMatrix aPsam = loadPsam(aPsamFileType, aPsamPathName);
                    psamsArrayList.add(aPsam);
                    System.out.println(" Done.\n");

                    if (displayMotifs) {
                        psamLogo = display("PSAM Logo: "+aPsamPathName, aPsam);
                    }
                }

                psamsArray = (WeightMatrix[])psamsArrayList.toArray(new WeightMatrix[0]);


                /////////////////////////////
                // Load the Pos Weights
                /////////////////////////////

                if (posWeightsPathNames != null) {

                    ArrayList<double[][]> psamsPosWeightsArrayList = new ArrayList<double[][]>();

                    for (String posWeightsPathName : posWeightsPathNames) {

                        // Load the PSAM from file
                        System.out.print("\nReading the positional-binding bias profile from "+posWeightsPathName+"...");
                        String[]  biasFileRows = FileTools.readStrings(posWeightsPathName, 0);
                        String[][] biasFileMatrix = StringTools.toStringMatrix(biasFileRows, null, true);
                        double[][] positionalBias = StringTools.toDoubleMatrix(biasFileMatrix, true, true);
                        psamsPosWeightsArrayList.add(positionalBias);
                        System.out.println(" Done.\n");
                        System.out.println("\nThe loaded positional-binding bias profile is..\n"+StringTools.toString(biasFileMatrix));
                        System.out.println("\nThe positional weights are "+positionalBias[0].length+"bp long.");
                    }

                    psamsPosWeights = (double[][][])psamsPosWeightsArrayList.toArray(new double[0][0][0]);
                }
                else {
                    psamsPosWeights = new double[psamPathNames.length][][];
                }

                fitPsams(displayMotifs, kmerLengthsArray, -1);
            }

            // if don't do more than just generate Input Files
            else if (genMatrixReduceInput.equalsIgnoreCase("Yes")) {
            // if (genMatrixReduceInput.equalsIgnoreCase("Yes")) {
                return;
            }
            else if (genWordCounts) {
                //else if ((kmerLengthsArray != null) && (dreamLabel == null)) {
                genWordCounts(kmerLengthsArray);
                return;
            }
            else if (intensitiesPathName != null) {
                genIntensities(intensitiesPathName);
                return;
            }
            else if (proteinAlignmentFileString != null) {
                fitModel(proteinAlignmentFileString, affinityModelListFileString);
            }
            else if (fsamPathName != null) {
                fsam = (FeaturedWeightMatrix)FileTools.readSerializedFile(fsamPathName);
                alphabet = (FiniteAlphabet)fsam.getAlphabet();
                fsam.reHashFeatures();
                //System.out.println("\nHi!\n");

                if (equilibrate) {
                    fsam.equilibrate();
                    //fsam.equilibrate();
                }

                // Display the FSAM if desired
                if (displayMotifs) {
                    fsamLogo = display("FSAM Logo: "+fsamPathName, fsam, featureThresh);

                    psamLogo = display("Old PSAM Logo: "+fsamPathName, fsam.oldAMposStrandPWM);

                    double[][] fsamPositionalWeights = fsam.getPositionalWeights();
                    if (fsamPositionalWeights != null) {
                        int probeSeqLengths = fsam.positionalWeights[0].length - fsam.columns() + 1;
                        int[] startPositions = getStartPositions(probeSeqLengths, fsam.columns);

                        String[] datasetLabels = {"Positive Strand", "Negative Strand"};
                        double[][] values = {ArrayTools.toDoubleArray(startPositions), fsamPositionalWeights[0], fsamPositionalWeights[1]};
                        BioJavaChart chart = new BioJavaChart("lines", null, "Start Position", "Normalized Coefficients", null, datasetLabels, true, values, -1);
                        posBiasChart = display("Positional Bias Profile inside "+fsamPathName, chart);
                    }

                }
            }
            else if (dreamInputFile != null) {
                fsam = fitModels(motif,
				 unifiedSeeds,
				 psam,
				 displayMotifs,
				 null,
				 Double.NaN,
				 ids,
				 seedType,
				 saveAll,
				 kmerLengthsArray,
				 dreamInputFile);

            }
            else if (unifiedMotifs != null) {
                //fsam = fitModels(motif, unifiedSeeds, psam, displayMotifs, refSymListArray[0], referenceKaArray[0], unifiedMotifs, saveAll, seedType);
                fsam = fitModels(motif, unifiedSeeds, psam, displayMotifs, null, Double.NaN, unifiedMotifs, saveAll, seedType);
            }
            else if (unifiedSeedFile != null) {
                String[] seedFileRows = FileTools.readStrings(unifiedSeedFile, 0);
                String[][] seedFileMatrix = StringTools.toStringMatrix(seedFileRows, null, false);

                //fsam = fitModels(motif, seedFileMatrix[1], psam, displayMotifs, refSymListArray[0], referenceKaArray[0], seedFileMatrix[0], saveAll, seedType);
                fsam = fitModels(motif, seedFileMatrix[1], psam, displayMotifs, null, Double.NaN, seedFileMatrix[0], saveAll, seedType);
            }
            else if (!anInitFile.startsWith("-") && (dream == false)) {
                //fsam = fitModels(motif, seedSymList, psam, displayMotifs, refSymListArray[0], referenceKaArray[0], FileTools.stripExtension(anInitFile), -1, saveAll, seedType);
                fsam = fitModels(motif, seedSymList, psam, displayMotifs, null, Double.NaN, FileTools.stripExtension(anInitFile), -1, saveAll, seedType);
            }
            else if (psamPathName != null) {
                // Load the PSAM from file
                System.out.print("\nReading the PSAM from "+psamFileType+" file "+psamPathName+"...");
                psam = loadPsam(psamFileType, psamPathName);
                System.out.println(" Done.\n");

                String name = anInitFile;
                if ((name == null) || (name.startsWith("-"))) {
                    name = "FSAM";
                }

                fsam = new FeaturedWeightMatrix(
                    //"FSAM",
                    //anInitFile,
                    name,
                    psam,
                    ReduceData.getStrand(strandString),
                    ReduceData.getBothStrandsCalc(calcString),
//                     WeightMatrixTools.BindingStrand.BOTH,
//                     WeightMatrixTools.BothStrandsCalc.NORMED_SUM,
                    complementTable);

                // load initPositionalWeights
                if (positionalBiasFileName != null) {
                    System.out.print("\nReading the positional-binding bias profile from "+positionalBiasFileName+"...");
                    String[]  biasFileRows = FileTools.readStrings(positionalBiasFileName, 0);
                    String[][] biasFileMatrix = StringTools.toStringMatrix(biasFileRows, null, true);
                    this.initPositionalBias = StringTools.toDoubleMatrix(biasFileMatrix, true, true);
                    this.initStartPositions = StringTools.toInts(StringTools.columnLabels(biasFileMatrix));
                    fsam.setPositionalWeights(this.initPositionalBias);
                    System.out.println(" Done.\n");
                    System.out.println("\nThe loaded positional-binding bias profile is..\n"+StringTools.toString(biasFileMatrix));
                    System.out.println("\nThe positional weights are "+this.initPositionalBias[0].length+"bp long.");
                }

                if (displayMotifs) {
                    psamLogo = display("PSAM Logo: "+psamPathName, psam);

                    if (this.initPositionalBias != null) {
                        // int probeSeqLengths = 35;
                        // int[] startPositions = getStartPositions(probeSeqLengths, fsam.columns);

                        String[] datasetLabels = {"Positive Strand", "Negative Strand"};
                        double[][] values = {ArrayTools.toDoubleArray(this.initStartPositions), this.initPositionalBias[0], this.initPositionalBias[1]};
                        BioJavaChart chart = new BioJavaChart("lines", null, "Start Position", "Normalized Coefficients", null, datasetLabels, true, values, -1);
                        posBiasChart = display("Positional Bias Profile: "+positionalBiasFileName, chart);
                    }
                }

            }


            // Display the FSAM if desired
//             if (displayMotifs) {
//                 GridBagLayoutFrame fsamLogo = display("FSAM Logo", fsam, featureThresh);
//                 GridBagLayoutFrame psamLogo = display("PSAM Logo", fsam.getPosStrandPWM());

// //                 if (saveAll) {
// //                     String intensitiesDir = (String)initFileCache.get("DATA", "ProbeIntensitiesDir");
// //                     String outputDir = (String)initFileCache.get("OUTPUT", "Directory");
// //                     String resultsDir = intensitiesDir + File.separator + outputDir;
// //                     GraphicsTools.writeImage(fsamLogo.getPanel(), resultsDir + File.separator + proteinLabel+".psam.png");
// //                 }
//             }

            //fsam.setRevCompSimilarity(1.0);
            //fsam.setNonSpecKa(0.0);

            ////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////

            if (fsam != null) {
                System.out.println("\nFSAM:\n"+fsam.toString(false));
                System.out.println("\nFSAM:\n"+fsam.toString(true));
            }

            //fsam.getHighestScore(featureThresh);

            // Save the FSAM
            if (saveFile != null) {
                save(fsam, saveFormat, saveFile);
            }

            // Save the PSAM
            if (savePsamFile != null) {
                save(fsam.getPosStrandPWM(), savePsamFormat, savePsamFile);
            }

            // Save the PosWeights
            if (savePosWeightsFile != null) {
                if (fsam.getPositionalWeights() != null) {
                    //int probeSeqLengths = reduceData.getSeqLength();
                    int probeSeqLengths = fsam.positionalWeights[0].length - fsam.columns() + 1;
                    int[] startPositions = getStartPositions(probeSeqLengths, fsam.columns);
                    save(fsam.getPositionalWeights(), startPositions, savePosWeightsFile);
                }
                else {
                    System.out.println("\nCannot Save File: Positional Weights are null!!");
                }
            }

            // Save the Logo
            if (saveLogoFile != null) {
                //save(psamLogo, fsamLogo, saveLogoFormat, saveLogoFile);
                save(psamLogo, fsamLogo, saveLogoFile);
            }


//             if (scoreArray != null) {
//                 // skip the first line
//                 scoreArray = FileTools.readStrings(scoreArray[0], 1);
//             }

            // score any sequences
            if (scoreArray != null) {
                //fsam.setStrand(WeightMatrixTools.BindingStrand.BOTH);
                //fsam.setCalc(WeightMatrixTools.BothStrandsCalc.SUM);

                WeightMatrixTools.BindingStrand scoreStrand = fsam.getStrand();
                if (scoreStrand != WeightMatrixTools.BindingStrand.BOTH) {
                    System.out.println("\nScoring will be performed using the "+scoreStrand+" strand.");
                }
                else {
                    System.out.println("\nScoring will be performed using "+scoreStrand+" strands with the "+fsam.getCalc()+" calculation.");
                }
                String scoreString = score(fsam, scoreType, scoreArray, refSymListArray, referenceKaArray, compType, compArray);

                if (scoreOutputFile != null) {
                    FileTools.write(scoreString, scoreOutputFile, false);
                }
                else {
                    System.out.println("\n"+scoreString);
                }
            }
            else if (dreamTemplateFile != null) {
                // score the sequences in the template file
                score(ids, dreamLabel, dreamTemplateFile, dreamOutputFile, kmerModelFlag, !noFeaturesFlag, !noPosBiasFlag);
            }

            if (affinitySphereFileName != null) {
                writeRelAffinitySphere(fsam, featureThresh, affinitySphereThresh, affinitySphereFileName);
            }

            // if batch mode then terminate
            if (batch) {
                System.exit(0);
            }
        }

        catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public void save(
        GridBagLayoutFrame psamLogoFrame,
        GridBagLayoutFrame fsamLogoFrame,
        String saveLogoFile) {

        if (psamLogoFrame != null) {
            GraphicsTools.writeImage(psamLogoFrame.getPanel(), saveLogoFile, false);
        }
        else {
            GraphicsTools.writeImage(fsamLogoFrame.getPanel(), saveLogoFile, false);
        }

    }

    public void save(
        double[][] posWeights,
        int[] startPositions,
        String filePathName)
    {
        String[] rowLabels = {"posStrandRelBias", "negStrandRelBias"};
        String[] columnLabels = StringTools.toStrings(startPositions);
        String outputTable = StringTools.toString(posWeights, rowLabels, columnLabels);
        FileTools.write(
            outputTable,
            filePathName,
            false);
    }


    public GridBagLayoutFrame display(String title, FeaturedWeightMatrix fsam, double featureThresh) {
        GridBagLayoutFrame logoFrame = new GridBagLayoutFrame(title, true);

        // default is get Mult Features
        ArrayList<WeightMatrixFeature> features = fsam.getFeatures(featureThresh);

        if (features != null) {
            // Sort the features first
            Table featuresTable = new Table();
            for (int index=0; index < features.size(); index++) {
                WeightMatrixFeature feature = features.get(index);
                double multKa = feature.getRelativeAffinity();
                double dddG = WeightMatrixTools.getGibbsFreeEnergy(multKa);

                java.util.List keyValueList = Arrays.asList(
                    feature,
                    new Double(multKa),
                    new Double(dddG),
                    new Double(Math.abs(dddG)));

                featuresTable.add(keyValueList);
            }

            featuresTable.sort(3); //sort by ascending abs(dddG)
            featuresTable.reverse(); // get descending order

            // this will not change features in the calling function
            features = (ArrayList)featuresTable.getColumn(0);
        }

        WeightMatrixLogo logoPanel = new WeightMatrixLogo(
            fsam.getPosStrandPWM(),
            true,
            true,
            true,
            0,
            WeightMatrixLogo.PositionLoc.TOP,
            features,
            false,
            .0001,
            0);

        logoFrame.add(logoPanel, 10, 10, false, true);
        return(logoFrame);
    }

    public GridBagLayoutFrame display(String title, WeightMatrix psam) {
        GridBagLayoutFrame logoFrame = new GridBagLayoutFrame(title, true);
        WeightMatrixLogo logoPanel = new WeightMatrixLogo(psam, true, true, true, 0, WeightMatrixLogo.PositionLoc.BOTTOM, null, false, .0001, 0);
        logoFrame.add(logoPanel, 10, 10, false, true);
        return(logoFrame);
    }

    public void add(GridBagLayoutFrame logoFrame, WeightMatrix psam) {
        WeightMatrixLogo logoPanel = new WeightMatrixLogo(psam, true, true, true, 0, WeightMatrixLogo.PositionLoc.BOTTOM, null, false, .0001, 0);
        logoFrame.add(logoPanel, 10, 10, false, true);
    }

    public GridBagLayoutFrame display(String frameTitle, BioJavaChart chart) {
        GridBagLayoutFrame logoFrame = new GridBagLayoutFrame(frameTitle, true);
        logoFrame.add(chart, 10, 10, false, true);
        return(logoFrame);
    }

    public GridBagLayoutFrame display(
        String frameTitle,
        String chartType,
        String chartTitle,
        String xAxisLabel,
        String yAxisLabel,
        String zAxisLabel,
        String[] datasetLabels,
        boolean legend,
        double[][] values,
        int bins) {

        //double[][] transformedValues = MathTools.transform(values);

        BioJavaChart chart = new BioJavaChart(
            chartType,
            chartTitle,
            xAxisLabel,
            yAxisLabel,
            zAxisLabel,
            datasetLabels,
            legend,
            //transformedValues,
            values,
            bins);

        return(display(frameTitle, chart));
    }

    public double getKaRef(
        Object motif,
        SymbolList refSymList,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc) {

        double refKaFitted = 1.0;
        if (refSymList != null) {
            if (motif instanceof WeightMatrix) {
                refKaFitted = WeightMatrixTools.score(
                    (WeightMatrix) motif,
                    refSymList,
                    strand,
                    calc,
                    complementTable,
                    0,
                    0);
            }
            else if (motif instanceof FeaturedWeightMatrix) {
                refKaFitted = ((FeaturedWeightMatrix)motif).score(refSymList, (FeatureKey)null);
            }
        }
        return(refKaFitted);
    }

    public double[] normalizeKaArray(double[] kAarray, double refKaFitted, double refKaReal) {
        double scaler = refKaReal / refKaFitted;
        double[] normalizedKaArray = MathTools.multiply(scaler, kAarray);
        return(normalizedKaArray);
    }

    public void save(FeaturedWeightMatrix fsam, String format, String filePathName) {
        if(format.equalsIgnoreCase("xml")) {
            // FIX ME!!
            //WeightMatrixTools.writeToXML(fsam, filePathName);
        }
        else {
            FileTools.writeSerializedFile(fsam, filePathName);
        }
    }

    public void save(WeightMatrix psam, String format, String filePathName) {
        save(psam, format, filePathName, false);
    }


    public void save(WeightMatrix psam, String format, String filePathName, boolean rnaBinding) {
        if (format.equalsIgnoreCase("xml")) {
            WeightMatrixTools.writeToXML(psam, filePathName, rnaBinding);
        }
        else if (format.equalsIgnoreCase("table")) {
            WeightMatrixTools.writeToTable(psam, filePathName, rnaBinding);
        }
        else if (format.equalsIgnoreCase("mxr")) {
            //WeightMatrixTools.writeToMatrixReduce(psam, filePathName);
        }
        else if (format.equalsIgnoreCase("uniprobe")) {
            //WeightMatrixTools.writeToUniprobe(psam, filePathName);
        }
        else {
            FileTools.writeSerializedFile(psam, filePathName);
        }
    }

    public void save(double[][] aPositionalWeights) {
    }

    private String format(double num, int minDecimalPlaces, int maxDecimalPlaces) {
        NumberFormat format = NumberFormat.getInstance();
        format.setMinimumFractionDigits(minDecimalPlaces);
        format.setMaximumFractionDigits(maxDecimalPlaces);
        return(format.format(num));
    }

    public String score(FeaturedWeightMatrix fsam,
        String scoreType,
        String[] aSeqsArray,
        SymbolList[] refSymListArray,
        double[] referenceKaArray,
        String aCompType,
        String[] aCompsArray)
    {
        double refKaReal = Double.NaN;
        //System.out.println("\naSeqsArray = "+StringTools.toString(aSeqsArray, "\t"));

        fsam.setNonSpecKa(0);

        StringBuffer stringBuffer = new StringBuffer();
        try {
            if (parser == null) {
                parser = alphabet.getTokenization("token");
            }

//             stringBuffer.append("label\tKa\trelKa\tddG");

            Table scoresTable = new Table();

            ArrayList<Double> psamScores = new ArrayList<Double>();
            ArrayList<Double> additiveScores = new ArrayList<Double>();
            ArrayList<Double> multScores = new ArrayList<Double>();
            ArrayList<Double> compScores = new ArrayList<Double>();

            //for (String aSeqString : aSeqsArray) {
            for (int i=0; i < aSeqsArray.length; i++) {
                String aSeqString = aSeqsArray[i];

                SymbolList refSymList = null;
                if (refSymListArray != null) {
                    refSymList = refSymListArray[i];
                }

                double referenceKa = Double.NaN;
                if (referenceKaArray != null) {
                    referenceKa = referenceKaArray[i];
                }

                double refKaPsamFitted = Double.NaN;
                double refKaMultFitted = Double.NaN;
                double refKaAddFitted = Double.NaN;
                if (refSymList != null) {
                    //double refKaFitted = getKaRef(fsam, refSymList, null, null);
                    refKaPsamFitted = fsam.psamScore(refSymList);
                    refKaMultFitted = fsam.score(refSymList, (FeatureKey)null);
                    refKaAddFitted = fsam.addModelScore(refSymList, (FeatureKey)null);
                }

                // get the sequences to score
                String[] seqStrings = null;
                String[] labelStrings = null;

                if (scoreType.equalsIgnoreCase("fasta")) {
                    System.out.println("\nLoading sequences to score from the fasta file " + aSeqString + ".");
                    BufferedReader bufferedReader = new BufferedReader(new FileReader(aSeqString));
                    SequenceIterator seqIter = SeqIOTools.readFasta(bufferedReader, parser);
                    ArrayList seqArrayList = new ArrayList();
                    ArrayList labelArrayList = new ArrayList();
                    while (seqIter.hasNext()) {
                        Sequence aSeq = seqIter.nextSequence();
                        seqArrayList.add(aSeq.seqString());
                        labelArrayList.add(aSeq.getName());
                    }
                    seqStrings = (String[])seqArrayList.toArray(new String[0]);
                    labelStrings = (String[])labelArrayList.toArray(new String[0]);
                }
                else if (scoreType.equalsIgnoreCase("table")) {
                    System.out.println("\nLoading sequences to score from the single-column table file " + aSeqString + ".");
                    // skip the first line
                    seqStrings = FileTools.readStrings(aSeqString, 1, "#");
                }
                else { // aSeqsArray contains actual sequences
                    seqStrings = aSeqsArray;
                }

                // get the comparison Ka's
                double[] compKas = null;
                if (aCompType != null) {
                    String aCompString = aCompsArray[i];
                    if (aCompType.equalsIgnoreCase("table")) {
                        System.out.println("\nLoading comparison Kas from the two-column table file " + aCompString + ".");
                        // skip the first line
                        String[] lineStrings = FileTools.readStrings(aCompString, 1, "#");
                        compKas = StringTools.toDoubles(StringTools.getColumn(lineStrings, this.whiteSpace, "#", 1).toArray(new String[0]));
                    }
                    else { // aSeqsArray contains actual sequences
                        compKas = StringTools.toDoubles(aCompsArray);
                    }

                    // add the comp scores
                    ArrayTools.addAll(compScores, compKas);
                }

                //for (String seqString : seqStrings) {
                for (int j=0; j < seqStrings.length;  j++) {
                    String seqString = seqStrings[j];

                    if ((seqString == null) || seqString.equals("")) {
                        continue;
                    }

                    //System.out.println(seqString.length());

                    java.util.List scoresList = new ArrayList();
                    scoresList.add(seqString);
                    SymbolList symbolList = new SimpleSymbolList (parser, seqString);
                    double kA;
                    double  relKa;

                    ///////////////////////
                    // get PSAM score
                    ///////////////////////
                    kA = fsam.psamScore(symbolList);

                    // perform necessary scaling
                    if (!Double.isNaN(refKaReal)) {
                        kA *=  (refKaReal / refKaPsamFitted);
                    }

                    if (!Double.isNaN(refKaPsamFitted)) {
                        relKa = kA / refKaPsamFitted;
                    }
                    else {
                        relKa = kA;
                    }
                    psamScores.add(relKa);
                    scoresList.add(relKa);

                    if (labelStrings == null) {
                        stringBuffer.append("\n"+relKa);
                    }
                    else {
                        stringBuffer.append("\n"+labelStrings[j]+"\t"+relKa);
                    }

                    ///////////////////////
                    // get mult score
                    ///////////////////////
                    //kA = fsam.score(symbolList, (FeatureKey)null);
                    kA = fsam.multModelScore(symbolList, (FeatureKey)null);

                    // perform necessary scaling
                    if (!Double.isNaN(refKaReal)) {
                        kA *=  (refKaReal / refKaMultFitted);
                    }

                    if (!Double.isNaN(refKaMultFitted)) {
                        relKa = kA / refKaMultFitted;
                    }
                    else {
                        relKa = kA;
                    }
                    multScores.add(relKa);
                    scoresList.add(relKa);

                    // if (labelStrings == null) {
                    //     stringBuffer.append("\n"+relKa);
                    // }
                    // else {
                    //     stringBuffer.append("\n"+labelStrings[j]+"\t"+relKa);
                    // }

                    //System.out.println(format(relKa,6,6));

                    ///////////////////////
                    // get Add score
                    ///////////////////////
                    kA = fsam.addModelScore(symbolList, (FeatureKey)null);

                    // perform necessary scaling
                    if (!Double.isNaN(refKaReal)) {
                        kA *=  (refKaReal / refKaAddFitted);
                    }

                    if (!Double.isNaN(refKaAddFitted)) {
                        relKa = kA / refKaAddFitted;
                    }
                    else {
                        relKa = kA;
                    }

                    additiveScores.add(relKa);
                    scoresList.add(relKa);

                    // if (labelStrings == null) {
                    //     stringBuffer.append("\n"+relKa);
                    // }
                    // else {
                    //     stringBuffer.append("\n"+labelStrings[j]+"\t"+relKa);
                    // }

                    if (compKas != null) {
                        // add the comp Ka
                        scoresList.add(compKas[j]);

                        //System.out.print(seqString+": compKa="+format(compKas[j],6,6)+"; ");
                    }

                    //System.out.println("relKa="+format(relKa,6,6)+"; kA(additive)="+format(kA,6,6)+"; refKaAddFitted="+format(refKaAddFitted,6,6));
                    scoresTable.add(scoresList);
                }
            }

            if (!compScores.isEmpty()) {
                // calculate R2s
                double psamR2 =     MathTools.rSquared(ArrayTools.toDoubleArray(psamScores),     ArrayTools.toDoubleArray(compScores));
                double multR2 =     MathTools.rSquared(ArrayTools.toDoubleArray(multScores),     ArrayTools.toDoubleArray(compScores));
                double additiveR2 = MathTools.rSquared(ArrayTools.toDoubleArray(additiveScores), ArrayTools.toDoubleArray(compScores));

                String fileName = fsam.getName()+".scores."+format(psamR2,3,3)+"."+format(multR2,3,3)+"."+format(additiveR2,3,3)+".table";
                FileTools.write("seq\tpsam\tmult\tadd\tcomp\n", fileName, false);
                FileTools.write(scoresTable.toString()+"\n", fileName, true);

                stringBuffer.append("\n\nPsam R2 = "+psamR2);
                stringBuffer.append("\nMult Model R2 = "+multR2);
                stringBuffer.append("\nAdditive Model R2 = "+additiveR2);

            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(stringBuffer.toString());
    }

    public void score(String[] ids,
        String dreamLabel,
        String templateFile,
        String outputFile,
        boolean kmerModelFlag,
        boolean includeFeaturesFlag,
        boolean includePosBiasFlag) {

        try {

            java.util.List<String> idsArrayList = null;
            LinkedHashSet<String> idsHashSet = null;
            if (ids != null) {
                idsArrayList = Arrays.asList(ids);
                idsHashSet = new LinkedHashSet<String>(idsArrayList);
            }

            String intensitiesDir = (String)initFileCache.get("DATA", "ProbeIntensitiesDir");
            //String probeLocsDir = (String)initFileCache.get("DATA", "ProbeLocsDir");
            String outputDir = (String)initFileCache.get("OUTPUT", "Directory");
            String resultsDir = intensitiesDir + File.separator + outputDir;

            double[] topPercentageKmersForScoring = initFileCache.getDoubles("AFFINITY MODELS", "TopPercentageKmersForScoring");
            boolean[] includeNegativeInTopPercentage = initFileCache.getBooleans("AFFINITY MODELS", "IncludeNegativeInTopPercentage");

            String aLine;
            String[] columnStrings;
            BufferedReader inBuffer = new BufferedReader(new InputStreamReader(new FileInputStream(templateFile)));
            BufferedWriter outBuffer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFile, false)));

            String oldID = null;
            FeaturedWeightMatrix fsam = null;

            // String modelString = null;
            // if (kmerModelFlag) {
            //     modelString = "\"All-Kmer Model\"";
            // }
            // else {
            //     modelString = "\"Single Motif Model\"";
            // }

            //System.out.print("\nScoring the sequences in the file " + templateFile + " using the "+ modelString +" and writing the output to "+ outputFile +"...");
            System.out.print("\nScoring the sequences in the file " + templateFile + " and writing the output to "+ outputFile +"...");
            System.out.print("\n\tScoring Options: IncludeKmers="+kmerModelFlag+", IncludeFeatures="+includeFeaturesFlag+", IncludePosBiasProfile="+includePosBiasFlag+"...");

            // skip the first row which is just column labels
            boolean firstRow = true;

            while((aLine = inBuffer.readLine()) != null) {

                if ((aLine != null) && !aLine.equals("")) {

                    if (firstRow) {
                        // write out the column header (first row) as is
                        outBuffer.write(aLine+"\n");
                        firstRow = false;
                    }
                    else {
                        columnStrings = tabPattern.split(aLine);

                        String newID = columnStrings[0];
                        if ( (idsHashSet != null) && (!idsHashSet.contains(newID)) ) {
                            continue;
                        }

                        if ( (oldID == null) || (!newID.equals(oldID)) ) {
                            System.out.print("\n\tfor protein " + newID + "...");

                            String fsamPathName = resultsDir+File.separator+dreamLabel+"."+newID+".fsam.6nt.ser";
                            if (!FileTools.exists(fsamPathName)) {
                                fsamPathName = resultsDir+File.separator+dreamLabel+"."+newID+".fsam.7nt.ser";
                            }
                            if (!FileTools.exists(fsamPathName)) {
                                fsamPathName = resultsDir+File.separator+dreamLabel+"."+newID+".fsam.8nt.ser";
                            }
                            if (!FileTools.exists(fsamPathName)) {
                                fsamPathName = resultsDir+File.separator+dreamLabel+"."+newID+".fsam.10nt.ser";
                            }

                            // String fsamPathName = resultsDir+File.separator+dreamLabel+"."+newID+".fsam.12nt.ser";
                            // if (!FileTools.exists(fsamPathName)) {
                            //     fsamPathName = resultsDir+File.separator+dreamLabel+"."+newID+".fsam.10nt.ser";
                            // }
                            // if (!FileTools.exists(fsamPathName)) {
                            //     fsamPathName = resultsDir+File.separator+dreamLabel+"."+newID+".fsam.8nt.ser";
                            // }

                            fsam = (FeaturedWeightMatrix)FileTools.readSerializedFile(fsamPathName);
                            fsam.reHashFeatures();
                            fsam.cleanKmerAffinities();
                            fsam.setScoringKmerAffinities(topPercentageKmersForScoring, includeNegativeInTopPercentage);

                            char[] kmerLengths = fsam.getKmerLengths();
                            //if (kmerLengths != null) {
                            if (fsam.getKmerToAffinityMatrix() != null) {
                                System.out.print("\tKmerLengths: "+StringTools.toString(ArrayTools.toIntArray(kmerLengths), ", "));

                                System.out.print("\tRev-comps Equivalent: "+StringTools.toString(fsam.includeRevComps, ", "));

                                System.out.print("    Non-negative: ");
                                for (int j=0; j < kmerLengths.length; j++) {
                                    double[] kmerAffinities = fsam.getKmerToAffinityMatrix(j);
                                    boolean hasNegatives = MathTools.hasLessThan(kmerAffinities, 0);
                                    if (hasNegatives) {
                                        System.out.print("F, ");
                                    }
                                    else {
                                        System.out.print("T, ");
                                    }
                                }

                                for (int j=0; j < kmerLengths.length; j++) {
                                    double[][] kmerPositionalWeights = fsam.getKmerPositionalWeights(j);
                                    if (kmerPositionalWeights != null) {
                                        out.println("\nIndex "+j+" positional weights are:\n\t"+StringTools.toString(kmerPositionalWeights, "\n\t", ", "));
                                    }
                                }


                            }


                            if (initPositionalBias != null) {
                                if (kmerModelFlag) {
                                    //fsam.kmerPositionalWeights = initPositionalBias;
                                }
                                else {
                                    fsam.setPositionalWeights(initPositionalBias);
                                }
                            }

                            oldID = newID;
                        }

                        //double score = score(fsam, columnStrings[2], 1.0, Double.NaN);


                        // numKmerWindows = probeSeqLengths + motifLength - 1;
                        int trainingProbeSeqLengths = fsam.positionalWeights[0].length - fsam.columns() + 1;
                        String seqToScore = null;

                        if (includePosBiasFlag) {
                            // if the scoringSeq is longer than the trainingSeq then truncate the scoringSeq
                            if (trainingProbeSeqLengths < columnStrings[2].length()) {
                                //seqToScore = columnStrings[2].substring(0, 35);
                                seqToScore = columnStrings[2].substring(0, trainingProbeSeqLengths);
                            }
                            // else use the scoring seq as is
                            else {
                                seqToScore = columnStrings[2];
                            }
                        }
                        // else use the scoring seq as is
                        else {
                            seqToScore = columnStrings[2];
                        }

                        SymbolList aSymList = DNATools.createDNA(seqToScore);

                        /////////////////////////////////////////////////
                        //fsam.setNonSpecKa(.00001);
                        fsam.setNonSpecKa(0);
                        /////////////////////////////////////////////////

                        double score = -1;

                        //score = fsam.kmerScore(aSymList);

                        //score = fsam.oldPsamScore(aSymList);
                        //score = fsam.psamScore(aSymList);
                        //score = fsam.psamScore(aSymList, fsam.getPositionalWeights());
                        //score = fsam.score(aSymList, fsam.getPositionalWeights(), (FeatureKey)null);
                        //score = fsam.addModelScore(aSymList, fsam.getPositionalWeights(), (FeatureKey)null);

                        score = fsam.modelScore(aSymList, true, kmerModelFlag, includeFeaturesFlag, includePosBiasFlag);

                        columnStrings[3] = Double.toString(score);
                        outBuffer.write(StringTools.toString(columnStrings, "\t")+"\n");
                    }
                }
            }
            inBuffer.close();
            outBuffer.close();
            System.out.println("Done.");

        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

// 	System.out.print("\nWriting the scores to the output file " + outputFile + "...");
// 	FileTools.write(StringTools.toString(templateFileMatrix, "\n", "\t"), outputFile, false);
// 	System.out.println("Done.");
    }

    public double[] score(FeaturedWeightMatrix fsam, String[] seqStrings, double refKaFitted, double refKaReal) {
        ArrayList arrayList = new ArrayList();

        for (int i=0; i < seqStrings.length; i++) {
            if ((seqStrings[i] != null) && !seqStrings[i].equals("")) {
                //stringBuffer.append("\n"+score(fsam, seqStrings[i], refKaFitted, refKaReal));
                arrayList.add(score(fsam, seqStrings[i], refKaFitted, refKaReal));
            }
        }
        return(ArrayTools.toDoubleArray(arrayList));
    }

    public double score(FeaturedWeightMatrix fsam, String seqString, double refKaFitted, double refKaReal) {
        try {
            if (parser == null) {
                parser = alphabet.getTokenization("token");
            }
            SymbolList symbolList = new SimpleSymbolList (parser, seqString);
            return(score(fsam, symbolList, seqString, refKaFitted, refKaReal));
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(Double.NaN);
    }

    public double score(FeaturedWeightMatrix fsam, SymbolList aSymList, String label, double refKaFitted, double refKaReal) {
        try {
            double kA = fsam.score(aSymList, (FeatureKey)null);

            // perform necessary scaling
            if (!Double.isNaN(refKaReal)) {
                kA *=  (refKaReal / refKaFitted);
            }

            double  relKa = kA / refKaFitted;
            //double ddG = WeightMatrixTools.getGibbsFreeEnergy(kA, refKaFitted);
            //return(label+"\t"+kA+"\t"+relKa+"\t"+ddG);
            return(relKa);
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        //return(null);
        return(Double.NaN);
    }

    public double[] score(FeaturedWeightMatrix fsam, BufferedReader fastaBufferedReader, double refKaFitted, double refKaReal) {
        try {
            if (parser == null) {
                parser = alphabet.getTokenization("token");
            }
            SequenceIterator seqIter = SeqIOTools.readFasta(fastaBufferedReader, parser);
            return(score(fsam ,seqIter, refKaFitted, refKaReal));
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(null);
    }

    public double[] score(FeaturedWeightMatrix fsam, SequenceIterator seqIter, double refKaFitted, double refKaReal) {
        ArrayList arrayList = new ArrayList();
        try {
            while (seqIter.hasNext()) {
                Sequence aSeq = seqIter.nextSequence();
                arrayList.add(score(fsam ,aSeq, aSeq.getName(), refKaFitted, refKaReal));
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(ArrayTools.toDoubleArray(arrayList));
    }

    // Used for both HammingREDUCE and columnREDUCE
    // pos = 0-based indexing
    // Now returns a Relative affinity PSAM (not a probability PWM)
    public Distribution getPsamColumn(
        RegressionFitData aRegressionFitData,
        int pos,
        double minAffinityWeight,
        Distribution background,
        double eToMu) {

        Distribution psamColumn = null;
        try {
            //boolean useRegressionIntercept = (boolean)initFileCache.getBoolean("AFFINITY MODELS", "UseRegressionIntercept");
            //double intercept = aRegressionFitData.getIntercept();

            psamColumn = DistributionFactory.DEFAULT.createDistribution(this.alphabet);

            double highestMagnitude = 0;
            double properSign = 0;
            // double highestMagnitude = getHighestMagnitude(aRegressionFitData);
            // double properSign = Math.signum(highestMagnitude);
            //int numWithHighestMagnitude = 0;

            //BigDecimal spKaSum = new BigDecimal(0.0);
            //BigDecimal spKaSlack = new BigDecimal(0.0);
            ArrayList<Symbol> symbols = new ArrayList<Symbol>(4);
            ArrayList<Double> weights = new ArrayList<Double>(4);

            // The seed should always be the first one with a coefficient
            Set<Object> hammingDistHashSet = aRegressionFitData.motifsWithCoeffsMap.keySet();
            Object seed = hammingDistHashSet.iterator().next();

            // Get the total sum of the coefficients in this column
            for (Iterator dnaIter = ((FiniteAlphabet)alphabet).iterator(); dnaIter.hasNext(); ) {
                Symbol symbol = (Symbol)dnaIter.next();
                //Symbol[] symbolArray = {symbol};

                //make a copy of the seed
                Object seedCopy;
                Object revComp;

                if (seed instanceof SymbolList) {
                    seedCopy = new SimpleSymbolList((SymbolList)seed);
                    //Edit edit = new Edit(pos+1, 1, new SimpleSymbolList(symbolArray, 1, alphabet));
                    Edit edit = new Edit(pos+1, alphabet, symbol);
                    ((SymbolList)seedCopy).edit(edit);
                    revComp = DNATools.reverseComplement((SymbolList)seedCopy);
                }
                else { // it's a WeightMatrix
                    seedCopy = WeightMatrixTools.mutate((WeightMatrix) seed, pos, symbol, true);
                    revComp = WeightMatrixTools.reverseComplement((WeightMatrix) seedCopy, complementTable);
                }

                Double coefObject = aRegressionFitData.motifsWithCoeffsMap.get(seedCopy);

                // if the coefficient is Null then check reverse complement
                // if it does have a rev comp in the set then grab those coefficients
                if (coefObject == null) {
                    coefObject = aRegressionFitData.motifsWithCoeffsMap.get(revComp);
                    if (coefObject == null) {
                        coefObject = new Double(0);
                    }
                }

                double coef = coefObject.doubleValue();

                // HACK!!!
                // Poisson Regression sometimes returns negative numbers of high magnitude
                // Forcing return of highest positive number!
                if (coef > highestMagnitude) {
                    highestMagnitude = coef;
                    properSign = Math.signum(highestMagnitude);
                }

                double spKa;
                if (eToMu == 0) {
                    spKa= coef;
                }
                else {
                    spKa = coef / (eToMu * (1 - coef));
                }

                symbols.add(symbol);
                weights.add(spKa);

            }

//                 if (useRegressionIntercept) {
//                     spKa += intercept;
//                 }

            // take the nth root of the normalized coefs for Poisson regression
            String method = (String)initFileCache.get("REGRESSION", "LinearAlgorithm");
            double numSelexRounds = 1;

            if (method.equalsIgnoreCase("glm.poisson")) {
                numSelexRounds = initFileCache.getInt("INTENSITIES", "NumSelexRounds");
            }

            // Calculate the weight = spKa / highestMagnitude
            for (int j=0; j < weights.size(); j++) {

                Symbol symbol = symbols.get(j);
                double spKa = weights.get(j);

                // if all the coefs are zero, then there's no specificity and set spKA=1
                if (highestMagnitude == 0) {
                    spKa = 1.0;
                }
                else if (Math.signum(spKa) != properSign) {
                    spKa = properSign * minAffinityWeight * background.getWeight(symbol);
                }
                else {
                    if (numSelexRounds == 1) {
                        spKa = spKa / highestMagnitude;
                    }
                    else { // numSelexRounds > 1
                        spKa = Math.pow(spKa / highestMagnitude, 1.0/numSelexRounds);
                    }

                    if (spKa < minAffinityWeight * background.getWeight(symbol)) {
                        spKa = properSign * minAffinityWeight * background.getWeight(symbol);
                    }
                }

                psamColumn.setWeight(symbol, spKa);

//                 if (spKa == highestMagnitude) {
//                     numWithHighestMagnitude++;
//                 }

//                 // If the spKa is NOT the proper sign or < minAffinityWeight, then it is noise so set it to minAffinityWeight * symbolWeight
//                 // also subract the spKa from the spKaSum
//                 if (Math.signum(spKa) != properSign) {
//                     spKa = properSign * minAffinityWeight * background.getWeight(symbol);
//                     psamColumn.setWeight(symbol, spKa);
//                     spKaSlack = spKaSlack.add(new BigDecimal(spKa));
//                 }
//                 else if ((spKa/highestMagnitude) < (minAffinityWeight * background.getWeight(symbol))) {
//                     spKa = properSign * minAffinityWeight * background.getWeight(symbol);
//                     psamColumn.setWeight(symbol, spKa);
//                     spKaSlack = spKaSlack.add(new BigDecimal(spKa));
//                 }
//                 else {
//                     // add the spKa to the spKa sum
//                     symbols.add(symbol);
//                     weights.add(spKa);
//                     spKaSum = spKaSum.add(new BigDecimal(spKa));
//                 }

            }

//             //double spKaSlackPerWeight = spKaSlack.divide(new BigDecimal(weights.size())).doubleValue();
//             //double spKaSlackPerWeight = spKaSlack.doubleValue() / weights.size();

//             double spKaSlackPerHighestMag = spKaSlack.doubleValue() / numWithHighestMagnitude;

//             // Calculate the weight = spKa / Sum(spKas)
//             for (int j=0; j < weights.size(); j++) {

//                 if (weights.get(j) < highestMagnitude) {
//                     psamColumn.setWeight(symbols.get(j), (weights.get(j) / spKaSum.doubleValue()));
//                 }
//                 else {
//                     psamColumn.setWeight(symbols.get(j), (weights.get(j) / spKaSum.doubleValue()) - spKaSlackPerHighestMag);
//                 }


//                 //psamColumn.setWeight(symbols.get(j), (weights.get(j) / spKaSum.doubleValue()) - spKaSlackPerWeight);

//                 //psamColumn.setWeight(symbols.get(j), (weights.get(j) + intercept) / (spKaSum.doubleValue() + intercept));
//                 //psamColun.setWeight(symbols.get(j), weights.get(j) / spKaSum);
//             }

        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(psamColumn);
    }

    public double getHighestMagnitude(RegressionFitData aRegressionFitData) {
            // Remove intercept from the hammingCoefs
            ArrayList hammingCoefs = new ArrayList(aRegressionFitData.motifsWithCoeffsMap.values());

            // get proper sign
            int properSign = 1;
            double min = MathTools.min(hammingCoefs);
            double max = MathTools.max(hammingCoefs);


            // HACK!!!
            // Poisson Regression somteimes returns negative numbers of high magnitude
            // Forcing return of highest positive number!
            return(max);


//             if (Math.abs(min) > max) {
//                 //properSign = -1;
//                 return(min);
//             }
//             //return(properSign);
//             return(max);
    }

    // perform columnREDUCE
    // If the PSAM is not self-reverse-complement then perform column regression on every column
    // else perform column regression on just the 1st half, and copy the complements for the 2nd half
    public WeightMatrix getPsam(
        String method,
        WeightMatrix seedPSAM,
        double[][] positionalWeights,
        double minAffinityWeight,
        Distribution background,
        WeightMatrixTools.BindingStrand effectiveStrand,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double eToMu,
        SymbolList refSymList,
        double referenceKa,
        int unifiedColumn,
        MutableDouble alpha,
        MutableDouble beta,
        MutableDouble nonSpecKa) {
        return(getPsam(
                method,
                seedPSAM,
                positionalWeights,
                minAffinityWeight,
                background,
                effectiveStrand,
                strand,
                calc,
                eToMu,
                refSymList,
                referenceKa,
                unifiedColumn,
                alpha,
                beta,
                nonSpecKa,
                null));

    }

    // called to re-fit a PSAM and step in that direction
    public WeightMatrix getPsam(
        String method,
        WeightMatrix seedPSAM,
        double[][] positionalWeights,
        double minAffinityWeight,
        Distribution background,
        WeightMatrixTools.BindingStrand effectiveStrand,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double eToMu,
        SymbolList refSymList,
        double referenceKa,
        int unifiedColumn,
        MutableDouble alpha,
        MutableDouble beta,
        MutableDouble nonSpecKa,
        boolean[] positions) {

        WeightMatrix newPsam = seedPSAM;

        int numIters = initFileCache.getInt("REGRESSION", "NumNewPsamIterations");
        if (numIters == 0) {
            numIters = 5;
        }

        for (int i=0; i<numIters; i++) {
        // for (int i=0; i<8; i++) {
            //for (int i=0; i<1; i++) {

            WeightMatrix iterPsam = getPsamIter(
                method,
                //seedPSAM,
                newPsam,
                positionalWeights,
                minAffinityWeight,
                background,
                effectiveStrand,
                strand,
                calc,
                eToMu,
                refSymList,
                referenceKa,
                unifiedColumn,
                alpha,
                beta,
                nonSpecKa,
                positions);

            //newPsam = iterPsam;

            // if positions!=null then set the starting PSAM to iterPSAM
            if ((i == 0) && (positions != null)) {
                //if (i == 0) {
                newPsam = iterPsam;
            }
            // else perform a step in the iterPSAM direction
            else {
                // weightedAverage(double weight1, double weight2, WeightMatrix pwm1, WeightMatrix pwm2)
                newPsam = WeightMatrixTools.weightedAverage(3, 1, newPsam, iterPsam);
                //newPsam = WeightMatrixTools.weightedAverage(2, 1, newPsam, iterPsam);
            }

            // Correspond Positions if they exist
            if (this.correspondenceStrings != null) {
                // out.println("\nBefore correspondMotif():");
                // out.println(WeightMatrixTools.toString(newPsam));

                WeightMatrixTools.correspondMotif(newPsam, this.correspondenceStrings, this.complementTable);
                newPsam = WeightMatrixTools.getRelativeAffinities(newPsam);

                // out.println("\nAfter correspondMotif():");
                // out.println(WeightMatrixTools.toString(newPsam));
            }

        }

        newPsam.setName("PSAM");
        return(newPsam);
    }


    // public WeightMatrix getPsam(
    //     String method,
    //     WeightMatrix seedPSAM,
    //     double[][] positionalWeights,
    //     double minAffinityWeight,
    //     Distribution background,
    //     WeightMatrixTools.BindingStrand effectiveStrand,
    //     WeightMatrixTools.BindingStrand strand,
    //     WeightMatrixTools.BothStrandsCalc calc,
    //     double eToMu,
    //     SymbolList refSymList,
    //     double referenceKa,
    //     int unifiedColumn,
    //     MutableDouble alpha,
    //     MutableDouble beta,
    //     MutableDouble nonSpecKa,
    //     boolean[] positions) {

    //     WeightMatrix newPsam = seedPSAM;

    //     //for (int i=0; i<8; i++) {
    //     for (int i=0; i<1; i++) {

    //         WeightMatrix iterPsam = getPsamIter(
    //             method,
    //             //seedPSAM,
    //             newPsam,
    //             positionalWeights,
    //             minAffinityWeight,
    //             background,
    //             effectiveStrand,
    //             strand,
    //             calc,
    //             eToMu,
    //             refSymList,
    //             referenceKa,
    //             unifiedColumn,
    //             alpha,
    //             beta,
    //             nonSpecKa,
    //             positions);

    //         newPsam = iterPsam;

    //         // // if positions!=null then set the starting PSAM to iterPSAM
    //         // if ((i == 0) && (positions != null)) {
    //         //     //if (i == 0) {
    //         //     newPsam = iterPsam;
    //         // }
    //         // // else perform a step in the iterPSAM direction
    //         // else {
    //         //     // weightedAverage(double weight1, double weight2, WeightMatrix pwm1, WeightMatrix pwm2)
    //         //     newPsam = WeightMatrixTools.weightedAverage(3, 1, newPsam, iterPsam);
    //         //     //newPsam = WeightMatrixTools.weightedAverage(2, 1, newPsam, iterPsam);
    //         // }
    //     }

    //     newPsam.setName("PSAM");
    //     return(newPsam);
    // }


    public WeightMatrix getPsamIter(
        String method,
        WeightMatrix seedPSAM,
        double[][] positionalWeights,
        double minAffinityWeight,
        Distribution background,
        WeightMatrixTools.BindingStrand effectiveStrand,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double eToMu,
        SymbolList refSymList,
        double referenceKa,
        int unifiedColumn,
        MutableDouble alpha,
        MutableDouble beta,
        MutableDouble nonSpecKa,
        boolean[] positions) {

        try {
            //String linearAlgorithm = (String)initFileCache.get("REGRESSION", "LinearAlgorithm");
            RegressionFitData regressionFitData = null;
            LinkedHashMap<Object, LinkedHashSet<SymbolList>> motifToHammingSphereMap = null;
            int seedPSAMlength = seedPSAM.columns();
            WeightMatrix newPSAM = WeightMatrixTools.clone(seedPSAM);

            boolean isRevCompPalindrome = WeightMatrixTools.isSelfReverseComplement(seedPSAM, complementTable);

            int endColumn;
            if (isRevCompPalindrome) {
                endColumn = seedPSAMlength / 2;
                /////////////////////////////////////////////////////////////
                // Find out why this is needed!!!
                //   if (effectiveStrand == WeightMatrixTools.BindingStrand.BOTH) {
                //       //strand = WeightMatrixTools.BindingStrand.NEG;
                //       strand = WeightMatrixTools.BindingStrand.POS;
                //   }
                /////////////////////////////////////////////////////////////
            }
            else {
                endColumn = seedPSAMlength;
            }

            for (int pos = 0; pos < endColumn; pos++) {
                //for (int pos = 0; pos < 1; pos++) {

                // look to see if we should perform columnREDUCE on this column
                if ((positions != null) && (positions[pos] == false)) {
                    continue;
                }

                motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();

                for (Iterator dnaIter = ((FiniteAlphabet)alphabet).iterator(); dnaIter.hasNext(); ) {
                    Symbol mutuateSymbol = (Symbol)dnaIter.next();

                    // create mutation in the PSAM
                    WeightMatrix singleMutationPSAM = WeightMatrixTools.mutate(newPSAM, pos, mutuateSymbol);
                    //WeightMatrix singleMutationPSAM = WeightMatrixTools.mutate(seedPSAM, pos, mutuateSymbol);

                    // name has 1-based indexing
                    //singleMutationPSAM.setName("[I("+ (pos+1) +"-"+(pos+2)+","+1+","+mutuateSymbol.getName().substring(0,1).toUpperCase()+")]");
                    singleMutationPSAM.setName("mutate_"+ (pos+1) +"_"+mutuateSymbol.getName().substring(0,1).toUpperCase()+"");

                    LinkedHashSet<SymbolList> hammingSphere = null;
                    motifToHammingSphereMap.put(singleMutationPSAM, hammingSphere);

                    //out.println((byte)0, "WeightMatrix: "+StringTools.toString(singleMutationPSAM));
                    out.println((byte)0, "WeightMatrix: "+WeightMatrixTools.toString(singleMutationPSAM));
                    out.println((byte)0, "**********************************************************");
                }

                // if method = nls.levMar then must populate startParams
                double[] startParams = null;
                if (method.equalsIgnoreCase("nls.levMar")) {
                    // startParams is in the format (Ka(A), Ka(C), Ka(G), Ka(T))
                    startParams = new double[4];
                    //startParams[0] = alpha.doubleValue();
                    //startParams[1] = beta.doubleValue();
                    //startParams[2] = nonSpecKa.doubleValue();
                    // initialize the Ka's all to 1.0
                    Arrays.fill(startParams, 0, startParams.length, 1.0);

                    this.globalAlpha = alpha.doubleValue();
                    this.globalBeta = beta.doubleValue();
                    this.globalNonSpecKa = nonSpecKa.doubleValue();;
                }

                // perform lm regression to find the coefficients for the column
                regressionFitData = fitModel(
                    motifToHammingSphereMap,
                    null, //startPositions
                    positionalWeights,
                    null, //motifsToMods
                    null, //motifsToMandatories
                    "PSAM column regression at column "+(pos+1),
                    method,
                    effectiveStrand,
                    strand,
                    calc,
                    startParams, //startParams
                    refSymList, //refSymList
                    referenceKa,
                    unifiedColumn); //refKaReal
//                     null, // alpha
//                     null); // nonSpecKa

                if (regressionFitData == null) {
                    out.println("\nRegressionFitData is null!! Skipping this column!");
                    continue;
                }

                // coefficients are as returned in R
                // First coefficient is the intercept, then the rest are in independent-var order
                out.print(regressionFitData.toString());

                if (MathTools.hasInfinite(regressionFitData.getCoefficients())) {
                    // regression failed!
                    out.println("\nNew Coefficients contain an infinite value, not changing the old ones!");
                }
                else {
                    //dists[pos] = getPsamColumn(regressionFitData, pos, minAffinityWeight, background, eToMu);
                    newPSAM.setColumn(pos, getPsamColumn(regressionFitData, pos, minAffinityWeight, background, eToMu));
                }
            }

            // if revCompPalidrome then need to copy the reverseComplement of the
            // second half of the newPSAM
            if (isRevCompPalindrome) {
                out.println("");
                for (int destPos = endColumn; destPos < seedPSAMlength; destPos++) {
                    int sourcePos = (seedPSAMlength - 1) - destPos;
                    out.println("Copying complement of column "+(sourcePos+1)+" distribution to column "+(destPos+1)+".");
                    //dists[destPos] = DistributionTools.complement(dists[sourcePos], complementTable);
                    newPSAM.setColumn(destPos, DistributionTools.complement(newPSAM.getColumn(sourcePos), complementTable));
                }
            }
            newPSAM.setName("PSAM");

            // Correspond Positions if they exist
            if (this.correspondenceStrings != null) {
                // out.println("\nBefore correspondMotif():");
                // out.println(WeightMatrixTools.toString(newPSAM));

                WeightMatrixTools.correspondMotif(newPSAM, this.correspondenceStrings, this.complementTable);
                newPSAM = WeightMatrixTools.getRelativeAffinities(newPSAM);

                // out.println("\nAfter correspondMotif():");
                // out.println(WeightMatrixTools.toString(newPSAM));
            }

            return(newPSAM);
            //return(new SimpleWeightMatrix("PSAM", dists));
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(null);

    }

//     public WeightMatrix getPsam(
//         WeightMatrix seedPSAM,
//         double[][] positionalWeights,
//         double minAffinityWeight,
//         Distribution background,
//         WeightMatrixTools.BindingStrand effectiveStrand,
//         WeightMatrixTools.BindingStrand strand,
//         WeightMatrixTools.BothStrandsCalc calc,
//         double eToMu) {

//         try {
//             String linearAlgorithm = (String)initFileCache.get("REGRESSION", "LinearAlgorithm");
//             RegressionFitData regressionFitData = null;
//             LinkedHashMap<Object, LinkedHashSet<SymbolList>> motifToHammingSphereMap = null;
//             int seedPSAMlength = seedPSAM.columns();
//             Distribution[] dists = new Distribution[seedPSAMlength];

//             boolean isRevCompPalindrome = WeightMatrixTools.isSelfReverseComplement(seedPSAM, complementTable);

//             int endColumn;
//             if (isRevCompPalindrome) {
//                 endColumn = seedPSAMlength / 2;
//                 /////////////////////////////////////////////////////////////
//                 // Find out why this is needed!!!
// //                 if (effectiveStrand == WeightMatrixTools.BindingStrand.BOTH) {
// //                     strand = WeightMatrixTools.BindingStrand.POS;
// //                 }
//                 /////////////////////////////////////////////////////////////
//             }
//             else {
//                 endColumn = seedPSAMlength;
//             }

//             for (int pos = 0; pos < endColumn; pos++) {

//                 motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();

//                 for (Iterator dnaIter = ((FiniteAlphabet)alphabet).iterator(); dnaIter.hasNext(); ) {
//                     Symbol mutuateSymbol = (Symbol)dnaIter.next();

//                     // create mutation in the PSAM
//                     WeightMatrix singleMutationPSAM = WeightMatrixTools.mutate(seedPSAM, pos, mutuateSymbol);

//                     // name has 1-based indexing
//                     //singleMutationPSAM.setName("[I("+ (pos+1) +"-"+(pos+2)+","+1+","+mutuateSymbol.getName().substring(0,1).toUpperCase()+")]");
//                     singleMutationPSAM.setName("mutate_"+ (pos+1) +"_"+mutuateSymbol.getName().substring(0,1).toUpperCase()+"");

//                     LinkedHashSet<SymbolList> hammingSphere = null;
//                     motifToHammingSphereMap.put(singleMutationPSAM, hammingSphere);

//                     out.println((byte)0, "WeightMatrix: "+StringTools.toString(singleMutationPSAM));
//                     out.println((byte)0, "**********************************************************");
//                 }

//                 // perform lm regression to find the coefficients for the column
//                 regressionFitData = fitModel(
//                     motifToHammingSphereMap,
//                     null,
//                     positionalWeights,
//                     null,
//                     null,
//                     "PSAM column regression at column "+(pos+1),
//                     linearAlgorithm,
//                     effectiveStrand,
//                     strand,
//                     calc,
//                     null);

//                 // coefficients are as returned in R
//                 // First coefficient is the intercept, then the rest are in independent-var order
//                 out.print(regressionFitData.toString());

//                 dists[pos] = getPsamColumn(regressionFitData, pos, minAffinityWeight, background, eToMu);
//             }

//             // if revCompPalidrome then need to copy the reverseComplement of the
//             // second half of the dists
//             if (isRevCompPalindrome) {
//                 out.println("");
//                 for (int destPos = endColumn; destPos < seedPSAMlength; destPos++) {
//                     int sourcePos = (seedPSAMlength - 1) - destPos;
//                     out.println("Copying complement of column "+(sourcePos+1)+" distribution to column "+(destPos+1)+".");
//                     dists[destPos] = DistributionTools.complement(dists[sourcePos], complementTable);
//                 }
//             }

//             return(new SimpleWeightMatrix("PSAM", dists));
//         }
//         catch (Exception ex) {
//             ex.printStackTrace();
//         }
//         return(null);

//     }

    // perform HammingREDUCE (and then columnREDUCE)
    public WeightMatrix getPsam(
        String method,
        RegressionFitData aRegressionFitData,
        double minAffinityWeight,
        Distribution background,
        WeightMatrixTools.BindingStrand effectiveStrand,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double eToMu,
        SymbolList refSymList,
        double referenceKa,
        int unifiedColumn,
        MutableDouble alpha,
        MutableDouble beta,
        MutableDouble nonSpecKa) {

        try {

            // The seed should always be the first one with a coefficient
            Set<Object> hammingDistHashSet = aRegressionFitData.motifsWithCoeffsMap.keySet();
            SymbolList seedSymList = (SymbolList)hammingDistHashSet.iterator().next();
            Distribution[] dists = new Distribution[seedSymList.length()];

            for (int i = 0; i < seedSymList.length(); i++) {

                dists[i] = getPsamColumn(aRegressionFitData, i, minAffinityWeight, background, eToMu);

            }

            // the hammingPSAM(seedSeq) is the seed for columnPSAM(seedPSAM)
            WeightMatrix seedPSAM = new SimpleWeightMatrix("PSAM", dists);

            // Correspond Positions if they exist
            if (this.correspondenceStrings != null) {
                // out.println("\nBefore correspondMotif():");
                // out.println(WeightMatrixTools.toString(seedPSAM));

                WeightMatrixTools.correspondMotif(seedPSAM, this.correspondenceStrings, this.complementTable);
                seedPSAM = WeightMatrixTools.getRelativeAffinities(seedPSAM);

                // out.println("\nAfter correspondMotif():");
                // out.println(WeightMatrixTools.toString(seedPSAM));
            }

            boolean iterateColumnwisePsamRegression = (boolean)initFileCache.getBoolean("REGRESSION", "IterateColumnwisePsamRegression");
            if (!iterateColumnwisePsamRegression) {
                return(seedPSAM);
            }
            else {
                // Now perform column-wise regression again with derived PSAM
                return(getPsam(
                        method,
                        seedPSAM,
                        null,
                        minAffinityWeight,
                        background,
                        effectiveStrand,
                        strand,
                        calc,
                        eToMu,
                        refSymList,
                        referenceKa,
                        unifiedColumn,
                        alpha,
                        beta,
                        nonSpecKa));
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(null);
    }

// hammingCoefs SHOULD NOT contain the intercept!!
//     public WeightMatrix getPsam(LinkedHashSet<SymbolList> hammingDistHashSet, double[] hammingCoefs, Distribution background) {

//         try {
//             // get proper sign
//             double noise = 0;
//             double min = MathTools.min(hammingCoefs);
//             double max = MathTools.max(hammingCoefs);
//             if (max>Math.abs(min) && min<0) {
//                 noise = Math.abs(min);
//             }
//             else if (Math.abs(min)>max && max>0) {
//                 noise = -max;
//             }

//             // create hammingCoefsMap
//             LinkedHashMap<SymbolList, Double> hammingCoefsMap = new LinkedHashMap<SymbolList, Double>(hammingDistHashSet.size());
//             for(int j=0; j < hammingDistHashSet.size(); j++) {
//                 hammingCoefsMap.put(hammingDistHashSet.get(j), hammingCoefs[j]);
//             }

//             // The seed should always be the first one
//             SymbolList seedSymList = hammingDistHashSet.get(0);
//             Distribution[] dists = new Distribution[seedSymList.length()];

//             for (int i = 1; i <= seedSymList.length(); i++) {

//                 dists[i-1] = DistributionFactory.DEFAULT.createDistribution(this.alphabet);
//                 double coefSum = 0;
//                 ArrayList<Symbol> symbols = new ArrayList<Symbol>(4);
//                 ArrayList<Double> weights = new ArrayList<Double>(4);

//                 // Get the total sum of the coefficients in this column
//                 for (Iterator dnaIter = ((FiniteAlphabet)alphabet).iterator(); dnaIter.hasNext(); ) {
//                     Symbol symbol = (Symbol)dnaIter.next();
//                     Symbol[] symbolArray = {symbol};
//                     //make a copy of the seedSymList
//                     SimpleSymbolList seqCopy = new SimpleSymbolList(seedSymList);
//                     Edit edit = new Edit(i, 1, new SimpleSymbolList(symbolArray, 1, alphabet));
//                     seqCopy.edit(edit);
//                     double coef = hammingCoefsMap.get(seqCopy);

//                     symbols.add(symbol);
//                     coefSum += coef + noise;
//                     weights.add(coef + noise);
//                 }

//                 // Calculate the weight = coeff / Sum(coeffs)
//                 for (int j=0; j < weights.size(); j++) {
//                     dists[i-1].setWeight(symbols.get(j), weights.get(j) / coefSum);
//                 }
//             }

//             return(new SimpleWeightMatrix(dists));
//         }
//         catch (Exception ex) {
//             ex.printStackTrace();
//         }
//         return(null);
//     }

    // method = nls.levMar, lm, nnls
    // alpha = gain factor = an overall signal intensity scaler (models the light source, reflection, light meter intensity)
    // beta = constant signal intercept
    public double[][] getPositionalBias(
        String method,
        Object motif,
        WeightMatrixTools.BothStrandsCalc bothStrandsCalc,
        int[] startPositions,
        WeightMatrixTools.BindingStrand effectiveStrand,
        WeightMatrixTools.BindingStrand strand,
        SymbolList refSymList,
        double referenceKa,
        int unifiedColumn,
        MutableDouble alpha,
        MutableDouble beta,
        MutableDouble nonSpecKa) {

        ///String linearAlgorithm = (String)initFileCache.get("REGRESSION", "LinearAlgorithm");
        int probeSeqLengths = reduceData.getSeqLength();

        int motifLength = -1;
        if (motif instanceof WeightMatrix) {
            motifLength = ((WeightMatrix)motif).columns();
        }
        else if (motif instanceof FeaturedWeightMatrix) {
            motifLength = ((FeaturedWeightMatrix)motif).getPosStrandPWM().columns();
        }

        int numKmerWindows = probeSeqLengths + motifLength - 1;
        boolean isRevCompPalindrome = isSelfReverseComplement(motif, complementTable);
        RegressionFitData regressionFitData = null;
        LinkedHashMap<Object, LinkedHashSet<SymbolList>> motifToHammingSphereMap = null;

        if ((effectiveStrand == WeightMatrixTools.BindingStrand.BOTH) && (!isRevCompPalindrome)) {
            this.globalNumStrands = 2;
        }
        else {
            this.globalNumStrands = 1;
        }

        double[][] positionalCoeffs = new double[2][];

        motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();
        motifToHammingSphereMap.put(motif, null);

        // if method = nls.levMar then must populate startParams
        double[] startParams = null;
        if (method.equalsIgnoreCase("nls.levMar")) {
            // startParams is in the format (alpha, beta, nonSpecKa, concGamma1pos, concGamma2pos, ...., concGamma1neg, ....)
            startParams = new double[3 + (this.globalNumStrands*numKmerWindows)];
            startParams[0] = alpha.doubleValue();
            startParams[1] = beta.doubleValue();
            startParams[2] = nonSpecKa.doubleValue();
            // initialize the startPosition biases all to 1.0
            Arrays.fill(startParams, 3, startParams.length, 1.0);
        }

        // perform lm regression to find the probe-position bias
        regressionFitData = fitModel(
            motifToHammingSphereMap,
            startPositions,
            null, //positionalWeights
            null, //motifsToMods
            null, //motifsToMandatories
            "probe-position bias on the positive and negative strands",
            //linearAlgorithm,
            method,
            effectiveStrand,
            strand,
            bothStrandsCalc,
            startParams, //startParams
            refSymList, //refSymList
            referenceKa,
            unifiedColumn); //refKaReal
//             alpha,
//             beta,
//             nonSpecKa);

        //regressionFitData = fitModel(motifToHammingSphereMap, startPositions, null, null, "probe-position bias on the positive strand", featureDiscoveryAlgorithm, WeightMatrixTools.BindingStrand.BOTH, null);

        // coefficients are as returned in R
        // First coefficient is the intercept, then the rest are in independent-var order
        out.print(regressionFitData.toString());

        // if method is nls.levMar then populate alpha and nonspecKa
        if (method.equalsIgnoreCase("nls.levMar")) {
            alpha.setValue(regressionFitData.getAlpha());
            beta.setValue(regressionFitData.getBeta());
            nonSpecKa.setValue(regressionFitData.getNonSpecKa());
        }

        WeightMatrixTools.BindingStrand[] strandsArray = WeightMatrixTools.BindingStrand.values();

        // populate the positive strand weights

        //                     int k;
        //                     if (includeIntercept) {
        //                         k = 1;
        //                     }
        //                     else {
        //                         k = 0;
        //                     }

        // There is always at least one set of weights
        positionalCoeffs[0] = Arrays.copyOfRange(regressionFitData.coefficients, 1 + (0*numKmerWindows), 1 + ((1+0)*numKmerWindows));
        positionalCoeffs[0] = ArrayTools.replace(positionalCoeffs[0], Double.NaN, 0);

        if (isRevCompPalindrome ||
            (strand == WeightMatrixTools.BindingStrand.POS) ||
            (strand == WeightMatrixTools.BindingStrand.NEG)) { // there is only one set of weights , copy them
            positionalCoeffs[1] = positionalCoeffs[0];
        }
        else { // populate the second set of weights
            positionalCoeffs[1] = Arrays.copyOfRange(regressionFitData.coefficients, 1 + (1*numKmerWindows), 1 + ((1+1)*numKmerWindows));
            positionalCoeffs[1] = ArrayTools.replace(positionalCoeffs[1], Double.NaN, 0);
        }

        // Now normalize all the weights by the maxCoeff
        double maxCoeff = MathTools.max(positionalCoeffs);
        positionalCoeffs = MathTools.multiply(1.0 / maxCoeff, positionalCoeffs);

        return(positionalCoeffs);

        //                     // look at both strands as default
        //                     int maxStrandIndex = 2;

        //                     // if motif is RevCompPalindrome then use only the positive strand
        //                     if (isRevCompPalindrome) {
        //                         maxStrandIndex = 1;
        //                     }

        //                     // populate the positionalWeights for both strands
        //                     // loop for the pos and neg strand
        //                     for (int strandIndex = 0; strandIndex < maxStrandIndex; strandIndex++) {

        //                         // loop for each startPosition
        //                         for (int posIndex=0; posIndex < numKmerWindows; posIndex++) {
        //                             positionalCoeffs[strandIndex] = Arrays.copyOfRange(regressionFitData.coefficients, 1 + (strandIndex*numKmerWindows), 1 + ((1+strandIndex)*numKmerWindows));
        //                             positionalCoeffs[strandIndex] = ArrayTools.replace(positionalCoeffs[strandIndex], Double.NaN, 0);
        //                             cleanedCoeffs[strandIndex] = cleanWeights(startPositions, positionalCoeffs[strandIndex], 28);
        //                             //cleanedCoeffs[strandIndex] = cleanWeights(startPositions, positionalCoeffs[strandIndex]);
        //                             smoothedCoeffs1[strandIndex] = MathTools.windowedTrimmedMean(cleanedCoeffs[strandIndex], 2, 1, 1);

        //                             //smoothedCoeffs[strandIndex] = MathTools.windowedMean(positionalCoeffs[strandIndex], 1);
        //                             //smoothedCoeffs1[strandIndex] = MathTools.windowedTrimmedMean(positionalCoeffs[strandIndex], 1, 1, 0);
        //                             //smoothedCoeffs2[strandIndex] = MathTools.windowedTrimmedMean(positionalCoeffs[strandIndex], 2, 3, 0);
        //                             //smoothedCoeffs3[strandIndex] = MathTools.windowedTrimmedMean(positionalCoeffs[strandIndex], 2, 2, 0);

        //                             positionalWeights[strandIndex] = smoothedCoeffs1[strandIndex];
        //                             //positionalWeights[strandIndex] = cleanedCoeffs[strandIndex];
        //                             //positionalWeights[strandIndex] = positionalCoeffs[strandIndex];
        //                         }
        //                     }

        //String[] datasetLabels = {"Positive Strand, ""Negative Strand"};
        //double[][] values = {MathTools.toDoubleArray(startPositions), positionalCoeffs[0], positionalCoeffs[1], smoothedCoeffs1[0], smoothedCoeffs1[1], smoothedCoeffs2[0], smoothedCoeffs2[1], smoothedCoeffs3[0], smoothedCoeffs3[1]};

        //String[] datasetLabels = {"Positive Strand Before Cleaning", "Negative Strand Before Cleaning", "Positive Strand After Cleaning", "Negative Strand After Cleaning", "Positive Strand After Smoothing", "Negative Strand After Smoothing"};
        //double[][] values = {MathTools.toDoubleArray(startPositions), positionalCoeffs[0], positionalCoeffs[1], cleanedCoeffs[0], cleanedCoeffs[1], smoothedCoeffs1[0], smoothedCoeffs1[1]};

    }

    public LinkedHashMap<SymbolList, WeightMatrix> getWeightMatrixMap(ArrayList<SymbolList> symbolListArrayList, double minAffinityWeight, Distribution background) {

        LinkedHashMap<SymbolList, WeightMatrix> weightMatrixMap = new LinkedHashMap<SymbolList, WeightMatrix>();
        try {
            for (SymbolList symbolList : symbolListArrayList ) {
                weightMatrixMap.put(symbolList, WeightMatrixTools.getWeightMatrix(symbolList, minAffinityWeight, background));
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

        return(weightMatrixMap);

    }

    public void initR() {

        String maxVsize = (String)initFileCache.get("REGRESSION", "MaxVsize");
        String maxMemSize = (String)initFileCache.get("REGRESSION", "MaxMemSize");

        ArrayList<String> rEngineArgs = new ArrayList<String>();
        rEngineArgs.add("--vanilla");
        rEngineArgs.add("--max-vsize="+maxVsize);
        rEngineArgs.add("--max-mem-size="+maxMemSize);

        REXP rexp;

        // Chech R version
        if (!Rengine.versionCheck()) {
            out.println("\nError: Version mismatch - Java files don't match library version.");
            System.exit(1);
        }
        out.println("\nCreating Rengine...");
        // 1) we pass the and command-line arguments to R
        // 2) we won't use the main loop at first, we'll start it later
        //    (that's the "false" as second argument)
        // 3) the callbacks are implemented by the TextConsole class above
        rengine = new Rengine(rEngineArgs.toArray(new String[0]), false, new TextConsole());
        //rengine = new Rengine(null, false, new TextConsole());
        out.println(" Done.");

        out.print("\nRengine successfully created, now waiting for R process...");
        // the engine creates R is a new thread, so we should wait until it's ready
        if (!rengine.waitForR()) {
            out.println(" Error: Cannot load R.");
            return;
        }
        out.println(" Done.");

        ArrayList<String> regressionStatements = new ArrayList<String>();
        regressionStatements.add("library(MASS)");
        regressionStatements.add("library(Matrix)");
        //regressionStatements.add("library(robustbase)");
        //regressionStatements.add("library(robust)");
        //regressionStatements.add("library(lars)");
        //regressionStatements.add("library(elasticnet)");
        regressionStatements.add("library(nnls)");
        //regressionStatements.add("library(minpack.lm)");
        regressionStatements.add("library(rJava)");
        regressionStatements.add(".jinit()");

        boolean displayReturnData = (boolean)initFileCache.getBoolean("REGRESSION", "DisplayReturnData");
        out.println("\nEvaluating the following regression statements in order in R:");
        for (String regressionStatement : regressionStatements) {
            out.println("\t"+regressionStatement);
            rexp = rengine.eval(regressionStatement);
            if (displayReturnData) {
                out.println("\t\treturn XT type is "+rexp.getType()+" : "+rexp.toString());
            }
        }
        out.println("Done.");
    }

    public void genWordCounts(String[] aKmerLengthsStringArray) {
        try {
            int numExperiments = reduceData.getNumExperiments();

            for (int experNum = 1; experNum <= numExperiments; experNum++) {

                //out.println("Hi7!");
                SymbolMatrix symbolMatrix = reduceData.createSymbolMatrixFile(experNum);
                if (symbolMatrix == null) {
                    //out.println("Hi8!");
                    symbolMatrix = reduceData.loadAllSymbolMatrixFile(experNum);
                }

                //out.println("Hi9!");
                for (String aKmerLengthString : aKmerLengthsStringArray) {
                    //out.println("Hi10!");
                    reduceData.createKmerMatrixFile(Integer.parseInt(aKmerLengthString), symbolMatrix, experNum);
                    //out.println("Hi11!");
                }
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    // get the weights for binding near the beginning of the probes
    public double[][] getRemovalWeights(int[] startPositions, double[][] weights) {
        int truncatePos = (int)initFileCache.getInt("AFFINITY MODELS", "TruncationPosition");

        //double[][] removalWeights = new double[2][startPositions.length];
        double[][] removalWeights = new double[2][];

        int firstStartPosition = startPositions[0];
        int startPosIsZeroIndex = Math.abs(firstStartPosition);

        for (int strandIndex = 0; strandIndex < 2; strandIndex++) {

            //removalWeights[strandIndex] = Arrays.copyOf(weights[strandIndex], weights[strandIndex].length);
            removalWeights[strandIndex] = new double[weights[strandIndex].length];
            Arrays.fill(removalWeights[strandIndex], 1.0);

            Arrays.fill(
                removalWeights[strandIndex],
                startPosIsZeroIndex+truncatePos,
                weights[strandIndex].length,
                0);

            // set first two positions to 0 as well;
            removalWeights[strandIndex][0] = 0;
            removalWeights[strandIndex][1] = 0;
            removalWeights[strandIndex][2] = 0;
            removalWeights[strandIndex][3] = 0.5;

        }

        // double max = MathTools.max(removalWeights);
        // removalWeights = MathTools.divide(removalWeights, max);

        return(removalWeights);
    }

    // get the weights for binding anywhere but near the beginning of the probes
    public double[][] getRemovedWeights(int[] startPositions, double[][] weights) {
        int truncatePos = (int)initFileCache.getInt("AFFINITY MODELS", "TruncationPosition");

        //double[][] removedWeights = new double[2][startPositions.length];
        double[][] removedWeights = new double[2][];

        int firstStartPosition = startPositions[0];
        int startPosIsZeroIndex = Math.abs(firstStartPosition);

        for (int strandIndex = 0; strandIndex < 2; strandIndex++) {

            removedWeights[strandIndex] = Arrays.copyOf(weights[strandIndex], weights[strandIndex].length);
            Arrays.fill(
                removedWeights[strandIndex],
                0,
                startPosIsZeroIndex+truncatePos,
                0);

        }

        return(removedWeights);
    }

    public double[][] noOverHangWeights(int[] startPositions, double[][] weights) {
        double[][] newWeights = new double[2][];
        int firstStartPosition = startPositions[0];
        int startPosIsZeroIndex = Math.abs(firstStartPosition);

        for (int strandIndex = 0; strandIndex < 2; strandIndex++) {
            newWeights[strandIndex] = Arrays.copyOfRange(
                weights[strandIndex],
                startPosIsZeroIndex,
                //weights[strandIndex].length - startPosIsZeroIndex);
                weights[strandIndex].length);
        }
        return(newWeights);
    }

    public int[] noOverHangStartPos(int[] startPositions) {
        int firstStartPosition = startPositions[0];
        int startPosIsZeroIndex = Math.abs(firstStartPosition);

        int[] newStartPos = Arrays.copyOfRange(
            startPositions,
            startPosIsZeroIndex,
            //startPositions.length - startPosIsZeroIndex);
            startPositions.length);

        return(newStartPos);
    }

    public double[][] splineWeights(int[] startPositions, double[][] rawWeights) {
        double[][] splineWeights = new double[2][];

        return(splineWeights);
    }

    // get the normalized L1 distance between the smoothed weights and the original weights
    public double getPosBiasWeightsL1(double[][] smoothWeights, double[][] rawWeights) {
        double L1 = 0;

        L1 += MathTools.L1Distance(smoothWeights[0], rawWeights[0]);
        L1 += MathTools.L1Distance(smoothWeights[1], rawWeights[1]);

        return(L1 / smoothWeights[0].length);
    }

    // cleaning up of the weights
    public double[][] cleanWeights(int[] startPositions, double[][] rawWeights) {
        double[][] cleanedWeights = new double[2][];

        cleanedWeights[0] = cleanWeights(startPositions, rawWeights[0]);
        cleanedWeights[1] = cleanWeights(startPositions, rawWeights[1]);

        double maxCoeff = MathTools.max(cleanedWeights);
        cleanedWeights = MathTools.multiply(1.0 / maxCoeff, cleanedWeights);

        return(cleanedWeights);
    }

    public double[][] smoothWeights(double[][] weights) {
        double[][] smoothedWeights = new double[2][];
        smoothedWeights[0] = MathTools.windowedTrimmedMean(weights[0], 2, 1, 1);
        smoothedWeights[1] = MathTools.windowedTrimmedMean(weights[1], 2, 1, 1);
        return(smoothedWeights);
    }

    public double[][] smoothWeights(int[] startPositions, double[][] weights) {
        //double[][] smoothedWeights = smoothWeights(weights);
        double[][] smoothedWeights = new double[2][startPositions.length];

        int firstStartPosition = startPositions[0];
        int startPosIsZeroIndex = Math.abs(firstStartPosition);

        for (int strandIndex = 0; strandIndex < 2; strandIndex++) {

            //double maxWeight = MathTools.max(weights[strandIndex]);
            //int maxIndex = ArrayTools.indexOf(weights[strandIndex], maxWeight);

//             int firstStartPosition = startPositions[0];
//             int startPosIsZeroIndex = Math.abs(firstStartPosition);

            // traverse backwards from startPosIsZeroIndex and copy all of the original values
            smoothedWeights[strandIndex] = ArrayTools.replaceRange(
                weights[strandIndex],
                0,
                startPosIsZeroIndex+1,
                smoothedWeights[strandIndex],
                0);


            // smooth the RHS of the startPosIsZeroIndex, extension=2, drop highest and lowest
            double[] rangeToSmooth = Arrays.copyOfRange(weights[strandIndex], startPosIsZeroIndex, startPositions.length);

            //double[] rangeSmoothed = MathTools.windowedTrimmedMean(rangeToSmooth, 4, 1, 1);
            //double[] rangeSmoothed = MathTools.windowedTrimmedMean(rangeToSmooth, 2, 1, 1);
            //double[] rangeSmoothed = MathTools.windowedTrimmedMean(rangeToSmooth, 1, 1, 0);
            //double[] rangeSmoothed = MathTools.windowedTrimmedMean(rangeToSmooth, 1, 0, 1);
            double[] rangeSmoothed = MathTools.windowedMean(rangeToSmooth, 1);

            smoothedWeights[strandIndex] = ArrayTools.replaceRange(
                    rangeSmoothed,
                    0,
                    rangeSmoothed.length,
                    smoothedWeights[strandIndex],
                    startPosIsZeroIndex);

//             // traverse backwards from startPos=0 and copy all of the original values
//             smoothedWeights[strandIndex] = ArrayTools.replaceRange(
//                 weights[strandIndex],
//                 0,
//                 startPosIsZeroIndex+1,
//                 smoothedWeights[strandIndex],
//                 0);

//             // traverse forwards from startPos=0 and set to 0 after first non-positive weights[posIndex]
//             int firstZeroValueOnRight = ArrayTools.indexOf(weights[strandIndex], 0, startPosIsZeroIndex);

//             if (firstZeroValueOnRight != -1) {
//                 smoothedWeights[strandIndex] = ArrayTools.replaceRange(
//                     weights[strandIndex],
//                     firstZeroValueOnRight,
//                     weights[strandIndex].length,
//                     smoothedWeights[strandIndex],
//                     firstZeroValueOnRight);
//             }

//             boolean setToZero = false;
//             for (int posIndex = startPosIsZeroIndex; posIndex < weights.length; posIndex++) {
//                 if (setToZero) {
//                     smoothedWeights[strandIndex][posIndex] = 0;
//                 }
//                 else if ((weights[strandIndex][posIndex] <= 0) || (Double.isNaN(weights[strandIndex][posIndex]))) {
//                     smoothedWeights[strandIndex][posIndex] = 0;
//                     setToZero = true;
//                 }
//             }

//             double maxWeight = MathTools.max(weights);
//             int maxIndex = ArrayTools.indexOf(weights[strandIndex], maxWeight);
//             smoothedWeights[strandIndex][maxIndex] = maxWeight;

//             for (int posIndex = 0; posIndex < weights.length; posIndex++) {
//                 if (weights[strandIndex][posIndex] == 0) {
//                     smoothedWeights[strandIndex][posIndex] = 0;
//                 }
//             }
        }

        return(smoothedWeights);
    }

    // performs in place cleaning up of the weights
    public double[] cleanWeights(int[] startPositions, double[] weights) {
        double[] cleanedWeights = Arrays.copyOf(weights, weights.length);
        int firstStartPosition = startPositions[0];
        int startPosIsZeroIndex = Math.abs(firstStartPosition);

        // first traverse backwards from startPos=0 and set to 0 after non-positive weight
        boolean setToZero = false;
        for (int posIndex = startPosIsZeroIndex; posIndex >= 0; posIndex--) {
            if (setToZero) {
                cleanedWeights[posIndex] = 0;
            }
            else if ((weights[posIndex] <= 0) || (Double.isNaN(weights[posIndex]))) {
                cleanedWeights[posIndex] = 0;
                setToZero = true;
            }
        }

        // second traverse forwards from (length - startPos=0) and set to 0 after first non-positive weight
        setToZero = false;
        for (int posIndex = (weights.length - startPosIsZeroIndex); posIndex < weights.length; posIndex++) {
            if (setToZero) {
                cleanedWeights[posIndex] = 0;
            }
            else if ((weights[posIndex] <= 0) || (Double.isNaN(weights[posIndex]))) {
                cleanedWeights[posIndex] = 0;
                setToZero = true;
            }
        }

        return(cleanedWeights);
    }

//     // performs in place cleaning up of the weights
//     public double[] cleanWeights(int[] startPositions, double[] weights, int maxPos) {
//         double[] cleanedWeights = Arrays.copyOf(weights, weights.length);
//         int firstStartPosition = startPositions[0];
//         int startPosIsZeroIndex = Math.abs(firstStartPosition);

//         // first traverse backwards from startPos=0 and set to 0 after first negative weight
//         boolean setToZero = false;
//         for (int posIndex = startPosIsZeroIndex; posIndex >= 0; posIndex--) {
//             if (setToZero) {
//                 cleanedWeights[posIndex] = 0;
//             }
//             else if ((weights[posIndex] <= 0) || (Double.isNaN(weights[posIndex]))) {
//                 cleanedWeights[posIndex] = 0;
//                 setToZero = true;
//             }
//         }

//         // second traverse forwards from startPos=0 and set to 0 after first negative weight
//         setToZero = false;
//         for (int posIndex = startPosIsZeroIndex; posIndex < weights.length; posIndex++) {
//             if (setToZero) {
//                 cleanedWeights[posIndex] = 0;
//             }
//             else if ((weights[posIndex] <= 0) || (Double.isNaN(weights[posIndex]))) {
//                 cleanedWeights[posIndex] = 0;
//                 setToZero = true;
//             }
//             else if (startPositions[posIndex] > maxPos) {
//                 cleanedWeights[posIndex] = 0;
//                 setToZero = true;
//             }
//         }

//         return(cleanedWeights);
//     }

    public void genIntensities(String intensitiesPathName) {
        try {
            //String intensitiesDir = (String)initFileCache.get("DATA", "ProbeIntensitiesDir");
            //String ProbeLocsDir = (String)initFileCache.get("DATA", "ProbeLocsDir");
            //String resultsDir = (String)initFileCache.get("OUTPUT", "Directory");
            //String motifsDir = intensitiesDir + File.separator + resultsDir;

            // get the strand from the init file
            WeightMatrixTools.BindingStrand strand = reduceData.getStrand();
            WeightMatrixTools.BothStrandsCalc calc = reduceData.getBothStrandsCalc();

            out.print("\nCalculating the counts...");
            //double[] motifCounts = wordCounts.getCountsArray(psam, strand);

            // Just use the first experiment for now
            double[] motifCounts = reduceData.getCountsArray(psam, null, null, null, strand, calc, 1, -1);
            out.println(" Done.\n");

//             out.print("\nCalculating the counts on the positive strand...");
//             double[] posMotifCounts = wordCounts.getCountsArray(psam, WeightMatrixTools.BindingStrand.POS);
//             out.println(" Done.\n");

//             out.print("\nCalculating the counts on the negative strand...");
//             double[] negMotifCounts = wordCounts.getCountsArray(psam, WeightMatrixTools.BindingStrand.NEG);
//             out.println(" Done.\n");

//             double[] motifCounts = MathTools.add(posMotifCounts, negMotifCounts);

            /////////////////////////////////////////////////////////////////////////
            //
            // make modifications to the counts here if desired
            //   - like adding noise
            //
            motifCounts = MathTools.multiply(1000, motifCounts);

            /////////////////////////////////////////////////////////////////////////

            // Just use the first experiment for now
            String[] keys = reduceData.getUsedProbeIDs(1, -1).toArray(new String[0]);

            out.print("\nWriting generated intensities to "+intensitiesPathName+"...");
            StringBuffer outBuffer = new StringBuffer();
            for (int i=0; i < keys.length; i++) {
                outBuffer.append(keys[i] + "\t" + motifCounts[i]+ "\n");
            }
            FileTools.write(outBuffer.toString(), intensitiesPathName, false);
            out.println(" Done.\n");

        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public GridBagLayoutFrame createLogoFrame() {
        GridBagLayoutFrame logoFrame = new GridBagLayoutFrame("Regression Results for different strands", true);
        logoFrame.add(new JLabel("Regression", SwingConstants.CENTER), 10, 10, false, false);
        logoFrame.add(new JLabel("Seed", SwingConstants.CENTER), 10, 10, false, false);
        logoFrame.add(new JLabel("Strand", SwingConstants.CENTER), 10, 10, false, false);
        logoFrame.add(new JLabel("Logo", SwingConstants.CENTER), 10, 10, false, false);
        logoFrame.add(new JLabel("R-squared", SwingConstants.CENTER), 10, 10, false, false);
        logoFrame.add(new JLabel("Occupancy Bias", SwingConstants.CENTER), 10, 10, false, false);
        logoFrame.add(new JLabel("Bias R-squared", SwingConstants.CENTER), 10, 10, false, true);

        return(logoFrame);
    }

    public GridBagLayoutFrame createFSAMFrame(String seed, String strand) {
        GridBagLayoutFrame logoFrame = new GridBagLayoutFrame("FSAM Logo - Regressions performed on "+strand+" strands with seed "+seed, true);
        //logoFrame.add(new JLabel("FSAM Logo"), 10, 10, false, true);
        //logoFrame.add(new JLabel("R-squared"), 10, 10, false, true);

        return(logoFrame);
    }

    public void addMotif(GridBagLayoutFrame logoFrame, WeightMatrix aPWM, String regressionTitle, String seed, String strand, double rSquared, BioJavaChart chart, double posWeightsRSquared) {
        WeightMatrixLogo logoPanel = new WeightMatrixLogo(aPWM);

        logoFrame.add(new JLabel(regressionTitle, SwingConstants.CENTER), 10, 10, false, false);
        logoFrame.add(new JLabel(seed, SwingConstants.CENTER), 10, 10, false, false);
        logoFrame.add(new JLabel(strand, SwingConstants.CENTER), 10, 10, false, false);
        logoFrame.add(logoPanel, 10, 10, false, false);

        //logoFrame.add(new JLabel(Double.toString(rSquared), SwingConstants.CENTER), 10, 10, false, false);
        logoFrame.add(new JLabel(format(rSquared, 4, 4), SwingConstants.CENTER), 10, 10, false, false);

        logoFrame.add(chart, 10, 10, false, false);

        //logoFrame.add(new JLabel(Double.toString(posWeightsRSquared), SwingConstants.CENTER), 10, 10, false, true);
        logoFrame.add(new JLabel(format(posWeightsRSquared, 4, 4), SwingConstants.CENTER), 10, 10, false, true);

        return;
    }

    //public void addMotif(GridBagLayoutFrame logoFrame, WeightMatrix aPWM, double rSquared) {
    public void addMotif(GridBagLayoutFrame logoFrame, WeightMatrix aPWM) {

        ArrayList<WeightMatrixFeature> features = new ArrayList<WeightMatrixFeature>();

        Symbol[] modsArray = new Symbol[(aPWM.columns()*2)-1];
        modsArray[0] = DNATools.g();
        modsArray[2] = DNATools.a();
        WeightMatrixFeature aFeature = new WeightMatrixFeature("dependency", aPWM, null, false, false, 0.2381, modsArray, complementTable);
        features.add(aFeature);

        modsArray = new Symbol[(aPWM.columns()*2)-1];
        modsArray[0] = DNATools.a();
        modsArray[2] = DNATools.t();
        aFeature = new WeightMatrixFeature("dependency", aPWM, null, true, false, 0.13, modsArray, complementTable);
        features.add(aFeature);

        modsArray = new Symbol[(aPWM.columns()*2)-1];
        modsArray[2] = DNATools.c();
        modsArray[4] = DNATools.n();
        modsArray[6] = DNATools.c();
        aFeature = new WeightMatrixFeature("deletion", aPWM, null, false, false, 0.08, modsArray, complementTable);
        features.add(aFeature);

        modsArray = new Symbol[(aPWM.columns()*2)-1];
        modsArray[2] = DNATools.t();
        modsArray[3] = DNATools.n();
        modsArray[4] = DNATools.g();
        aFeature = new WeightMatrixFeature("insertion", aPWM, null, true, false, 0.07, modsArray, complementTable);
        features.add(aFeature);

        modsArray = new Symbol[(aPWM.columns()*2)-1];
        modsArray[2] = DNATools.g();
        modsArray[3] = DNATools.g();
        modsArray[4] = DNATools.g();
        aFeature = new WeightMatrixFeature("insertion", aPWM, null, true, false, 0.05, modsArray, complementTable);
        features.add(aFeature);

        //WeightMatrixLogo logoPanel = new WeightMatrixLogo(aPWM);
        WeightMatrixLogo logoPanel = new WeightMatrixLogo(aPWM, true, true, true, 0, WeightMatrixLogo.PositionLoc.TOP, features, false, .0001,0);

        logoFrame.add(logoPanel, 10, 10, false, true);
        //logoFrame.add(new JLabel(Double.toString(rSquared)), 10, 10, false, true);

        return;
    }

    public static boolean isSelfReverseComplement(Object aMotif, ReversibleTranslationTable aComplementTable) {
        try {
            if (aMotif instanceof WeightMatrix) {
                return(WeightMatrixTools.isSelfReverseComplement((WeightMatrix)aMotif, aComplementTable));
            }
            else if (aMotif instanceof SymbolList) {
                return(SymbolListTools.isSelfReverseComplement((SymbolList)aMotif, aComplementTable));
            }
            else if (aMotif instanceof FeaturedWeightMatrix) {
                return(WeightMatrixTools.isSelfReverseComplement(((FeaturedWeightMatrix)aMotif).getPosStrandPWM(), aComplementTable));
            }
        }
        catch(Exception ex) {
            ex.printStackTrace(System.err);
        }
        return(false);
    }

    public RegressionFitData getRegressionFit(
        FeaturedWeightMatrix anFSAM,
        String linearAlg,
        int unifiedColumn) {

        return(getRegressionFit(
                anFSAM,
                anFSAM.getPositionalWeights(),
                anFSAM.getName(),
                linearAlg,
                anFSAM.getStrand(),
                anFSAM.getStrand(),
                anFSAM.getCalc(),
                unifiedColumn));
    }

    public RegressionFitData getRegressionFit(
        Object motif,
        double[][] positionalWeights,
        String label,
        String linearAlg,
        WeightMatrixTools.BindingStrand effectiveStrand,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        int unifiedColumn) {

        LinkedHashMap<Object, LinkedHashSet<SymbolList>> motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();
        motifToHammingSphereMap.put(motif, null);

        RegressionFitData regressionFitData = fitModel(
            motifToHammingSphereMap,
            null, //startPositions
            positionalWeights,
            null, // aMotifsToModsMap
            null, // aMotifsToMandatoriesMap
            label,
            linearAlg,
            effectiveStrand,
            strand,
            calc,
            null, //startParams
            null, //refSymList
            Double.NaN, //refKaReal
            unifiedColumn);

        // coefficients are as returned in R
        // First coefficient is the intercept, then the rest are in independent-var order
        out.print(regressionFitData.toString());

        return(regressionFitData);
    }

    public double getRsquared(
        Object motif,
        double[][] positionalWeights,
        String label,
        String linearAlg,
        WeightMatrixTools.BindingStrand effectiveStrand,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        int unifiedColumn) {

        RegressionFitData regressionFitData = getRegressionFit(
            motif,
            positionalWeights,
            label,
            linearAlg,
            effectiveStrand,
            strand,
            calc,
            unifiedColumn);

        return(regressionFitData.rSquared);
    }

    // The positional weights should be the same as those used to determine the sequence-specific PSAM (Ka's)
    public double getRelativeKaNonSpec(
        double intercept,
        double motifCoeff,
        double[][] positionalWeights,
        double eToMu,
        int probeSeqLengths,
        int motifLength,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double revCompSimilarity) {

        double relKaNonSpec;
        double relOccupancySum = 0;

        ///////////////////////////////////////////////
        ///////////////////////////////////////////////
        if (true) {
            return(0);
        }
        ///////////////////////////////////////////////
        ///////////////////////////////////////////////

        // relOccupancySum = numKmerWindows;
        if (positionalWeights == null) {
            // calculate the number of kmer windows
            int numKmerWindows = probeSeqLengths - motifLength + 1;
            switch (strand) {
            case POS:
            case NEG:
                relOccupancySum = numKmerWindows;
                break;
            case BOTH:
                switch (calc) {
                case MAX:
                    relOccupancySum = numKmerWindows;
                    break;
                case SUM:
                    relOccupancySum = numKmerWindows*2;
                    break;
                case UNION:
                    relOccupancySum = (numKmerWindows*2) - (numKmerWindows*revCompSimilarity);
                    break;
                case NORMED_SUM:
                    relOccupancySum = (numKmerWindows*2) / (revCompSimilarity + 1);
                    break;
                }
                break;
            }
        }

        // relOccupancySum = sum(positionalWeights);
        else {
            switch (strand) {
            case POS:
                relOccupancySum = MathTools.sum(positionalWeights[0]);
                break;
            case NEG:
                relOccupancySum = MathTools.sum(positionalWeights[1]);
                break;
            case BOTH:
                switch (calc) {
                case MAX:
                    for (int i=0; i < positionalWeights[0].length; i++) {
                        relOccupancySum += Math.max(positionalWeights[0][i], positionalWeights[1][i]);
                    }
                    break;
                case SUM:
                    for (int i=0; i < positionalWeights[0].length; i++) {
                        relOccupancySum += positionalWeights[0][i] + positionalWeights[1][i];
                    }
                    break;
                case UNION:
                    for (int i=0; i < positionalWeights[0].length; i++) {
                        relOccupancySum += positionalWeights[0][i] + positionalWeights[1][i] - (positionalWeights[0][i] * positionalWeights[1][i]);
                    }
                    break;
                case NORMED_SUM:
                    for (int i=0; i < positionalWeights[0].length; i++) {
                        relOccupancySum += (positionalWeights[0][i] + positionalWeights[1][i]) / (revCompSimilarity+1);
                    }
                    break;
                }
                break;
            }
        }

        if (eToMu == 0) {
            relKaNonSpec = (intercept / relOccupancySum) / motifCoeff;
        }
        else {
            double normalizedInt = intercept / relOccupancySum;
            double seedKaSeqSpec = motifCoeff / (eToMu * (1 - motifCoeff));
            double aKaNonSpec = normalizedInt / (eToMu * (1 - normalizedInt));
            relKaNonSpec = aKaNonSpec / seedKaSeqSpec;
        }

        out.println("\nRelative Ka non-specific is "+relKaNonSpec);
        return(relKaNonSpec);
    }

    public int[] getStartPositions(int probeSeqLengths, int motifLength) {
        int numKmerWindows = probeSeqLengths + motifLength - 1;
        int[] startPositions = new int[numKmerWindows];

        int firstStartPos = -1 * (motifLength - 1);
        for (int posIndex= 0; posIndex < numKmerWindows; posIndex++) {
            startPositions[posIndex] = firstStartPos + posIndex;
        }
        return(startPositions);
    }

    public void writeRelAffinitySphere(
        FeaturedWeightMatrix fsam,
        double featureThreshold,
        double relAffinityThreshold,
        String affinitySphereFileName)
    {
        System.out.print("\n"+"Creating the Relative Affinity Sphere that contains all sequences with relative affinity >= "+relAffinityThreshold+"...");

//         LinkedHashMap<SymbolList, Double> relAffinitySphere = fsam.getRelAffinitySphere(
//             //ScoreType.PROBABILITY,
//             null,
//             featureThreshold,
//             relAffinityThreshold,
//             DNATools.complementTable());

//         Table affinitySphereTable = new Table(relAffinitySphere);
//         affinitySphereTable.sort(1);
//         affinitySphereTable.reverse();

        Table affinitySphereTable = fsam.getRelAffinitySphere(
            //ScoreType.PROBABILITY,
            null,
            featureThreshold,
            relAffinityThreshold,
            DNATools.complementTable());

        System.out.println( " Done.");

        //FileTools.write(affinitySphereTable.toString("\t"), affinitySphereFileName, false);
        System.out.print("\n"+"Writing Relative Affinity Sphere to file "+affinitySphereFileName+"...");
        StringBuffer stringBuffer = new StringBuffer();
        stringBuffer.append("Kmer" +"\t"+ "relAffinity");

        for (Iterator rowIter = affinitySphereTable.iterator(); rowIter.hasNext(); ){
            java.util.List currentRow = (java.util.List) rowIter.next();
            String kmer = ((SymbolList)currentRow.get(0)).seqString().toUpperCase();
            Double relAffinity = (Double)currentRow.get(1);
            stringBuffer.append("\n"+kmer +"\t"+ relAffinity);
        }
        FileTools.write(stringBuffer.toString(), affinitySphereFileName, false);
        System.out.println( " Done.");


    }

    // Loop through the ids calling fitModels on each experiment
    public FeaturedWeightMatrix fitModels(String motif,
        String[] seeds,
        WeightMatrix psam,
        boolean displayMotifs,
        SymbolList refSymList,
        double referenceKa,
        String[] ids,
        String seedType,
        boolean saveAll,
        String[] kmerLengthsArray,
        String dreamInputFile) {

        FeaturedWeightMatrix fsam = null;

        try {

            if (ids == null) {
                // if ids is null then construct a list of all of them
                // 		String[] lineStrings = FileTools.readStrings(dreamInputFile, 0);
                // 		ArrayList<String> idsArrayList = StringTools.getColumn(lineStrings, tabPattern, null, initFileCache.getInt("DREAM", "LabelColumn"));
                // 		LinkedHashSet<String> idsHashSet = new LinkedHashSet<String>(idsArrayList);

                LinkedHashSet<String> idsHashSet = FileTools.getColumnUnique(dreamInputFile,
                    1,
                    (String)null,
                    -1,
                    tab,
                    initFileCache.getInt("DREAM", "LabelColumn"));

                ids = (String[])idsHashSet.toArray(new String[0]);
            }

            // loop through all the ids and generate FSAMs
            for (int j=0; j < ids.length; j++) {

                // RESET ANY GLOBALS
                removeProbesArray = null;

                // create new SymbolMatrix and KmerMatrix
                initFileCache.put("DREAM", "ID", ids[j]);
                reduceData = new ReduceData(initFileCache);

                // create new symbol matrix and kmer-table
                //out.println("Hi1!");
                genWordCounts(kmerLengthsArray);

                String dreamLabel = (String)initFileCache.get("DREAM", "Label");

                SymbolList seedSymList = null;
                if (seeds != null) {
                    //out.println("Hi5!");
                    seedSymList = seedSymList = DNATools.createDNA(seeds[j]);
                }

                //out.println("Hi6!");
                // create the fsam from the top pearson-correlation 6mer
                fsam = fitModels(motif,
                    seedSymList,
                    //null,
                    psam,
                    displayMotifs,
                    null, // refSymList
                    Double.NaN, // refKa
                    dreamLabel+"."+ids[j], // proteinLabel
                    //FileTools.stripExtension(anInitFile),
                    -1, // unifiedColumn
                    saveAll,
                    seedType);

            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

        //out.println("Hi4!");
        return(fsam);
    }

    // Loop through the unifiedMotifs calling fitModels on each motif
    public FeaturedWeightMatrix fitModels(String motif, String[] seeds, WeightMatrix psam, boolean displayMotifs, SymbolList refSymList, double referenceKa, String[] unifiedMotifs, boolean saveAll, String seedType) {
        FeaturedWeightMatrix fsam = null;

        try {
            int experimentNumber = 1;

            String inputFileType = (String)initFileCache.get("INTENSITIES", "InputFileType");
            String unifiedFile = (String)initFileCache.get("SEQUENCE", Integer.toString(experimentNumber));
            String[] columnLabels = FileTools.columnLabels(unifiedFile, null, "\t");

            // get rid of newline
            columnLabels[columnLabels.length - 1] = columnLabels[columnLabels.length - 1].replace("\n", "");

            int startColumn = -1;
            out.println("Input File Type is :"+inputFileType);
            if (inputFileType.equalsIgnoreCase(".unified")) {
                // sequence   intensity1  intensity2 ....
                // 0          1           2
                startColumn = 1;
            }
            else if (inputFileType.equalsIgnoreCase(".unifiedLabels")) {
                // label   sequence   intensity1  intensity2 ....
                // 0          1           2
                startColumn = 2;
            }
            else {
                out.println("Error: Unified table format unknown!!");
            }

            //out.println("Hi1!");

            // if no explicit motifs specified then loop through them all
            if ((unifiedMotifs == null) || (unifiedMotifs.length == 0)) {
                for (int i = startColumn; i < columnLabels.length; i++) {

                    SymbolList seedSymList = null;
                    if (seeds != null) {
                        seedSymList = seedSymList = DNATools.createDNA(seeds[i-2]);
                    }




                    String intensitiesDir = (String)initFileCache.get("DATA", "ProbeIntensitiesDir");
                    //String probeLocsDir = (String)initFileCache.get("DATA", "ProbeLocsDir");
                    String outputDir = (String)initFileCache.get("OUTPUT", "Directory");
                    String resultsDir = intensitiesDir + File.separator + outputDir;


                    // load initPositionalWeights
                    String positionalBiasFileName = resultsDir + File.separator + columnLabels[i]+".psam."+10+"nt.positionalBias.table";
                    //if (positionalBiasFileName != null) {
                    if (FileTools.exists(positionalBiasFileName)) {
                        out.print("\nReading the positional-binding bias profile from "+positionalBiasFileName+"...");
                        String[]  biasFileRows = FileTools.readStrings(positionalBiasFileName, 0);
                        String[][] biasFileMatrix = StringTools.toStringMatrix(biasFileRows, null, true);
                        this.initPositionalBias = StringTools.toDoubleMatrix(biasFileMatrix, true, true);
                        out.println(" Done.\n");
                        out.println("\nThe loaded positional-binding bias profile is..\n"+StringTools.toString(biasFileMatrix));
                    }
                    else {
                        this.initPositionalBias = null;
                    }


                    String psamPathName = resultsDir + File.separator + columnLabels[i] +".psam."+10+"nt.xml";
                    //if (psamPathName != null) {
                    if (FileTools.exists(psamPathName)) {
                        // Load the PSAM from file
                        if (psamPathName.endsWith(".xml") || psamPathName.endsWith(".XML")) {
                            out.print("\nReading the PSAM XML file from "+psamPathName+"...");
                            psam = WeightMatrixTools.readFromXML(psamPathName);
                        }
                        else {
                            out.print("\nReading the PSAM MatrixReduce XML file from "+psamPathName+"...");
                            psam = WeightMatrixTools.readFromMatrixReduce(psamPathName);
                        }
                        psam.setName("PSAM");
                        out.println(" Done.\n");
                    }
                    else {
                        psam = null;
                    }



                    //String proteinLabel = columnLabels[i].replace("\\n", "");
                    //fsam = fitModels(motif, seedSymList, psam, displayMotifs, refSymList, referenceKa, proteinLabel, i, saveAll);
                    fsam = fitModels(motif, seedSymList, psam, displayMotifs, refSymList, referenceKa, columnLabels[i], i, saveAll, seedType);
                    //out.println("Hi2!");
                }
            }
            // generate motifs only for the columns specified
            else {
                for (int j=0; j < unifiedMotifs.length; j++) {
                    for (int i = startColumn; i < columnLabels.length; i++) {
                        //out.println("unifiedMotif="+unifiedMotifs[j]+"  columnLabel="+columnLabels[i]);
                        if (unifiedMotifs[j].equals(columnLabels[i])) {
                            //out.println("Hi3!");

                            SymbolList seedSymList = null;
                            if (seeds != null) {
                                //out.println("Hi5!");
                                seedSymList = seedSymList = DNATools.createDNA(seeds[j]);
                            }

                            fsam = fitModels(motif, seedSymList, psam, displayMotifs, refSymList, referenceKa, unifiedMotifs[j], i, saveAll, seedType);
                            break;
                        }
                    }
                }
            }
        }

        catch (Exception ex) {
            ex.printStackTrace();
        }

        //out.println("Hi4!");
        return(fsam);
    }


    public double[] fitPsams(boolean displayMotifs, String[] kmerLengthsArray, int unifiedColumn) {

        LinkedHashMap<Object, boolean[][]> motifsToMandatoriesMap;
        RegressionFitData regressionFitData = null;
        LinkedHashMap<Object, LinkedHashSet<SymbolList>> motifToHammingSphereMap = null;
        double[][] positionalWeights = initPositionalBias;

        reduceData = new ReduceData(initFileCache);

        // init the R engine
        if (rengine == null) {
            initR();
        }

        // create new symbol matrix and kmer-table
        genWordCounts(kmerLengthsArray);

        String intensitiesDir = (String)initFileCache.get("DATA", "ProbeIntensitiesDir");
        //String probeLocsDir = (String)initFileCache.get("DATA", "ProbeLocsDir");
        String outputDir = (String)initFileCache.get("OUTPUT", "Directory");
        String resultsDir = intensitiesDir + File.separator + outputDir;

        // get the strand from the init file
        WeightMatrixTools.BindingStrand effectiveStrand = reduceData.getStrand();
        WeightMatrixTools.BindingStrand strand = effectiveStrand;
        boolean isRevCompPalindrome = false;
        double featuresScalar = 1.0;
        double eToMu = 0;

        MutableDouble alpha = new MutableDouble(1.0);
        MutableDouble beta = new MutableDouble(1.0);
        MutableDouble nonSpecKa = new MutableDouble(Math.pow(10, -12));

        WeightMatrixTools.BothStrandsCalc calc = reduceData.getBothStrandsCalc();

        // Get the motif length from the motif (since the seedSymList may be null)
        String[] correspondSymbols = null;

        char[] kmerModelLengths = initFileCache.getChars("AFFINITY MODELS", "KmerModelLengths");
        String kmerModelPosBiases = (String) initFileCache.get("AFFINITY MODELS", "KmerModelPosBiases");

        double[][][] kmerModelWeights = null;
        boolean[] isLearnedPosBias = null;
        boolean[] isUniformPosBias = null;

        // the normalized L1 distance between the raw and smooth positional bias weights
        double posBiasWeightsL1 = 0;

        double rSquared = 0;
        double noPosWeightsRSquared = -1;

        String spacerLenRegressionAlgorithm = (String)initFileCache.get("REGRESSION", "spacerLenRegressionAlgorithm");

        motifsToMandatoriesMap = new LinkedHashMap<Object, boolean[][]>();
        motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();

        //                     if (useHammingDistForApprox) {
        //                         motifToHammingSphereMap.put(hamming1PSAM, hammingSpheres.get(hammingDistForApprox));
        //                     }
        //                     else {
        //                         motifToHammingSphereMap.put(hamming1PSAM, null);
        //                     }

        //Distribution allOnesDist = HMMTools.createUniformDistribution(alphabet,true);
        Distribution allOnesDist = new SimpleDistribution(alphabet);
        DistributionTools.allWeightsTo(allOnesDist, 1.0);

        // Add all spacer-lengths between Min and Max inclusive
        // 0-based indexing
        int spacerLenMin = (int)initFileCache.getInt("AFFINITY MODELS", "SpacerLengthMin");
        int spacerLenMax = (int)initFileCache.getInt("AFFINITY MODELS", "SpacerLengthMax");

        boolean getSpacerLengths = ((String)initFileCache.get("AFFINITY MODELS", "GetSpacerLengthAffinities")).equalsIgnoreCase("Yes");

        double[][][] newPsamsPosWeights = new double[psamsArray.length * (spacerLenMax - spacerLenMin + 2)][][];

        // loop for each psam
        for (int varNum = 0; varNum < psamsArray.length; varNum++) {

            WeightMatrix aPsam = psamsArray[varNum];
            //double[][] aPositionalWeights = psamsPosWeights[varNum];

            aPsam.setName("PSAM_"+varNum);

            // Don't send in a hamming sphere ever
            motifToHammingSphereMap.put(aPsam, null);

            int aMotifLength = aPsam.columns();

            if (getSpacerLengths) {

                int spacerPos = (int)Math.round(Math.floor((float) aMotifLength / 2));

                newPsamsPosWeights[varNum * (spacerLenMax - spacerLenMin + 2)] = psamsPosWeights[varNum];

                for (int spacerLen = spacerLenMin; spacerLen <= spacerLenMax; spacerLen++) {

                    newPsamsPosWeights[(varNum * (spacerLenMax - spacerLenMin + 2)) + (spacerLen - spacerLenMin +1)] = psamsPosWeights[varNum];

                    boolean[][] mandatoryColumns = new boolean[2][aMotifLength+spacerLen];
                    // include the positions neighboring the spacer as mandatory
                    Arrays.fill(mandatoryColumns[0], spacerPos-1, spacerPos+spacerLen+2, true);
                    mandatoryColumns[1] = ArrayTools.reverse(mandatoryColumns[0]);

                    // create uniform-dist spacer in the PSAM
                    WeightMatrix addedSpacerPSAM = WeightMatrixTools.insert(aPsam, spacerPos, spacerLen, allOnesDist);

                    // name has 1-based indexing
                    //addedSpacerPSAM.setName("[I("+ spacerPos +"-"+(spacerPos+1)+","+spacerLen+","+"N"+")]");
                    //addedSpacerPSAM.setName("insert_"+ spacerPos +"_"+(spacerPos+1)+"_"+spacerLen+"_"+"N"+"");
                    addedSpacerPSAM.setName("PSAM_"+varNum+"_insert_"+ spacerPos +"_"+(spacerPos+1)+"_"+spacerLen+"_"+"N"+"");

                    //                         out.println((byte)0, "**********************************************************");
                    //                         out.println((byte)0, "WeightMatrix: "+WeightMatrixTools.toString(addedSpacerPSAM));
                    //                         out.println((byte)0, "**********************************************************");

                    // Don't send in a hamming sphere ever
                    motifsToMandatoriesMap.put(addedSpacerPSAM, mandatoryColumns);
                    motifToHammingSphereMap.put(addedSpacerPSAM, null);

                } // end for (each spacer)

            } // if getSpacerLengths

        } // end for (each motif)

        if (getSpacerLengths) {
            psamsPosWeights = newPsamsPosWeights;
        }

        // Add Positions to the RHS
        //this.addWeightsToLeft = false;
        this.subtractOtherAffinities = true;

        // perform linear regression to find coefficients for each spacer-length
        regressionFitData = fitModel(
            motifToHammingSphereMap,
            null, //startPositions
            positionalWeights,
            null, //motifsToModsMap
            motifsToMandatoriesMap,
            "nucleotide-independent spacer-length features",
            //linearAlgorithm,
            spacerLenRegressionAlgorithm,
            effectiveStrand,
            strand,
            calc,
            null, //startParams
            null, //refSymList
            Double.NaN, //refKaReal);
            unifiedColumn);

        this.addWeightsToLeft = true;
        this.subtractOtherAffinities = false;

        // coefficients are as returned in R
        // First coefficient is the intercept, then the rest are in independent-var order
        out.print(regressionFitData.toString());

        double[] coefs = regressionFitData.getCoefficients();
        double[] coefsWithoutIntercept = Arrays.copyOfRange(coefs, 1, coefs.length);
        double maxCoef = MathTools.max(coefsWithoutIntercept);
        double[] relAffinities = MathTools.divide(coefsWithoutIntercept, maxCoef);

        out.println("\nRelative Affiniities are: "+StringTools.toString(relAffinities, ", "));

        return(relAffinities);
    }



    // Fits regression models according to the INI file and builds as FSAM
    public FeaturedWeightMatrix fitModels(String motif, SymbolList seedSymList, WeightMatrix psam, boolean displayMotifs, SymbolList refSymList, double referenceKa, String proteinLabel, int unifiedColumn, boolean saveAll, String seedType) {
        try {

            boolean cacheCounts = (boolean)initFileCache.getBoolean("REGRESSION", "CacheCounts");
            boolean useHammingDistForApprox = (boolean)initFileCache.getBoolean("REGRESSION", "UseHammingDistForApprox");
            int hammingDistForApprox = (int)initFileCache.getInt("REGRESSION", "HammingDistForApprox");
            double minAffinityWeight = initFileCache.getDouble("REGRESSION", "MinAffinityWeight");

            boolean includeRevComps = true;

            String linearAlgorithm = (String)initFileCache.get("REGRESSION", "LinearAlgorithm");
            String nonlinearAlgorithm = (String)initFileCache.get("REGRESSION", "NonlinearAlgorithm");

            //boolean includeIntercept = (boolean)initFileCache.getBoolean("REGRESSION", "IncludeIntercept");
            String featureDiscoveryAlgorithm = (String)initFileCache.get("REGRESSION", "FeatureDiscoveryAlgorithm");

            String positionalBiasType = (String)initFileCache.get("AFFINITY MODELS", "PositionalBiasType");

            //boolean performStrandDetection = true;
            boolean performStrandDetection = (boolean)initFileCache.getBoolean("AFFINITY MODELS", "PerformStrandDetection");

            String intensitiesDir = (String)initFileCache.get("DATA", "ProbeIntensitiesDir");
            //String probeLocsDir = (String)initFileCache.get("DATA", "ProbeLocsDir");
            String outputDir = (String)initFileCache.get("OUTPUT", "Directory");
            String resultsDir = intensitiesDir + File.separator + outputDir;

            // get the strand from the init file
            WeightMatrixTools.BindingStrand effectiveStrand = reduceData.getStrand();
            WeightMatrixTools.BindingStrand strand = effectiveStrand;
            boolean isRevCompPalindrome = false;
            double featuresScalar = 1.0;
            double eToMu = 0;

            MutableDouble alpha = new MutableDouble(1.0);
            MutableDouble beta = new MutableDouble(1.0);
            MutableDouble nonSpecKa = new MutableDouble(Math.pow(10, -12));

            WeightMatrixTools.BothStrandsCalc calc = reduceData.getBothStrandsCalc();

            if ((proteinLabel == null) || proteinLabel.equals("Count") || proteinLabel.equals("relAffinity")) {
                proteinLabel = FileTools.stripExtension(this.initFile);
            }

            String posBiasAlgorithm;
            if ((nonlinearAlgorithm == null) || nonlinearAlgorithm.equals("")) {
                posBiasAlgorithm = linearAlgorithm;
            }
            else {
                posBiasAlgorithm = nonlinearAlgorithm;
            }

            // Get the motif length from the motif (since the seedSymList may be null)
            String[] correspondSymbols = null;

            char[] kmerModelLengths = initFileCache.getChars("AFFINITY MODELS", "KmerModelLengths");
            String kmerModelPosBiases = (String) initFileCache.get("AFFINITY MODELS", "KmerModelPosBiases");

            double[][][] kmerModelWeights = null;
            boolean[] isLearnedPosBias = null;
            boolean[] isUniformPosBias = null;

            int motifLength = -1;

            // the normalized L1 distance between the raw and smooth positional bias weights
            double posBiasWeightsL1 = 0;

            double rSquared = 0;
            double noPosWeightsRSquared = -1;

            //else if (psam != null) {
            if (psam != null) {
                motifLength = psam.columns();
            }
            else if (seedSymList != null) {
                motifLength = seedSymList.length();
            }
            else if (motif != null) {
                correspondSymbols = whiteSpacePattern.split(motif);
                motifLength = correspondSymbols.length;
            }
            else {
                System.out.println("\nError: I don't have a motif length!!");
                System.exit(-1);
            }

            if (kmerModelLengths != null) {

                kmerModelWeights = new double[kmerModelLengths.length][][];
                isLearnedPosBias = new boolean[kmerModelLengths.length];
                isUniformPosBias = new boolean[kmerModelLengths.length];

                String[] kmerModelPosBiasesArray  = whiteSpacePattern.split(kmerModelPosBiases);
                for (int i=0; i < kmerModelLengths.length; i++) {
                    if (kmerModelPosBiasesArray[i].startsWith("l") || kmerModelPosBiasesArray[i].equalsIgnoreCase("learned")) {
                        isLearnedPosBias[i] = true;
                    }
                    if (kmerModelPosBiasesArray[i].startsWith("u") || kmerModelPosBiasesArray[i].equalsIgnoreCase("uniform")) {
                        isUniformPosBias[i] = true;
                    }
                }
            }

            System.out.println("\nunifiedColumn = "+unifiedColumn);
            int probeSeqLengths = reduceData.getSeqLength(1, unifiedColumn);
            System.out.println("\nprobeSeqLengths = "+probeSeqLengths);

            // cache Count files for each motif (SymbolList or a weight matrix) (creates if they don't exist)
            //reduceData.getCountsMaps(motifLength, 1);
            // Distribution backgroundDist = HMMTools.createUniformDistribution(this.alphabet, false);
            Distribution backgroundDist = new SimpleDistribution(alphabet);
            DistributionTools.allWeightsTo(backgroundDist, 1.0);

            // do not cache counts arrays yet
            reduceData.setCacheCountsArraysFlag(false);

            // init the R engine
            if (rengine == null) {
                initR();
            }

            ////////////////////////////////////////////////////////////////////////////////
            // perform seed search if needed
            ////////////////////////////////////////////////////////////////////////////////

            // if there is no psam and no seedSymList, then we must find the best seed(s)
            if ((psam == null) && (seedSymList == null)) {

                out.print("Seed for motif not given. Now retrieving kmers of length "+motifLength+" found in the probe-associated sequences ...");
                KmerMatrix kmerMatrix = reduceData.getKmerMatrix(motifLength, 1, unifiedColumn);
                out.println(" Done.\n");

                String seedRegressionAlgorithm = (String)initFileCache.get("REGRESSION", "SeedRegressionAlgorithm");

                Table kmerFitTable = new Table();
                Table repeatKmerFitTable = new Table();
                Table palinKmerFitTable = new Table();
                Table rcPalinKmerFitTable = new Table();
                Table nonRcPalinKmerFitTable = new Table();

                String seedRegressionText = null;
                String sortText = null;
                String kmerType = "";
                int numberOfSeedSequences = 0;

                boolean repeat = false;
                boolean palindrome = false;
                boolean rcPalindrome = false;

                if ((seedType != null) && seedType.equalsIgnoreCase("repeat")) {
                    kmerType = "repeat";
                    repeat = true;
                }
                else if ((seedType != null) && seedType.equalsIgnoreCase("palindrome")) {
                    kmerType = "palindromic";
                    palindrome = true;
                }
                else if ((seedType != null) && seedType.equalsIgnoreCase("rcPalindrome")) {
                    kmerType = "reverse-complement palindromic";
                    rcPalindrome = true;
                }

                if (seedRegressionAlgorithm.equalsIgnoreCase("sparseRegression")
                    || seedRegressionAlgorithm.equalsIgnoreCase("robustGradDescent")) {

                    seedRegressionText = "sparse K-mer regression";
                    sortText = "affinity coefficients";

                    char[] allKmers = kmerMatrix.makeAllKmers(motifLength);

                    if (kmerMatrix.wordToProbesMatrix == null) {
                        System.out.print("Creating K-mer to Probes lookup table...");
                        kmerMatrix.wordToProbesMatrix = kmerMatrix.makeKmerToProbesMatrix(
                            motifLength,
                            allKmers,
                            false); // includeRevComps
                        System.out.println("Done.");
                    }

                    double[] intensities = reduceData.getUsedIntensities(1, unifiedColumn);

                    int[] psamStartPositions = getStartPositions(probeSeqLengths, motifLength);
                    int[] kmerStartPositions = noOverHangStartPos(psamStartPositions);

                    // create uniform 1.0 weights
                    double[][] noOverHangWeights = null;
                    noOverHangWeights = new double[2][kmerStartPositions.length];
                    Arrays.fill(noOverHangWeights[0], 1.0);
                    Arrays.fill(noOverHangWeights[1], 1.0);

                    double[] kmerToAffinityMatrix = null;
                    if (seedRegressionAlgorithm.equalsIgnoreCase("sparseRegression")) {

                        // perform sprase K-mer regression
                        kmerToAffinityMatrix = fitSparseKmerModel(
                            motifLength,
                            noOverHangWeights,
                            allKmers,
                            kmerMatrix.wordToProbesMatrix,
                            intensities,
                            kmerMatrix,
                            4095, // maxColumns
                            .20, // fractionDrops,
                            false, // fsam.includeRevComps[z],
                            true); // intercept
                    }
                    else {
                        // perform robust gradient descent
                        kmerToAffinityMatrix = fitGradDescentKmerModel(
                            motifLength,
                            noOverHangWeights,
                            allKmers,
                            kmerMatrix.wordToProbesMatrix,
                            intensities,
                            kmerMatrix,
                            .30, // fractionDrops,
                            false, // fsam.includeRevComps[z],
                            "avg", // initValFlag,
                            100, // numIters,
                            .15, // initTopDrop,
                            .15, // initBottomDrop,
                            .10,  // percentStep,
                            //1.3, // probesPerKmerMultFactor,
                            .5, // probesPerKmerMultFactor,
                            true, // nonNegativeFlag,
                            true); // intercept

                    }

                    // put the K-mer affinites in a table and sort
                    out.println("\nSorting the "+motifLength+"-mers by their estimated seed affinities...");
                    for (int i=0; i < kmerMatrix.wordToProbesMatrix.length; i++) {
                        java.util.List keyValueList = Arrays.asList(
                            new Long(i),
                            new Double(kmerToAffinityMatrix[i]));

                        kmerFitTable.add(keyValueList);
                    }
                    kmerFitTable.sort(1); // ascending
                    kmerFitTable.reverse(); // now descending
                    out.println(" Done.");

                    // Go thru all the K-mers and divvy out palindromes
                    for (int i=0; i < kmerFitTable.rows(); i++) {
                        java.util.List keyValueList = kmerFitTable.getRow(i);
                        long word = ((Long)keyValueList.get(0)).longValue();

                        if (kmerMatrix.isRepeat(word, motifLength)) {
                            repeatKmerFitTable.add(keyValueList);
                        }

                        if (kmerMatrix.isSelfReverse(word, motifLength)) {
                            palinKmerFitTable.add(keyValueList);
                        }

                        // else
                        if (kmerMatrix.isSelfReverseComplement(word, motifLength)) {
                            rcPalinKmerFitTable.add(keyValueList);
                        }
                        else {
                            nonRcPalinKmerFitTable.add(keyValueList);
                        }

                    }

                }
                else { // univariate pearson correlation or regression
                    // get all non-zero kmers from word Counts
                    // perform motifReduce or Pure

                    String numberOfSeedSequencesString = (String)initFileCache.get("REGRESSION", "NumberOfSeedSequences");
                    if (numberOfSeedSequencesString.equalsIgnoreCase("ALL")) {
                        // use just the first experiment for now
                        numberOfSeedSequences = reduceData.numberOfUsedProbes(1, unifiedColumn);
                        out.print("Retrieving all kmers of length "+motifLength+" found in the probe-associated sequences ...");
                    }
                    else {
                        numberOfSeedSequences = Integer.parseInt(numberOfSeedSequencesString);
                        out.print("Retrieving all kmers of length "+motifLength+" found in the top "+numberOfSeedSequences+" intensity-ordered probe-associated sequences ...");
                    }

                    long[] wordIndexes = kmerMatrix.getWordIndexes(numberOfSeedSequences);
                    out.println(" Done.\n");

                    double[] intensities = null;
                    double[] deltaIntensities = null;
                    double stdevIntensities = -1;

                    if (seedRegressionAlgorithm.equalsIgnoreCase("pearson")) {
                        seedRegressionText = "Pearson r correlations";
                        sortText = "Pearson r correlations";
                        // Get the intensities and set them to y
                        intensities = reduceData.getUsedIntensities(1, unifiedColumn);
                        //FileTools.write(intensities, "usedIntensities.txt", false);

                        deltaIntensities = MathTools.getDelta(intensities);
                        stdevIntensities = MathTools.getStdev(deltaIntensities);
                    }
                    else {
                        seedRegressionText = "univariate regressions";
                        sortText = "R-squared";
                    }

                    if ((seedType != null) && seedType.equalsIgnoreCase("palindrome")) {
                        out.println("Will calculate "+seedRegressionText+" for only Kmers that are "+kmerType+".");
                    }
                    else if ((seedType != null) && seedType.equalsIgnoreCase("rcPalindrome")) {
                        out.println("Will calculate "+seedRegressionText+" for only Kmers that are "+kmerType+".");
                    }

                    out.println("Calculating "+seedRegressionText+" on "+wordIndexes.length+" retrieved "+motifLength+"-mers with the model (intensity ~ kmer-count) ...");
                    int wordCounter = 0;
                    for (long wordIndex : wordIndexes) {

                        wordCounter++;
                        if((wordCounter % 500) == 0 ) {
                            out.println("\n\tCalculated "+wordCounter+" out of "+wordIndexes.length+" "+seedRegressionText+" so far...");

                            //out.print("Now sorting the kmers by their R-squared ...");
                            //kmerFitTable.sort(1); //sort by ascending R-squared values
                            //kmerFitTable.reverse(); // get descending order
                            //kmerFitTable.removeRows(10, kmerFitTable.rows()); // remove all but the top 10
                            //out.println(" Done.\n");
                            //break;
                        }

                        SymbolList wordSymList = kmerMatrix.getSymbolList(wordIndex, motifLength);
                        //SymbolList wordSymList = DNATools.createDNA("ACGTCA");
                        //SymbolList wordSymList = DNATools.createDNA("CCCGTG");

                        boolean isRepeat = SymbolListTools.isRepeat(wordSymList);
                        boolean isSelfReverse = SymbolListTools.isSelfReverse(wordSymList);
                        boolean isSelfReverseComplement = SymbolListTools.isSelfReverseComplement(wordSymList, complementTable);

                        if ((repeat) && (!isRepeat)) {
                            continue;
                        }
                        else if ((palindrome) && (!isSelfReverse)) {
                            continue;
                        }
                        else if ((rcPalindrome) && (!isSelfReverseComplement)) {
                            continue;
                        }

                        double result = 0;
                        if (seedRegressionAlgorithm.equalsIgnoreCase("pearson")) {

                            double[] motifCounts =  reduceData.getCountsArray(
                                wordSymList,
                                null, // hammingSphere,
                                null, // startPositions[posIndex],
                                null, // mandatoryColumns,
                                null, //positionalWeights,
                                WeightMatrixTools.BindingStrand.BOTH,
                                calc,
                                1,
                                unifiedColumn);

                            //result = MathTools.pearsonCorrelation(intensities, motifCounts);
                            result = MathTools.pearsonCorrelation(deltaIntensities, stdevIntensities, motifCounts);
                        }
                        else {
                            LinkedHashSet<SymbolList> aMotifHashSet = new LinkedHashSet<SymbolList>(1);
                            aMotifHashSet.add(wordSymList);
                            out.setMinLevel((byte)101);
                            RegressionFitData fitData = fitModel(
                                aMotifHashSet,
                                null, //startPositions
                                null, //positionalWeights
                                "", //motifsFamily
                                //linearAlgorithm,
                                seedRegressionAlgorithm,
                                WeightMatrixTools.BindingStrand.BOTH,
                                WeightMatrixTools.BindingStrand.BOTH,
                                calc,
                                null, //startParams
                                null, //refSymList
                                Double.NaN, //refKaReal
                                unifiedColumn);
                            out.setMinLevel((byte)25);
                            result = fitData.rSquared;
                        }

                        java.util.List keyValueList = Arrays.asList(new Long(wordIndex), new Double(result));
                        kmerFitTable.add(keyValueList);

                        if (isRepeat) {
                            repeatKmerFitTable.add(keyValueList);
                        }
                        if (isSelfReverse) {
                            palinKmerFitTable.add(keyValueList);
                        }

                        // else
                        if (isSelfReverseComplement) {
                            rcPalinKmerFitTable.add(keyValueList);
                        }
                        else {
                            nonRcPalinKmerFitTable.add(keyValueList);
                        }

                        //out.print("\n"+wordSymList.seqString()+"\t"+result);


                    }
                    out.println(" Done.\n");

                    out.print("Now sorting the kmers by their "+sortText+" ...");
                    kmerFitTable.sort(1); //sort by ascending R-squared values
                    kmerFitTable.reverse(); // get descending order

                    repeatKmerFitTable.sort(1); //sort by ascending R-squared values
                    repeatKmerFitTable.reverse(); // get descending order

                    palinKmerFitTable.sort(1); //sort by ascending R-squared values
                    palinKmerFitTable.reverse(); // get descending order

                    rcPalinKmerFitTable.sort(1); //sort by ascending R-squared values
                    rcPalinKmerFitTable.reverse(); // get descending order

                    nonRcPalinKmerFitTable.sort(1); //sort by ascending R-squared values
                    nonRcPalinKmerFitTable.reverse(); // get descending order

                    out.println(" Done.\n");

                }

                StringBuffer topSeedsBuffer = new StringBuffer();
                Table printKmerFitTable = null;
                String topSeedsString = null;

                kmerFitTable.removeRows(100, kmerFitTable.rows()); // remove all but the top 100
                repeatKmerFitTable.removeRows(100, repeatKmerFitTable.rows()); // remove all but the top 100
                palinKmerFitTable.removeRows(100, palinKmerFitTable.rows()); // remove all but the top 100
                rcPalinKmerFitTable.removeRows(100, rcPalinKmerFitTable.rows()); // remove all but the top 100
                nonRcPalinKmerFitTable.removeRows(100, nonRcPalinKmerFitTable.rows()); // remove all but the top 100


                ///////////////////////////////////////////////////////////////////////////
                // Display and save the the different type of seeds in ranked order
                ///////////////////////////////////////////////////////////////////////////
                String[] kmerTypes = {"",    "repeat", "palindromic", "reverse-complement palindromic"};
                String[] fileNames = {"all", "repeat", "palindromic", "rcPalindromic"};

                for (int i=0; i < kmerTypes.length; i++) {

                    topSeedsBuffer = new StringBuffer();

                    if (seedRegressionAlgorithm.equalsIgnoreCase("sparseRegression")) {
                        topSeedsBuffer.append("The top 100 "+kmerTypes[i]+" k-mer seeds with the highest multivariate affinity-coefficients for protein "+proteinLabel+" are:\n");
                    }
                    else if (seedRegressionAlgorithm.equalsIgnoreCase("robustGradDescent")) {
                        topSeedsBuffer.append("The top 100 "+kmerTypes[i]+" k-mer seeds with the highest multivariate affinity-coefficients for protein "+proteinLabel+" are:\n");
                    }
                    else {
                        topSeedsBuffer.append("The top 100 "+kmerTypes[i]+" kmer seeds (of those found in the top "+numberOfSeedSequences+" intensity-ordered probe-associated sequences) that best fit the intensity data (using "+sortText+") for "+proteinLabel+" are:\n");
                    }

                    if (fileNames[i].equals("all")) {
                        printKmerFitTable = kmerFitTable;
                    }
                    else if (fileNames[i].equals("repeat")) {
                        printKmerFitTable = repeatKmerFitTable;
                    }
                    else if (fileNames[i].equals("palindromic")) {
                        printKmerFitTable = palinKmerFitTable;
                    }
                    else if (fileNames[i].equals("rcPalindromic")) {
                        printKmerFitTable = rcPalinKmerFitTable;
                    }

                    //Iterator iter = printKmerPrintTable.iterator();
                    int j = 1;
                    for (Iterator iter = printKmerFitTable.iterator(); iter.hasNext(); ) {
                        //for (java.util.List<Object> rowList : printKmerPrintTable) {
                        java.util.List rowList = (java.util.List)iter.next();
                        String seedString = kmerMatrix.getWordString(((Long)rowList.get(0)).longValue(), motifLength);
                        topSeedsBuffer.append("\t"+j+"\t"+seedString+"\t"+((Double)rowList.get(1)).doubleValue()+"\n");
                        j++;
                    }

                    topSeedsString = topSeedsBuffer.toString();
                    out.println("\n"+topSeedsString);
                    FileTools.write(topSeedsString, resultsDir + File.separator + proteinLabel+"."+fileNames[i]+"TopSeeds.txt", false);

                }

                ///////////////////////////////////////////////////////////////////////////
                // get the top seed
                ///////////////////////////////////////////////////////////////////////////
                int rankedKmerForSeed = initFileCache.getInt("REGRESSION", "RankedKmerForSeed"); // 1-based
                out.println("\nWill be using top ranked rc-palindromic or non-rc-palindromic K-mer number "+rankedKmerForSeed+" as the seed.");

                rankedKmerForSeed--; // 0-based

                java.util.List firstRowList;
                //firstRowList = (java.util.List)kmerFitTable.iterator().next();
                firstRowList = (java.util.List)kmerFitTable.getRow(rankedKmerForSeed);
                SymbolList anyKmerTopSymList = kmerMatrix.getSymbolList(((Long)firstRowList.get(0)).longValue(), motifLength);
                double anyKmerTopCorrelation = ((Double)firstRowList.get(1)).doubleValue();

                // SymbolList palinKmerTopSymList = null;
                // double palinKmerTopCorrelation = 0;
                // if ((palinKmerFitTable != null) && (!palinKmerFitTable.isEmpty())) {
                //     //firstRowList = (java.util.List)palinKmerFitTable.iterator().next();
                //     firstRowList = (java.util.List)palinKmerFitTable.getRow(rankedKmerForSeed);
                //     palinKmerTopSymList = kmerMatrix.getSymbolList(((Long)firstRowList.get(0)).longValue(), motifLength);
                //     palinKmerTopCorrelation = ((Double)firstRowList.get(1)).doubleValue();
                // }

                SymbolList rcPalinKmerTopSymList = null;
                double rcPalinKmerTopCorrelation = 0;
                if ((rcPalinKmerFitTable != null) && (!rcPalinKmerFitTable.isEmpty())) {
                    //firstRowList = (java.util.List)rcPalinKmerFitTable.iterator().next();
                    firstRowList = (java.util.List)rcPalinKmerFitTable.getRow(rankedKmerForSeed);
                    rcPalinKmerTopSymList = kmerMatrix.getSymbolList(((Long)firstRowList.get(0)).longValue(), motifLength);
                    rcPalinKmerTopCorrelation = ((Double)firstRowList.get(1)).doubleValue();
                }

                SymbolList nonRcPalinKmerTopSymList = null;
                double nonRcPalinKmerTopCorrelation = 0;
                if ((nonRcPalinKmerFitTable != null) && (!nonRcPalinKmerFitTable.isEmpty())) {
                    //firstRowList = (java.util.List)nonRcPalinKmerFitTable.iterator().next();
                    firstRowList = (java.util.List)nonRcPalinKmerFitTable.getRow(rankedKmerForSeed);
                    nonRcPalinKmerTopSymList = kmerMatrix.getSymbolList(((Long)firstRowList.get(0)).longValue(), motifLength);
                    nonRcPalinKmerTopCorrelation = ((Double)firstRowList.get(1)).doubleValue();
                }

                double palindromicSeedThresh = initFileCache.getDouble("REGRESSION", "PalindromicSeedThresh");

                double rcPalinRelAffinity = rcPalinKmerTopCorrelation / nonRcPalinKmerTopCorrelation;
                if (rcPalinKmerTopCorrelation >= (nonRcPalinKmerTopCorrelation * palindromicSeedThresh)) {
                    out.println("\nTop ranked rc-palindromic seed "+rcPalinKmerTopSymList.seqString()+" has relative affinity of "+rcPalinRelAffinity+", which is over the threshold of "+palindromicSeedThresh+".");
                    seedSymList = rcPalinKmerTopSymList;
                }
                // else if (palinKmerTopCorrelation >= (nonRcPalinKmerTopCorrelation * palindromicSeedThresh)) {
                //     seedSymList = palinKmerTopSymList;
                // }
                else {
                    out.println("\nTop ranked rc-palindromic seed "+rcPalinKmerTopSymList.seqString()+" has relative affinity of "+rcPalinRelAffinity+", which is NOT over the threshold of "+palindromicSeedThresh+".");
                    seedSymList = nonRcPalinKmerTopSymList;
                }

                out.println("\nWill be using top-ranked K-mer "+seedSymList.seqString()+" as the seed.");

                ///////////////////////////////////////////////////////////////////////
                //seedSymList = DNATools.reverseComplement(seedSymList);
                ///////////////////////////////////////////////////////////////////////

            }


            ////////////////////////////////////////////////////////////////////////////////
            // Get seeds only!
            ////////////////////////////////////////////////////////////////////////////////
//             if (true) {
//                 return(null);
//             }

            // We are done with the BIG table so free the memory
            //kmerFitTable = null;

            ////////////////////////////////////////////////////////////////////////////////
            // Initialize for HammingDistance spheres (sets) and HammingPSAMs
            ////////////////////////////////////////////////////////////////////////////////
            LinkedHashSet<SymbolList> hammingDist1WithRevCompsSet = null;
            LinkedHashSet<SymbolList> hammingDist1WithoutRevCompsSet = null;

            // 0 = seed alone
            // 1 = Hamming-1 sphere
            // 2 = hamming-2 sphere .....
            ArrayList<LinkedHashSet<SymbolList>> hammingSpheres = new ArrayList<LinkedHashSet<SymbolList>>(4);

            WeightMatrix hamming1PSAM = null;
            WeightMatrix hamming1PosStrandPSAM = null;
            WeightMatrix hamming1NegStrandPSAM = null;
            WeightMatrix hamming1BothStrandsPSAM = null;

            WeightMatrix hamming2PSAM = null;
            WeightMatrix hamming2PosStrandPSAM = null;
            WeightMatrix hamming2NegStrandPSAM = null;
            WeightMatrix hamming2BothStrandsPSAM = null;

            double posStrandRSquared = 0;
            double negStrandRSquared = 0;
            double bothStrandsRSquared = 0;

            RegressionFitData regressionFitDataPosStrand = null;
            RegressionFitData regressionFitDataNegStrand = null;
            RegressionFitData regressionFitDataBothStrands = null;

            GridBagLayoutFrame logoFrame;

            ////////////////////////////////////////////////////////////////////////////////
            // Initialize Logo Frame
            ////////////////////////////////////////////////////////////////////////////////
            GridBagLayoutFrame strandLogoFrame = null;
            if (displayMotifs) {
                strandLogoFrame = createLogoFrame();
            }

            double[][] positionalWeights = initPositionalBias;

            // kmerModelWeights = positionalWeights;
            if (kmerModelLengths != null) {
                for(int z=0; z < kmerModelLengths.length; z++) {
                    if ((isLearnedPosBias[z]) && (motifLength == kmerModelLengths[z])) {
                        kmerModelWeights[z] = positionalWeights;
                    }
                }
            }

//             double[][] positionalWeights = new double[2][16 + 2*(motifLength-1)];
//             MathTools.fill(positionalWeights, 1.0);

            if (psam == null) {
                // perform Hamming Reduce to get the PSAM

                ////////////////////////////////////////////////////////////////////////////////
                // perform Hamming Reduce on positive strand if desired
                ////////////////////////////////////////////////////////////////////////////////

                if (effectiveStrand == WeightMatrixTools.BindingStrand.POS || performStrandDetection) {

                    if (cacheCounts && effectiveStrand == WeightMatrixTools.BindingStrand.POS) { // start caching counts arrays
                        reduceData.setCacheCountsArraysFlag(true);
                    }

                    // get all the sequences within some hamming distance from the seed sequence
                    // including revComps
                    if (hammingDist1WithRevCompsSet == null) {
                        hammingDist1WithRevCompsSet = SymbolListTools.getHammingSphere(seedSymList, 1, true, complementTable);
                    }
                    //out.println("hammingDistSeqMap = "+symListCollectionDump(hammingDistanceSeqMap));

                    // send the intensities and counts to R to perform regression
                    regressionFitDataPosStrand = fitModel(
                        hammingDist1WithRevCompsSet,
                        null, //startPositions
                        positionalWeights,
                        "single point-mutations (Hamming-1 sphere) about the seed sequence on the positive strand",
                        linearAlgorithm,
                        WeightMatrixTools.BindingStrand.POS,
                        WeightMatrixTools.BindingStrand.POS,
                        calc,
                        null, //startParams
                        null, //refSymList
                        Double.NaN, //refKaReal
                        unifiedColumn);

                    // coefficients are as returned in R
                    // First coefficient is the intercept, then the rest are in independent var order
                    out.print(regressionFitDataPosStrand.toString());

                    // Create the weight matrix
                    hamming1PosStrandPSAM = getPsam(
                        linearAlgorithm,
                        regressionFitDataPosStrand,
                        minAffinityWeight,
                        backgroundDist,
                        WeightMatrixTools.BindingStrand.POS,
                        WeightMatrixTools.BindingStrand.POS,
                        calc,
                        eToMu,
                        refSymList,
                        referenceKa,
                        unifiedColumn,
                        alpha,
                        beta,
                        nonSpecKa);

                }

                ////////////////////////////////////////////////////////////////////////////////
                // perform Hamming Reduce on negative strand if desired
                ////////////////////////////////////////////////////////////////////////////////

                if (effectiveStrand == WeightMatrixTools.BindingStrand.NEG || performStrandDetection) {

                    if (cacheCounts && effectiveStrand == WeightMatrixTools.BindingStrand.NEG) { // start caching counts arrays
                        reduceData.setCacheCountsArraysFlag(true);
                    }

                    // get all the sequences within some hamming distance from the seed sequence
                    // including revComps
                    if (hammingDist1WithRevCompsSet == null) {
                        hammingDist1WithRevCompsSet = SymbolListTools.getHammingSphere(seedSymList, 1, true, complementTable);
                    }
                    //out.println("hammingDistSeqMap = "+symListCollectionDump(hammingDistanceSeqMap));

                    // send the intensities and counts to R to perform regression
                    regressionFitDataNegStrand = fitModel(
                        hammingDist1WithRevCompsSet,
                        null, //startPositions
                        positionalWeights,
                        "single point-mutations (Hamming-1 sphere) about the seed sequence on the negative strand",
                        linearAlgorithm,
                        WeightMatrixTools.BindingStrand.NEG,
                        WeightMatrixTools.BindingStrand.NEG,
                        calc,
                        null, //startPositions
                        null, //refSymList
                        Double.NaN, //refKaReal
                        unifiedColumn);

                    // coefficients are as returned in R
                    // First coefficient is the intercept, then the rest are in independent var order
                    out.print(regressionFitDataNegStrand.toString());

                    // Create the weight matrix
                    hamming1NegStrandPSAM = getPsam(
                        linearAlgorithm,
                        regressionFitDataNegStrand,
                        minAffinityWeight,
                        backgroundDist,
                        WeightMatrixTools.BindingStrand.NEG,
                        WeightMatrixTools.BindingStrand.NEG,
                        calc,
                        eToMu,
                        refSymList,
                        referenceKa,
                        unifiedColumn,
                        alpha,
                        beta,
                        nonSpecKa);
                }


                ////////////////////////////////////////////////////////////////////////////////
                // perform Hamming Reduce on both strands if desired
                ////////////////////////////////////////////////////////////////////////////////

                if (effectiveStrand == WeightMatrixTools.BindingStrand.BOTH || performStrandDetection) {

                    // get all the sequences within some hamming distance from the seed sequence
                    // including revComps
                    if (hammingDist1WithRevCompsSet == null) {
                        hammingDist1WithRevCompsSet = SymbolListTools.getHammingSphere(seedSymList, 1, true, complementTable);
                    }
                    //out.println("hammingDistSeqMap = "+symListCollectionDump(hammingDistanceSeqMap));

                    // get all the sequences within some hamming distance from the seed sequence
                    // NOT including revComps
//                     if (hammingDist1WithoutRevCompsSet == null) {
//                         hammingDist1WithoutRevCompsSet = SymbolListTools.getHammingSphere(seedSymList, 1, false, complementTable);
//                     }
                    //out.println("hammingDistSeqMap = "+symListCollectionDump(hammingDistanceSeqMap));

                    // send the intensities and counts to R to perform regression
                    regressionFitDataBothStrands = fitModel(
                        hammingDist1WithRevCompsSet,
                        //hammingDist1WithoutRevCompsSet,
                        null, //startPositions
                        positionalWeights,
                        "single point-mutations (Hamming-1 sphere) about the seed sequence on both strands",
                        linearAlgorithm,
                        WeightMatrixTools.BindingStrand.BOTH,
                        WeightMatrixTools.BindingStrand.BOTH,
                        calc,
                        null, //startParams
                        null, //refSymList
                        Double.NaN, //refKaReal
                        unifiedColumn);

                    // coefficients are as returned in R
                    // First coefficient is the intercept, then the rest are in independent var order
                    out.print(regressionFitDataBothStrands.toString());

                    // Create the weight matrix
                    hamming1BothStrandsPSAM = getPsam(
                        linearAlgorithm,
                        regressionFitDataBothStrands,
                        minAffinityWeight,
                        backgroundDist,
                        WeightMatrixTools.BindingStrand.BOTH,
                        WeightMatrixTools.BindingStrand.BOTH,
                        calc,
                        eToMu,
                        refSymList,
                        referenceKa,
                        unifiedColumn,
                        alpha,
                        beta,
                        nonSpecKa);
                }

            } // if (PSAM == null) - done performing hammingReduce

            // use the PSAM that was given
            else {
//                 hamming1PosStrandPSAM = getPsam(psam, null, minAffinityWeight, backgroundDist, WeightMatrixTools.BindingStrand.POS, calc, eToMu);
//                 hamming1NegStrandPSAM = getPsam(psam, null, minAffinityWeight, backgroundDist, WeightMatrixTools.BindingStrand.NEG, calc, eToMu);
//                 hamming1BothStrandsPSAM = getPsam(psam, null, minAffinityWeight, backgroundDist, WeightMatrixTools.BindingStrand.BOTH, calc, eToMu);
                hamming1PosStrandPSAM = psam;
                hamming1NegStrandPSAM = psam;
                hamming1BothStrandsPSAM = psam;
                seedSymList = WeightMatrixTools.getMostLikelySymList(psam);

            }

            if (displayMotifs) {
                if (hamming1PosStrandPSAM != null) {
                    logoFrame = display(proteinLabel+": "+"Initial Hamming-1 positive-strand PSAM", hamming1PosStrandPSAM);
                }
                if (hamming1NegStrandPSAM != null) {
                    logoFrame = display(proteinLabel+": "+"Initial Hamming-1 negative-strand PSAM", hamming1NegStrandPSAM);
                }
                if (hamming1BothStrandsPSAM != null) {
                    logoFrame = display(proteinLabel+": "+"Initial Hamming-1 both-strands PSAM", hamming1BothStrandsPSAM);
                }
            }

            ////////////////////////////////////////////////////////////////////////////////
            // Now perform fits on the PSAMs
            ////////////////////////////////////////////////////////////////////////////////
            RegressionFitData regressionFitData = null;
            LinkedHashMap<Object, LinkedHashSet<SymbolList>> motifToHammingSphereMap = null;
            int[] startPositions = getStartPositions(probeSeqLengths, motifLength);

            if (effectiveStrand == WeightMatrixTools.BindingStrand.POS || performStrandDetection) {

                double revCompSimilarity = WeightMatrixTools.revCompSimilarity(hamming1PosStrandPSAM);

                // fit with this PSAM
                regressionFitData = getRegressionFit(
                    hamming1PosStrandPSAM,
                    positionalWeights,
                    "fit with the Positive-Strand PSAM",
                    "lm",
                    WeightMatrixTools.BindingStrand.POS,
                    WeightMatrixTools.BindingStrand.POS,
                    calc,
                    unifiedColumn);

                rSquared = regressionFitData.rSquared;

                nonSpecKa.setValue(getRelativeKaNonSpec(
                    regressionFitData.getIntercept(),
                    regressionFitData.getCoefficient(hamming1PosStrandPSAM),
                    positionalWeights,
                    eToMu,
                    probeSeqLengths,
                    motifLength,
                    WeightMatrixTools.BindingStrand.POS,
                    calc,
                    revCompSimilarity));

                // get the positional bias
                double[][] rawWeights = getPositionalBias(
                    //"nnls",
                    //"nls.levMar",
                    posBiasAlgorithm,
                    hamming1PosStrandPSAM,
                    calc,
                    startPositions,
                    WeightMatrixTools.BindingStrand.POS,
                    WeightMatrixTools.BindingStrand.POS,
                    refSymList,
                    referenceKa,
                    unifiedColumn,
                    alpha,
                    beta,
                    nonSpecKa);

                double[][] cleanedWeights = cleanWeights(startPositions, rawWeights);

                double[][] smoothedWeights = smoothWeights(startPositions, cleanedWeights);
                //double[][] smoothedWeights = smoothWeights(startPositions, rawWeights);

                if (kmerModelLengths != null) {
                    for(int z=0; z < kmerModelLengths.length; z++) {
                        if ((isLearnedPosBias[z]) && (motifLength == kmerModelLengths[z])) {
                            kmerModelWeights[z] = smoothedWeights;
                        }
                    }
                }

                // Get the fit for the new hamming1PSAM with the positional Weights
                regressionFitData = getRegressionFit(
                    hamming1PosStrandPSAM,
                    //cleanedWeights,
                    smoothedWeights,
                    "fit with the Positive-Strand PSAM with weights",
                    "lm",
                    WeightMatrixTools.BindingStrand.POS,
                    WeightMatrixTools.BindingStrand.POS,
                    calc,
                    unifiedColumn);

                double posWeightsRSquared = regressionFitData.rSquared;

                nonSpecKa.setValue(getRelativeKaNonSpec(
                    regressionFitData.getIntercept(),
                    regressionFitData.getCoefficient(hamming1PosStrandPSAM),
                    //cleanedWeights,
                    smoothedWeights,
                    eToMu,
                    probeSeqLengths,
                    motifLength,
                    WeightMatrixTools.BindingStrand.POS,
                    calc,
                    revCompSimilarity));

                posStrandRSquared = posWeightsRSquared;

                if (displayMotifs) {
                    String[] datasetLabels = {"Raw Positive Strand", "Cleaned Positive Strand", "Smoothed Positive Strand"};
                    double[][] values = {ArrayTools.toDoubleArray(startPositions), rawWeights[0], cleanedWeights[0], smoothedWeights[0]};
                    BioJavaChart chart = new BioJavaChart("lines", null, "Start Position", "Normalized Coefficients", null, datasetLabels, true, values, -1);

                    addMotif(strandLogoFrame, hamming1PosStrandPSAM, "Hamming PSAM", seedSymList.seqString(), "POS", rSquared, chart, posWeightsRSquared);
                }

            }
            if (effectiveStrand == WeightMatrixTools.BindingStrand.NEG || performStrandDetection) {

                double revCompSimilarity = WeightMatrixTools.revCompSimilarity(hamming1NegStrandPSAM);

                // fit with this PSAM
                regressionFitData = getRegressionFit(
                    hamming1NegStrandPSAM,
                    positionalWeights,
                    "fit with the Negative-Strand PSAM",
                    "lm",
                    WeightMatrixTools.BindingStrand.NEG,
                    WeightMatrixTools.BindingStrand.NEG,
                    calc,
                    unifiedColumn);

                rSquared = regressionFitData.rSquared;

                nonSpecKa.setValue(getRelativeKaNonSpec(
                    regressionFitData.getIntercept(),
                    regressionFitData.getCoefficient(hamming1NegStrandPSAM),
                    positionalWeights,
                    eToMu,
                    probeSeqLengths,
                    motifLength,
                    WeightMatrixTools.BindingStrand.NEG,
                    calc,
                    revCompSimilarity));


                // get the positional bias
                double[][] rawWeights = getPositionalBias(
                    //"nnls",
                    posBiasAlgorithm,
                    hamming1NegStrandPSAM,
                    calc,
                    startPositions,
                    WeightMatrixTools.BindingStrand.NEG,
                    WeightMatrixTools.BindingStrand.NEG,
                    refSymList,
                    referenceKa,
                    unifiedColumn,
                    alpha,
                    beta,
                    nonSpecKa);

                double[][] cleanedWeights = cleanWeights(startPositions, rawWeights);

                double[][] smoothedWeights = smoothWeights(startPositions, cleanedWeights);
                //double[][] smoothedWeights = smoothWeights(startPositions, rawWeights);

                if (kmerModelLengths != null) {
                    for(int z=0; z < kmerModelLengths.length; z++) {
                        if ((isLearnedPosBias[z]) && (motifLength == kmerModelLengths[z])) {
                            kmerModelWeights[z] = smoothedWeights;
                        }
                    }
                }

                // Get the fit for the new hamming1PSAM with the positional Weights
                regressionFitData = getRegressionFit(
                    hamming1NegStrandPSAM,
                    //cleanedWeights,
                    smoothedWeights,
                    "fit with the Negative-Strand PSAM with weights",
                    "lm",
                    WeightMatrixTools.BindingStrand.NEG,
                    WeightMatrixTools.BindingStrand.NEG,
                    calc,
                    unifiedColumn);

                double posWeightsRSquared = regressionFitData.rSquared;

                nonSpecKa.setValue(getRelativeKaNonSpec(
                    regressionFitData.getIntercept(),
                    regressionFitData.getCoefficient(hamming1NegStrandPSAM),
                    //cleanedWeights,
                    smoothedWeights,
                    eToMu,
                    probeSeqLengths,
                    motifLength,
                    WeightMatrixTools.BindingStrand.NEG,
                    calc,
                    revCompSimilarity));

                negStrandRSquared = posWeightsRSquared;

                if (displayMotifs) {
                    String[] datasetLabels = {"Raw Negative Strand", "Cleaned Negative Strand", "Smoothed Negative Strand"};
                    double[][] values = {ArrayTools.toDoubleArray(startPositions), rawWeights[1], cleanedWeights[1], smoothedWeights[1]};
                    BioJavaChart chart = new BioJavaChart("lines", null, "Start Position", "Normalized Coefficients", null, datasetLabels, true, values, -1);

                    addMotif(strandLogoFrame, hamming1NegStrandPSAM, "Hamming PSAM", seedSymList.seqString(), "NEG", rSquared, chart, posWeightsRSquared);
                }

            }
            if (effectiveStrand == WeightMatrixTools.BindingStrand.BOTH || performStrandDetection) {
                isRevCompPalindrome = WeightMatrixTools.isSelfReverseComplement(hamming1BothStrandsPSAM, complementTable);

                WeightMatrixTools.BindingStrand tempStrand;

//                 if (isRevCompPalindrome) {
//                     out.println("\n\nSince the motif is self reverse-complement, we can look only at the positive strand.");
//                     tempStrand = WeightMatrixTools.BindingStrand.POS;
//                 }
//                 else {
//                     tempStrand = WeightMatrixTools.BindingStrand.BOTH;
//                 }
                tempStrand = WeightMatrixTools.BindingStrand.BOTH;


                double revCompSimilarity = WeightMatrixTools.revCompSimilarity(hamming1BothStrandsPSAM);

                // fit with this PSAM
                regressionFitData = getRegressionFit(
                    hamming1BothStrandsPSAM,
                    positionalWeights,
                    "fit with Both-Strands PSAM",
                    "lm",
                    WeightMatrixTools.BindingStrand.BOTH,
                    //WeightMatrixTools.BindingStrand.BOTH,
                    tempStrand,
                    calc,
                    unifiedColumn);

                rSquared = regressionFitData.rSquared;

                nonSpecKa.setValue(getRelativeKaNonSpec(
                    regressionFitData.getIntercept(),
                    regressionFitData.getCoefficient(hamming1BothStrandsPSAM),
                    positionalWeights,
                    eToMu,
                    probeSeqLengths,
                    motifLength,
                    WeightMatrixTools.BindingStrand.BOTH,
                    calc,
                    revCompSimilarity));

                // get the positional bias
                //double[][] rawWeights = getPositionalBias(hamming1BothStrandsPSAM, calc, startPositions, WeightMatrixTools.BindingStrand.BOTH, tempStrand);

//                 double[][] rawWeights = null;
//                 double[][] cleanedWeights = null;
//                 double[][] smoothedWeights = null;
                double[][] rawWeights = positionalWeights;
                double[][] cleanedWeights = positionalWeights;
                double[][] smoothedWeights = positionalWeights;
                double posWeightsRSquared = 0;
                if (((String)initFileCache.get("AFFINITY MODELS", "GetPositionalBias")).equalsIgnoreCase("Yes")) {

                    rawWeights = getPositionalBias(
                        //"nnls",
                        //"nls.levMar",
                        posBiasAlgorithm,
                        hamming1BothStrandsPSAM,
                        calc,
                        startPositions,
                        WeightMatrixTools.BindingStrand.BOTH,
                        //WeightMatrixTools.BindingStrand.BOTH,
                        tempStrand,
                        refSymList,
                        referenceKa,
                        unifiedColumn,
                        alpha,
                        beta,
                        nonSpecKa);

                    cleanedWeights = cleanWeights(startPositions, rawWeights);

                    smoothedWeights = smoothWeights(startPositions, cleanedWeights);
                    //smoothedWeights = smoothWeights(startPositions, rawWeights);

                    posBiasWeightsL1 = getPosBiasWeightsL1(smoothedWeights, cleanedWeights);

                    if (kmerModelLengths != null) {
                        for(int z=0; z < kmerModelLengths.length; z++) {
                            if ((isLearnedPosBias[z]) && (motifLength == kmerModelLengths[z])) {
                                kmerModelWeights[z] = smoothedWeights;
                            }
                        }
                    }

                    // Get the fit for the new hamming1PSAM with the positional Weights
                    regressionFitData = getRegressionFit(
                        hamming1BothStrandsPSAM,
                        //cleanedWeights,
                        smoothedWeights,
                        "fit with the Both-Strands PSAM with weights",
                        "lm",
                        WeightMatrixTools.BindingStrand.BOTH,
                        //WeightMatrixTools.BindingStrand.BOTH,
                        tempStrand,
                        calc,
                        unifiedColumn);

                    posWeightsRSquared = regressionFitData.rSquared;

                    nonSpecKa.setValue(getRelativeKaNonSpec(
                            regressionFitData.getIntercept(),
                            regressionFitData.getCoefficient(hamming1BothStrandsPSAM),
                            //cleanedWeights,
                            smoothedWeights,
                            eToMu,
                            probeSeqLengths,
                            motifLength,
                            WeightMatrixTools.BindingStrand.BOTH,
                            calc,
                            revCompSimilarity));

                    bothStrandsRSquared = posWeightsRSquared;
                }

                if (displayMotifs) {
                    BioJavaChart chart = null;
                    //if (((String)initFileCache.get("AFFINITY MODELS", "GetPositionalBias")).equalsIgnoreCase("Yes")) {
                    if ((positionalWeights != null) || ((String)initFileCache.get("AFFINITY MODELS", "GetPositionalBias")).equalsIgnoreCase("Yes")) {
                        String[] datasetLabels = {"Raw Positive Strand", "Raw Negative Strand", "Cleaned Positive Strand", "Cleaned Negative Strand", "Smoothed Positive Strand", "Smoothed Negative Strand"};
                        double[][] values = {ArrayTools.toDoubleArray(startPositions), rawWeights[0], rawWeights[1], cleanedWeights[0], cleanedWeights[1], smoothedWeights[0], smoothedWeights[1]};
                        chart = new BioJavaChart("lines", null, "Start Position", "Normalized Coefficients", null, datasetLabels, true, values, -1);
                    }
                    String seedString = null;
                    if (seedSymList != null) {
                        seedString = seedSymList.seqString();
                    }
                    addMotif(strandLogoFrame, hamming1BothStrandsPSAM, "Hamming PSAM", seedString, "BOTH", rSquared, chart, posWeightsRSquared);
                }

            }

            ////////////////////////////////////////////////////////////////////////////////
            // If effectiveStrand==Both:
            // Due to strand-binding biases in PBM data.  We cannot use R^2 of the regressions for strand determination
            // Instead, we should use the JSEntropy between the positive and negative strand motifs.
            // FIX THIS!
            ////////////////////////////////////////////////////////////////////////////////

            out.println("\n*******************************************************************************");
            out.println("\nThe R-squared using the positive strand is "+posStrandRSquared+".");
            out.println("\nThe R-squared using the negative strand is "+negStrandRSquared+".");
            out.println("\nThe R-squared using both strands is "+bothStrandsRSquared+".");

            if ((effectiveStrand == WeightMatrixTools.BindingStrand.POS) ||
                (performStrandDetection && (posStrandRSquared > 1000*negStrandRSquared))) {

                out.println("\n\nThe R-squareds show that this is a positive-strand motif.");

                effectiveStrand = WeightMatrixTools.BindingStrand.POS;
                strand = WeightMatrixTools.BindingStrand.POS;
                includeRevComps = true;
                hamming1PSAM = hamming1PosStrandPSAM;
                //regressionFitData = regressionFitDataPosStrand;
                //hamming2PSAM = hamming2PosStrandPSAM;

            }

            else if ((effectiveStrand == WeightMatrixTools.BindingStrand.NEG) ||
                (performStrandDetection && (negStrandRSquared > 1000*posStrandRSquared))) {
                out.println("\n\nThe R-squareds show that this is a negative-strand motif.");

                effectiveStrand = WeightMatrixTools.BindingStrand.NEG;
                strand = WeightMatrixTools.BindingStrand.NEG;
                includeRevComps = true;
                hamming1PSAM = hamming1NegStrandPSAM;
                //regressionFitData = regressionFitDataNegStrand;
                //hamming2PSAM = hamming2NegStrandPSAM;

            }

            else { // effectiveStrand == WeightMatrixTools.BindingStrand.BOTH;

                //                     if ((effectiveStrand == WeightMatrixTools.BindingStrand.BOTH) &&
                //                     (regressionFitDataBothStrands.rSquared >= (1.35 * regressionFitDataPosStrand.rSquared)) &&
                //                     (regressionFitDataBothStrands.rSquared >= (1.35 * regressionFitDataNegStrand.rSquared))) {

                out.println("\n\nThe R-squareds show that this is a dual-strand motif (present on both strands).");

                // If the PSAM is a reverse-complement palindrome, then can look at just the positive
                // or negative strand
//                 if (isRevCompPalindrome) {
//                     out.println("Since the motif is self reverse-complement, we can look just at just the positive strand.");
//                     strand = WeightMatrixTools.BindingStrand.POS;
//                     //strand = WeightMatrixTools.BindingStrand.NEG;
//                     //featuresScalar = 0.5; // If you build with only pos, but plan to score with both
//                 }
//                 else {
//                     strand = WeightMatrixTools.BindingStrand.BOTH;
//                 }

                String seedString = seedSymList.seqString();

                if (seedString.contains("ccccc") || seedString.contains("CCCCC")) {
                    out.println("\n\nWARNING! - The seed contains CCCCC, which is heavily depleted on the NEGATIVE strand!!  Must discard the NEGATIVE strand, and use only the POSITIVE strand.");

                    effectiveStrand = WeightMatrixTools.BindingStrand.POS;
                    strand = WeightMatrixTools.BindingStrand.POS;
                    includeRevComps = true;
                    hamming1PSAM = hamming1PosStrandPSAM;
                }
                else if (seedString.contains("ggggg") || seedString.contains("GGGGG")) {
                    out.println("\n\nWARNING! - The seed contains GGGGG, which is heavily depleted on the POSITIVE strand!!  Must discard the POSITIVE strand, and use only the NEGATIVE strand.");

                    effectiveStrand = WeightMatrixTools.BindingStrand.NEG;
                    strand = WeightMatrixTools.BindingStrand.NEG;
                    includeRevComps = true;
                    hamming1PSAM = hamming1NegStrandPSAM;
                }
                else {
                    strand = WeightMatrixTools.BindingStrand.BOTH;

                    includeRevComps = false;
                    hamming1PSAM = hamming1BothStrandsPSAM;
                    //regressionFitData = regressionFitDataBothStrands;
                    //hamming2PSAM = hamming2BothStrandsPSAM;
                }

            }

            out.println("\n*******************************************************************************");

            // save the weight matrix as PNG and XML
            //GraphicsTools.writeImage(logoPanel, resultsDir + File.separator + "hammingPsam.png");
            WeightMatrixTools.writeToXML(hamming1PSAM, resultsDir + File.separator + "hammingPsamBeforeExt.xml");

            //                 Object revComp1 = (Object)WeightMatrixTools.reverseComplement(hamming1PSAM, complementTable());
            //                 WeightMatrix hamming1PSAMmod = WeightMatrixTools.insert((WeightMatrix)hamming1PSAM, 1, DNATools.a());
            //                 out.println((byte)100, "\nhamming1PSAMmod = "+StringTools.toString(hamming1PSAMmod));

            //                 LinkedHashSet<Object> LHS1 = new LinkedHashSet<Object>();
            //                 LHS1.add((Object)hamming1PSAMmod);

            //                 WeightMatrix revComp1Mod = WeightMatrixTools.insert((WeightMatrix)revComp1, 5, DNATools.t());
            //                 Object revComp2 = (Object)WeightMatrixTools.reverseComplement((WeightMatrix)revComp1Mod, complementTable());
            //                 out.println((byte)100, "\nrevComp2 = "+StringTools.toString((WeightMatrix)revComp2));

            //                 if (LHS1.contains((Object)revComp2)) {
            //                     out.println("\nFound the reverseComplement!");
            //                 }
            //                 else {
            //                     out.println("\nDid not find the reverseComplement!");
            //                 }

            //                 if (((Object)hamming1PSAMmod).equals((Object)revComp2)) {
            //                     out.println("\nThey are equal!");
            //                 }
            //                 else {
            //                     out.println("\nThey are not equal!");
            //                 }


//             GridBagLayoutFrame psamsFrame1 = null;
//             if (displayMotifs) {
//                 psamsFrame1 = new GridBagLayoutFrame(proteinLabel+": "+"ColumnREDUCE PSAMs before Extension", true);
//             }

            //strand = WeightMatrixTools.BindingStrand.POS;

            LinkedHashMap<Object, boolean[][]> motifsToMandatoriesMap;

            // ////////////////////////////////////////////////////////////////////////////////
            // // Perform regression to compare Motifs (currently must be psams)
            // ////////////////////////////////////////////////////////////////////////////////
            // if (psamsArray != null) {

            //     String spacerLenRegressionAlgorithm = (String)initFileCache.get("REGRESSION", "spacerLenRegressionAlgorithm");

            //     motifsToMandatoriesMap = new LinkedHashMap<Object, boolean[][]>();
            //     motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();

            //     //                     if (useHammingDistForApprox) {
            //     //                         motifToHammingSphereMap.put(hamming1PSAM, hammingSpheres.get(hammingDistForApprox));
            //     //                     }
            //     //                     else {
            //     //                         motifToHammingSphereMap.put(hamming1PSAM, null);
            //     //                     }

            //     //Distribution allOnesDist = HMMTools.createUniformDistribution(alphabet,true);
            //     Distribution allOnesDist = new SimpleDistribution(alphabet);
            //     DistributionTools.allWeightsTo(allOnesDist, 1.0);

            //     // Add all spacer-lengths between Min and Max inclusive
            //     // 0-based indexing
            //     int spacerLenMin = (int)initFileCache.getInt("AFFINITY MODELS", "SpacerLengthMin");
            //     int spacerLenMax = (int)initFileCache.getInt("AFFINITY MODELS", "SpacerLengthMax");

            //     boolean getSpacerLengths = ((String)initFileCache.get("AFFINITY MODELS", "GetSpacerLengthAffinities")).equalsIgnoreCase("Yes");

            //     double[][][] newPsamsPosWeights = new double[psamsArray.length * (spacerLenMax - spacerLenMin + 2)][][];

            //     // loop for each psam
            //     for (int varNum = 0; varNum < psamsArray.length; varNum++) {

            //         WeightMatrix aPsam = psamsArray[varNum];
            //         //double[][] aPositionalWeights = psamsPosWeights[varNum];

            //         aPsam.setName("PSAM_"+varNum);

            //         // Don't send in a hamming sphere ever
            //         motifToHammingSphereMap.put(aPsam, null);

            //         int aMotifLength = aPsam.columns();

            //         if (getSpacerLengths) {

            //             int spacerPos = (int)Math.round(Math.floor((float) aMotifLength / 2));

            //             newPsamsPosWeights[varNum * (spacerLenMax - spacerLenMin + 2)] = psamsPosWeights[varNum];

            //             for (int spacerLen = spacerLenMin; spacerLen <= spacerLenMax; spacerLen++) {

            //                 newPsamsPosWeights[(varNum * (spacerLenMax - spacerLenMin + 2)) + (spacerLen - spacerLenMin +1)] = psamsPosWeights[varNum];

            //                 boolean[][] mandatoryColumns = new boolean[2][aMotifLength+spacerLen];
            //                 // include the positions neighboring the spacer as mandatory
            //                 Arrays.fill(mandatoryColumns[0], spacerPos-1, spacerPos+spacerLen+2, true);
            //                 mandatoryColumns[1] = ArrayTools.reverse(mandatoryColumns[0]);

            //                 // create uniform-dist spacer in the PSAM
            //                 WeightMatrix addedSpacerPSAM = WeightMatrixTools.insert(aPsam, spacerPos, spacerLen, allOnesDist);

            //                 // name has 1-based indexing
            //                 //addedSpacerPSAM.setName("[I("+ spacerPos +"-"+(spacerPos+1)+","+spacerLen+","+"N"+")]");
            //                 //addedSpacerPSAM.setName("insert_"+ spacerPos +"_"+(spacerPos+1)+"_"+spacerLen+"_"+"N"+"");
            //                 addedSpacerPSAM.setName("PSAM_"+varNum+"_insert_"+ spacerPos +"_"+(spacerPos+1)+"_"+spacerLen+"_"+"N"+"");

            //                 //                         out.println((byte)0, "**********************************************************");
            //                 //                         out.println((byte)0, "WeightMatrix: "+WeightMatrixTools.toString(addedSpacerPSAM));
            //                 //                         out.println((byte)0, "**********************************************************");

            //                 // Don't send in a hamming sphere ever
            //                 motifsToMandatoriesMap.put(addedSpacerPSAM, mandatoryColumns);
            //                 motifToHammingSphereMap.put(addedSpacerPSAM, null);

            //             } // end for (each spacer)

            //         } // if getSpacerLengths

            //     } // end for (each motif)

            //     if (getSpacerLengths) {
            //         psamsPosWeights = newPsamsPosWeights;
            //     }

            //     // Add Positions to the RHS
            //     //this.addWeightsToLeft = false;
            //     this.subtractOtherAffinities = true;

            //     // perform linear regression to find coefficients for each spacer-length
            //     regressionFitData = fitModel(
            //         motifToHammingSphereMap,
            //         null, //startPositions
            //         positionalWeights,
            //         null, //motifsToModsMap
            //         motifsToMandatoriesMap,
            //         "nucleotide-independent spacer-length features",
            //         //linearAlgorithm,
            //         spacerLenRegressionAlgorithm,
            //         effectiveStrand,
            //         strand,
            //         calc,
            //         null, //startParams
            //         null, //refSymList
            //         Double.NaN, //refKaReal);
            //         unifiedColumn);

            //     this.addWeightsToLeft = true;
            //     this.subtractOtherAffinities = false;

            //     // coefficients are as returned in R
            //     // First coefficient is the intercept, then the rest are in independent-var order
            //     out.print(regressionFitData.toString());

            //     double[] coefs = regressionFitData.getCoefficients();
            //     double[] coefsWithoutIntercept = Arrays.copyOfRange(coefs, 1, coefs.length);
            //     double maxCoef = MathTools.max(coefsWithoutIntercept);
            //     double[] relAffinities = MathTools.divide(coefsWithoutIntercept, maxCoef);

            //     out.println("\nRelative Affiniities are: "+StringTools.toString(relAffinities, ", "));

            // } //end - Perform regression to compare motifs


            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            // Hone the PSAM before extension
            ////////////////////////////////////////////////////////////////////////////////
            // Refit PSAM
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////

            int numHoneIterations = (int)initFileCache.getInt("REGRESSION", "NumHoneIterations");

            //for (int psamIter=1; psamIter <= 1; psamIter++) {
            //for (int psamIter=1; psamIter <= 0; psamIter++) {
            for (int psamIter=1; psamIter <= numHoneIterations; psamIter++) {

                hamming1PSAM = getPsam(
                    linearAlgorithm,
                    hamming1PSAM,
                    positionalWeights,
                    minAffinityWeight,
                    backgroundDist,
                    effectiveStrand,
                    strand,
                    calc,
                    eToMu,
                    refSymList,
                    referenceKa,
                    unifiedColumn,
                    alpha,
                    beta,
                    nonSpecKa);

//                 if (displayMotifs) {
//                     add(psamsFrame1, hamming1PSAM);
//                 }
            }

            if (saveAll) {
                WeightMatrixTools.writeToXML(hamming1PSAM, resultsDir + File.separator + proteinLabel+".psam."+motifLength+"nt.xml");
                FileTools.writeSerializedFile(hamming1PSAM, resultsDir + File.separator + proteinLabel+".psam."+motifLength+"nt.ser");

//                 if (displayMotifs) {
//                     GridBagLayoutFrame saveLogoFrame = display(proteinLabel+" PSAM before extension", hamming1PSAM);
//                     GraphicsTools.writeImage(saveLogoFrame.getPanel(), resultsDir + File.separator + proteinLabel+".psam."+motifLength+"nt.png", false);
//                 }

//                 FileTools.write(proteinLabel +"\t"+ seedSymList.seqString() +"\t"+ noPosWeightsRSquared +"\t"+ rSquared +"\n", resultsDir + File.separator + "unifiedMotifs."+motifLength+"nt.info", true);
//                 FileTools.write(seedSymList.seqString() +"\t"+ noPosWeightsRSquared +"\t"+ rSquared +"\n", resultsDir + File.separator + proteinLabel+"."+motifLength+"nt.info", false);

                //FileTools.write(""+rSquared, resultsDir + File.separator + proteinLabel+".final.rSquared", false);
                //GraphicsTools.writeImage(psamsFrame.getPanel(), resultsDir + File.separator + proteinLabel+".psams.png");
            }

            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            // Perform regression to extend the motif on one or both sides (if desired)
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            int seedMotifMinLeftExtensionLength = initFileCache.getInt("AFFINITY MODELS", "SeedMotifMinLeftExtensionLength");
            int seedMotifMaxLeftExtensionLength = initFileCache.getInt("AFFINITY MODELS", "SeedMotifMaxLeftExtensionLength");
            int seedMotifMinRightExtensionLength = initFileCache.getInt("AFFINITY MODELS", "SeedMotifMinRightExtensionLength");
            int seedMotifMaxRightExtensionLength = initFileCache.getInt("AFFINITY MODELS", "SeedMotifMaxRightExtensionLength");
            boolean fitAllPositionsWhenExtend = (boolean)initFileCache.getBoolean("REGRESSION", "FitAllPositionsWhenExtend");

            int extensionLength = 0;
            double lastRsquared = 0;

            //positionalWeights = null;
            double[][] rawWeights = null;
            double[][] cleanedWeights =null;
            double[][] smoothedWeights = null;

            int start;
            int increment;

            int extLeftIndex=0;
            int extRightIndex=0;

            // extend the PSAM on both sides seedMotifMaxExtensionLength times
            // then re-Train 3 times using the positionalWeights
            for (int iter=1; iter <= Math.max(seedMotifMaxLeftExtensionLength, seedMotifMaxRightExtensionLength); iter++) {
            //for (int iter=0; iter < seedMotifMaxExtensionLength+1; iter++) {
            //for (int iter=0; iter < seedMotifMaxExtensionLength+3; iter++) {

                // newPSAM grows as we add to the left-side then the right-side
                UniformDistribution uniformDist = new UniformDistribution(alphabet);
                WeightMatrix newPSAM = hamming1PSAM;
                double revCompSimilarity = WeightMatrixTools.revCompSimilarity(newPSAM);
                boolean[] positionsToFit = new boolean[hamming1PSAM.columns()];

                // Extend the PSAM and don't use the positionalWeights for training
                // extend the PSAM on both sides seedMotifMaxExtensionLength times
                if ((extLeftIndex < seedMotifMaxLeftExtensionLength) || (extRightIndex < seedMotifMaxRightExtensionLength)) {

                    // Perform minExtensionLength extensions immediately
                    if ((extLeftIndex < seedMotifMinLeftExtensionLength) || (extRightIndex < seedMotifMinRightExtensionLength)) {

                        while ((extLeftIndex < seedMotifMinLeftExtensionLength) || (extRightIndex < seedMotifMinRightExtensionLength)) {

                            if (extLeftIndex < seedMotifMinLeftExtensionLength) {

                                newPSAM = WeightMatrixTools.insert(
                                    newPSAM,
                                    0,
                                    1,
                                    uniformDist);

                                positionsToFit = ArrayUtils.add(positionsToFit, 0, true);

                                if (positionalWeights != null) {
                                    // add 0 to the end of the positionalWeights
                                    double[][] newPositionalWeights = new double[2][];
                                    newPositionalWeights[0] = ArrayUtils.add(positionalWeights[0], positionalWeights[0].length, 0);
                                    //newPositionalWeights[1] = ArrayUtils.add(positionalWeights[1], positionalWeights[1].length, 0);
                                    newPositionalWeights[1] = ArrayUtils.add(positionalWeights[1], 0, 0);
                                    positionalWeights = newPositionalWeights;

                                }
                                extLeftIndex++;
                            }
                            if (extRightIndex < seedMotifMinRightExtensionLength) {
                                newPSAM = WeightMatrixTools.insert(
                                    newPSAM,
                                    newPSAM.columns(),
                                    1,
                                    uniformDist);

                                    positionsToFit = ArrayUtils.add(positionsToFit, positionsToFit.length, true);

                                if (positionalWeights != null) {
                                    // add 0 to the beginning of the positionalWeights
                                    double[][] newPositionalWeights = new double[2][];
                                    newPositionalWeights[0] = ArrayUtils.add(positionalWeights[0], 0, 0);
                                    //newPositionalWeights[1] = ArrayUtils.add(positionalWeights[1], 0, 0);
                                    newPositionalWeights[1] = ArrayUtils.add(positionalWeights[1], positionalWeights[1].length, 0);
                                    positionalWeights = newPositionalWeights;
                                }
                                extRightIndex++;
                            }

                        }
                    }

                    // if over the minExtensionLength then perform one extension
                    else {

                        if (extLeftIndex < seedMotifMaxLeftExtensionLength) {
                            newPSAM = WeightMatrixTools.insert(
                                newPSAM,
                                0,
                                1,
                                uniformDist);

                            positionsToFit = ArrayUtils.add(positionsToFit, 0, true);

                            if (positionalWeights != null) {
                                // add 0 to the end of the positionalWeights
                                double[][] newPositionalWeights = new double[2][];
                                newPositionalWeights[0] = ArrayUtils.add(positionalWeights[0], positionalWeights[0].length, 0);
                                //newPositionalWeights[1] = ArrayUtils.add(positionalWeights[1], positionalWeights[1].length, 0);
                                newPositionalWeights[1] = ArrayUtils.add(positionalWeights[1], 0, 0);
                                positionalWeights = newPositionalWeights;
                            }
                            extLeftIndex++;

                        }
                        if (extRightIndex < seedMotifMaxRightExtensionLength) {
                            newPSAM = WeightMatrixTools.insert(
                                newPSAM,
                                newPSAM.columns(),
                                1,
                                uniformDist);

                            positionsToFit = ArrayUtils.add(positionsToFit, positionsToFit.length, true);

                            if (positionalWeights != null) {
                                // add 0 to the beginning of the positionalWeights
                                double[][] newPositionalWeights = new double[2][];
                                newPositionalWeights[0] = ArrayUtils.add(positionalWeights[0], 0, 0);
                                //newPositionalWeights[1] = ArrayUtils.add(positionalWeights[1], 0, 0);
                                newPositionalWeights[1] = ArrayUtils.add(positionalWeights[1], positionalWeights[1].length, 0);
                                positionalWeights = newPositionalWeights;
                            }
                            extRightIndex++;

                        }
                    }

                }

                // get new startPositions
                startPositions = getStartPositions(probeSeqLengths, newPSAM.columns());

//                 GridBagLayoutFrame psamsFrame2 = null;
//                 if (displayMotifs) {
//                     psamsFrame2 = new GridBagLayoutFrame(proteinLabel+": "+"ColumnREDUCE PSAMs after Extension "+extLeftIndex, true);
//                 }

//                 positionalWeights = new double[2][16 + 2*(newPSAM.columns()-1)];
//                 MathTools.fill(positionalWeights, 1.0);

//                 for (int ii = 0; ii < positionalWeights.length; ii++) {
//                     for (int jj = 0; jj < (newPSAM.columns()-3); jj++) {
//                         positionalWeights[ii][jj] = 0;
//                         positionalWeights[ii][positionalWeights[0].length - jj - 1] = 0;
//                     }
//                 }

                ////////////////////////////////////////////////////////////////////////////////
                // Refit PSAM
                ////////////////////////////////////////////////////////////////////////////////

                if (fitAllPositionsWhenExtend) {
                    positionsToFit = null;
                }

                for (int psamIter=1; psamIter <= 1; psamIter++) {
                    //for (int psamIter=1; psamIter <= 1; psamIter++) {
                    newPSAM = getPsam(
                        linearAlgorithm,
                        newPSAM,
                        positionalWeights,
                        minAffinityWeight,
                        backgroundDist,
                        effectiveStrand,
                        strand,
                        calc,
                        eToMu,
                        refSymList,
                        referenceKa,
                        unifiedColumn,
                        alpha,
                        beta,
                        nonSpecKa,
                        positionsToFit);

//                     if (displayMotifs) {
//                         add(psamsFrame2, newPSAM);
//                     }
                }

                newPSAM.setName("Extended"+extLeftIndex+"PSAM");

                // Get the fit for the new hamming1PSAM with no positional Weights
                regressionFitData = getRegressionFit(
                    newPSAM,
                    null,
                    "fit with the extended PSAM",
                    "lm",
                    effectiveStrand,
                    strand,
                    calc,
                    unifiedColumn);

                noPosWeightsRSquared = regressionFitData.rSquared;

                revCompSimilarity = WeightMatrixTools.revCompSimilarity(newPSAM);

                nonSpecKa.setValue(getRelativeKaNonSpec(
                    regressionFitData.getIntercept(),
                    regressionFitData.getCoefficient(newPSAM),
                    null,
                    eToMu,
                    probeSeqLengths,
                    motifLength,
                    strand,
                    calc,
                    revCompSimilarity));

                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                // Get Positional Bias for this strand
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                double posWeightsRSquared = 0;
                if (((String)initFileCache.get("AFFINITY MODELS", "GetPositionalBias")).equalsIgnoreCase("Yes")) {

                    // get occupancy bias
                    rawWeights = getPositionalBias(
                        //"nnls",
                        posBiasAlgorithm,
                        newPSAM,
                        calc,
                        startPositions,
                        effectiveStrand,
                        strand,
                        refSymList,
                        referenceKa,
                        unifiedColumn,
                        alpha,
                        beta,
                        nonSpecKa);

                    cleanedWeights = cleanWeights(startPositions, rawWeights);

                    smoothedWeights = smoothWeights(startPositions, cleanedWeights);
                    //smoothedWeights = smoothWeights(startPositions, rawWeights);

                    posBiasWeightsL1 = getPosBiasWeightsL1(smoothedWeights, cleanedWeights);

                    //if (((String)initFileCache.get("AFFINITY MODELS", "GetPositionalBias")).equalsIgnoreCase("Yes")) {

                    if (positionalBiasType.equalsIgnoreCase("Smoothed")) {
                        out.println("\nUsing smoothed Positional Bias Coefficients for normalization.");
                        positionalWeights = smoothedWeights;
                    }
                    else if (positionalBiasType.equalsIgnoreCase("Cleaned")) {
                        out.println("\nUsing cleaned Positional Bias Coefficients for normalization.");
                        positionalWeights = cleanedWeights;
                    }
                    else { // Use Raw
                        out.println("\nUsing raw Positional Bias Coefficients for normalization.");
                        positionalWeights = rawWeights;
                    }

                    if (kmerModelLengths != null) {
                        for(int z=0; z < kmerModelLengths.length; z++) {
                            if ((isLearnedPosBias[z]) && (motifLength == kmerModelLengths[z])) {
                                kmerModelWeights[z] = smoothedWeights;
                            }
                        }
                    }

                    out.println("\nThe new positional weights are:\n"+StringTools.toString(positionalWeights));
                    //}

                    // Get the fit for the new hamming1PSAM with the positional Weights
                    regressionFitData = getRegressionFit(
                        newPSAM,
                        positionalWeights,
                        "fit with the extended PSAM with weights",
                        "lm",
                        effectiveStrand,
                        strand,
                        calc,
                        unifiedColumn);

                    posWeightsRSquared = regressionFitData.rSquared;

                    revCompSimilarity = WeightMatrixTools.revCompSimilarity(newPSAM);

                    nonSpecKa.setValue(getRelativeKaNonSpec(
                            regressionFitData.getIntercept(),
                            regressionFitData.getCoefficient(newPSAM),
                            positionalWeights,
                            eToMu,
                            probeSeqLengths,
                            motifLength,
                            strand,
                            calc,
                            revCompSimilarity));

                }

                if (extLeftIndex < seedMotifMaxLeftExtensionLength+1) {
                    WeightMatrixTools.writeToXML(newPSAM, resultsDir + File.separator + "hammingPsamAfterExt"+extLeftIndex+".xml");
                }

                if (displayMotifs) {
                    BioJavaChart chart = null;
                    if (((String)initFileCache.get("AFFINITY MODELS", "GetPositionalBias")).equalsIgnoreCase("Yes")) {
                        String[] datasetLabels = {"Raw Positive Strand", "Raw Negative Strand", "Cleaned Positive Strand", "Cleaned Negative Strand", "Smoothed Positive Strand", "Smoothed Negative Strand"};
                        double[][] values = {ArrayTools.toDoubleArray(startPositions), rawWeights[0], rawWeights[1], cleanedWeights[0], cleanedWeights[1], smoothedWeights[0], smoothedWeights[1]};
                        chart = new BioJavaChart("lines", null, "Start Position", "Normalized Coefficients", null, datasetLabels, true, values, -1);
                    }
                    addMotif(strandLogoFrame, newPSAM, "Hamming PSAM", (extLeftIndex+0)+" + "+seedSymList.seqString()+" + "+(extLeftIndex+0), effectiveStrand.toString(), noPosWeightsRSquared, chart, posWeightsRSquared);
                }

                hamming1PSAM = newPSAM;
                motifLength = hamming1PSAM.columns();

//                 // If beyond minimum motifExtensionLength and the R2 has increased then continue
//                 //if ((extensionLength < seedMotifMinExtensionLength) || ((posWeightsRSquared * 1.0) > lastRsquared)) {
//                 if ((extensionLength < seedMotifMinExtensionLength) && !(extensionLength > seedMotifMaxExtensionLength)) {
//                     //out.println("extensionLength="+extensionLength+"  seedMotifMinExtensionLength="+seedMotifMinExtensionLength+"\n");
//                     extensionLength++;
//                     lastRsquared = posWeightsRSquared;
//                     hamming1PSAM = newPSAM;
//                     motifLength = hamming1PSAM.columns();
//                 }
//                 // else the motif has grown too long so stop
//                 else {
//                     break;
//                 }

//                 if (extLeftIndex < (seedMotifMaxExtensionLength-1)) {
//                     extLeftIndex++;
//                 }


                if (saveAll) {
                    WeightMatrixTools.writeToXML(hamming1PSAM, resultsDir + File.separator + proteinLabel+".psam."+motifLength+"nt.xml");
                    FileTools.writeSerializedFile(hamming1PSAM, resultsDir + File.separator + proteinLabel+".psam."+motifLength+"nt.ser");
                    //FileTools.writeSerializedFile(fsam, resultsDir + File.separator + proteinLabel+".fsam."+motifLength+"nt.ser");

                    if (displayMotifs) {

                        GraphicsTools.writeImage(strandLogoFrame.getPanel(), resultsDir + File.separator + proteinLabel+".results."+motifLength+"nt.png", false);

                        //                 GridBagLayoutFrame finalLogoFrame = new GridBagLayoutFrame(proteinLabel+" Final PSAM", true);
                        //                 WeightMatrixLogo finalLogo = new WeightMatrixLogo(hamming1PSAM, true, true, true, 0, WeightMatrixLogo.PositionLoc.BOTTOM, null);
                        //                 finalLogoFrame.add(finalLogo, 10, 10, false, true);
                        GridBagLayoutFrame saveLogoFrame = display(proteinLabel+" Final PSAM", hamming1PSAM);
                        GraphicsTools.writeImage(saveLogoFrame.getPanel(), resultsDir + File.separator + proteinLabel+".psam."+motifLength+"nt.png", false);

                    }

                    if (positionalWeights != null) {
                        String[] rowLabels = {"posStrandRelBias", "negStrandRelBias"};
                        String[] columnLabels = StringTools.toStrings(startPositions);
                        String outputTable = StringTools.toString(positionalWeights, rowLabels, columnLabels);
                        FileTools.write(
                            outputTable,
                            resultsDir + File.separator + proteinLabel+".psam."+motifLength+"nt.positionalBias.table",
                            false);
                    }

                    //GridBagLayoutFrame fsamLogo = display(proteinLabel+" FSAM Logo", fsam, featureThresh);
                    //GraphicsTools.writeImage(fsamLogo.getPanel(), resultsDir + File.separator + proteinLabel+".fsam."+motifLength+"nt.png", false);

                    FileTools.write(proteinLabel +"\t"+ seedSymList.seqString() +"\t"+ noPosWeightsRSquared +"\t"+ posWeightsRSquared +"\n", resultsDir + File.separator + "unifiedMotifs."+motifLength+"nt.psam.info", true);


                    FileTools.write(                    "seedMotif"+"\t"+ "noPosWeightsRSquared" +"\t"+ "posWeightsRSquared" +"\n", resultsDir + File.separator + proteinLabel+"."+motifLength+"nt.psam.info", false);
                    FileTools.write(                    seedSymList.seqString() +"\t"+ noPosWeightsRSquared +"\t"+ posWeightsRSquared +"\n", resultsDir + File.separator + proteinLabel+"."+motifLength+"nt.psam.info", true);

                    //FileTools.write(""+rSquared, resultsDir + File.separator + proteinLabel+".final.rSquared", false);
                    //GraphicsTools.writeImage(psamsFrame.getPanel(), resultsDir + File.separator + proteinLabel+".psams.png");
                }


            } // end - for (int iter=0; iter < seedMotifMaxExtensionLength+3; iter++) {


            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            // Hone the PSAM after extension
            ////////////////////////////////////////////////////////////////////////////////
            // Refit PSAM
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////

            int numHoneIterationsAfter = (int)initFileCache.getInt("REGRESSION", "NumHoneIterationsAfterExtend");

            //for (int psamIter=1; psamIter <= 1; psamIter++) {
            //for (int psamIter=1; psamIter <= 0; psamIter++) {
            for (int psamIter=1; psamIter <= numHoneIterationsAfter; psamIter++) {

                hamming1PSAM = getPsam(
                    linearAlgorithm,
                    hamming1PSAM,
                    positionalWeights,
                    minAffinityWeight,
                    backgroundDist,
                    effectiveStrand,
                    strand,
                    calc,
                    eToMu,
                    refSymList,
                    referenceKa,
                    unifiedColumn,
                    alpha,
                    beta,
                    nonSpecKa);

//                 if (displayMotifs) {
//                     add(psamsFrame1, hamming1PSAM);
//                 }
            }

            extensionLength = extLeftIndex;

            // out.println("\n*****************************************************************************************");
            // out.println("The R-squareds show that the seed PSAM should be extended "+extensionLength+" base pairs on "+seedMotifExtensionSide+" sides.");
            // out.println("******************************************************************************************");

            // save the weight matrix as PNG and XML
            //GraphicsTools.writeImage(logoPanel.getPanel(), resultsDir + File.separator + "hammingPsam.png");
            WeightMatrixTools.writeToXML(hamming1PSAM, resultsDir + File.separator + "hammingPsamAfterPosWeights.xml");

//             GridBagLayoutFrame FSAMFrame = null;
//             if (displayMotifs) {
//                 // create new FSAM logo
//                 FSAMFrame = createFSAMFrame(seedSymList.seqString(), strand);
//                 //addMotif(FSAMFrame, hamming1PSAM, rSquared);
//                 addMotif(FSAMFrame, hamming1PSAM);
//             }

//             hamming1PSAM = getPsam(linearAlgorithm, hamming1PSAM, null, minAffinityWeight, backgroundDist, strand, calc, eToMu);
//             logoFrame = display("PSAM after regression on every column (after extension)", hamming1PSAM);

//             hamming1PSAM = getPsam(linearAlgorithm, hamming1PSAM, null, minAffinityWeight, backgroundDist, strand, calc, eToMu);
//             logoFrame = display("PSAM after second regression on every column (after extension)", hamming1PSAM);

            ////////////////////////////////////////////////////////////////////////////////
            // Loop through
            //   1. Fitting for probe-position bias
            //   2. Fitting for psam-column specificity
            //   3. Fitting for features
            // until R^2 no longer increases
            ////////////////////////////////////////////////////////////////////////////////

            FeaturedWeightMatrix fsam = new FeaturedWeightMatrix("FSAM", hamming1PSAM, strand, calc, complementTable);
            fsam.setPositionalWeights(positionalWeights);
            fsam.kmerLengths = kmerModelLengths;
            fsam.kmerPositionalWeights = kmerModelWeights;

            LinkedHashMap<Object, Symbol[]> motifsToModsMap;
            //LinkedHashMap<Object, boolean[][]> motifsToMandatoriesMap;
            ArrayList<FeatureKey> aFeatureRemovalList = null;

            int minIterations = (int)initFileCache.getInt("AFFINITY MODELS", "MinIterations");
            startPositions = getStartPositions(probeSeqLengths, motifLength);

            for (int iterNum = 1; iterNum <= minIterations; iterNum++) {


                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                // Perform regression to find any positional bias
                // if positionalWeights haven't already been calculated for iter 1
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
//                 if (((positionalWeights == null) || (iterNum > 1 ))
//                     && (((String)initFileCache.get("AFFINITY MODELS", "GetPositionalBias")).equalsIgnoreCase("Yes"))) {

                // if (((positionalWeights == null) && (iterNum == 1 ))
                //     && (((String)initFileCache.get("AFFINITY MODELS", "GetPositionalBias")).equalsIgnoreCase("Yes"))) {

                if (((true) && (iterNum == 1 ))
                    && (((String)initFileCache.get("AFFINITY MODELS", "GetPositionalBias")).equalsIgnoreCase("Yes"))) {

                // get gamma biases (strand bias and local concentration scalars)
                    rawWeights = getPositionalBias(
                        //"nnls",
                        posBiasAlgorithm,
                        fsam,
                        calc,
                        startPositions,
                        effectiveStrand,
                        strand,
                        refSymList,
                        referenceKa,
                        unifiedColumn,
                        alpha,
                        beta,
                        nonSpecKa);

                    cleanedWeights = cleanWeights(startPositions, rawWeights);

                    smoothedWeights = smoothWeights(startPositions, cleanedWeights);
                    //smoothedWeights = smoothWeights(startPositions, rawWeights);

                    posBiasWeightsL1 = getPosBiasWeightsL1(smoothedWeights, cleanedWeights);

                    if (displayMotifs) {
                        String[] datasetLabels = {"Raw Positive Strand", "Raw Negative Strand", "Cleaned Positive Strand", "Cleaned Negative Strand", "Smoothed Positive Strand", "Smoothed Negative Strand"};
                        double[][] values = {ArrayTools.toDoubleArray(startPositions), rawWeights[0], rawWeights[1], cleanedWeights[0], cleanedWeights[1], smoothedWeights[0], smoothedWeights[1]};
                        //BioJavaChart chart = new BioJavaChart("lines", "Raw Weights", "Start Position", "Normalized Coefficients", null, datasetLabels, true, values, -1);
                        logoFrame = display(proteinLabel+": "+"Probe Bias Regression Results for iteration "+iterNum, "lines", "Positional Weights", "Start Position", "Normalized Coefficients", null, datasetLabels, true, values, -1);
                        //logoFrame = display("lines", "Positional Weights", "Start Position", "Normalized Coefficients", null, null, false, values, -1);

                        //BioJavaChart chart = new BioJavaChart("lines", null, "Start Position", "Normalized Coefficients", null, datasetLabels, true, values, -1);
                        //addMotif(strandLogoFrame, fsame.getPosStrandPWM(), "Hamming PSAM", (extLeftIndex+0)+" + "+seedSymList.seqString()+" + "+(extLeftIndex+0), effectiveStrand.toString(), rSquared, chart, posWeightsRSquared);
                    }

                    if (positionalBiasType.equalsIgnoreCase("Smoothed")) {
                        out.println("\nUsing smoothed Positional Bias Coefficients for normalization.");
                        positionalWeights = smoothedWeights;
                    }
                    else if (positionalBiasType.equalsIgnoreCase("Cleaned")) {
                        out.println("\nUsing cleaned Positional Bias Coefficients for normalization.");
                        positionalWeights = cleanedWeights;
                    }
                    else { // Use Raw
                        out.println("\nUsing raw Positional Bias Coefficients for normalization.");
                        positionalWeights = rawWeights;
                    }
                    out.println("\nThe new positional weights are:\n"+StringTools.toString(positionalWeights));

                    fsam.setPositionalWeights(positionalWeights);

                    if (kmerModelLengths != null) {
                        for(int z=0; z < kmerModelLengths.length; z++) {
                            if ((isLearnedPosBias[z]) && (motifLength == kmerModelLengths[z])) {
                                fsam.kmerPositionalWeights[z] = smoothedWeights;
                            }
                        }
                    }

                    String psamAlgorithm;
                    if ((nonlinearAlgorithm == null) || nonlinearAlgorithm.equals("")) {
                        psamAlgorithm = linearAlgorithm;
                    }
                    else {
                        psamAlgorithm = nonlinearAlgorithm;
                    }

                    ////////////////////////////////////////////////////////////////////////////////
                    // Perform regression to calculate a new PSAM using the positional bias
                    ////////////////////////////////////////////////////////////////////////////////
                    GridBagLayoutFrame psamsFrame = null;
                    if (displayMotifs) {
                        psamsFrame = new GridBagLayoutFrame(proteinLabel+": "+"PSAMs after iteration "+iterNum, true);
                    }

                    int numPosBiasIterations = (int)initFileCache.getInt("REGRESSION", "NumPosBiasIterations");

                    for (int psamIter=1; psamIter <= numPosBiasIterations; psamIter++) {

                        hamming1PSAM = getPsam(
                            psamAlgorithm,
                            fsam.getPosStrandPWM(),
                            positionalWeights,
                            minAffinityWeight,
                            backgroundDist,
                            effectiveStrand,
                            strand,
                            calc,
                            eToMu,
                            refSymList,
                            referenceKa,
                            unifiedColumn,
                            alpha,
                            beta,
                            nonSpecKa);

                        fsam.putPWMs(hamming1PSAM, complementTable);

                        if (displayMotifs) {
                            add(psamsFrame, hamming1PSAM);
                        }
                    }

                    if((boolean)initFileCache.getBoolean("AFFINITY MODELS", "RemoveProbes")) {

                        int truncatePos = (int)initFileCache.getInt("AFFINITY MODELS", "TruncationPosition");
                        out.println("Truncation Position to remove leading nucleotide positions is "+truncatePos+".");

                        double[][] removalWeights = getRemovalWeights(startPositions, positionalWeights);
                        //positionalWeights = removalWeights;

                        if (removalWeights != null) {
                            String[] rowLabels = {"posStrandRemovalWeights", "negStrandRemovalWeights"};
                            String[] columnLabels = StringTools.toStrings(startPositions);
                            String outputTable = StringTools.toString(removalWeights, rowLabels, columnLabels);
                            out.println("\n"+outputTable);
                        }

                        // get relative affinities for all the probes using the removalWeights
                        out.print("Retrieving the counts with the current PSAM using the "+effectiveStrand.toString()+" strand ...");

                        double[] motifCounts =  reduceData.getCountsArray(
                            hamming1PSAM,
                            null, //hammingSphere,
                            null, //startPosition
                            null, //mandatoryColumns,
                            removalWeights,
                            strand,
                            calc,
                            1, //experimentNumber,
                            unifiedColumn);

                        out.println(" Done. The sum is "+MathTools.sum(motifCounts)+".");

                        // remove probes above the threshold
                        double threshold = (double)initFileCache.getDouble("AFFINITY MODELS", "ProbeRemovalThresh");
                        removeProbesArray = MathTools.isGreaterThan(motifCounts, threshold);
                        int numTrue = MathTools.sum(true, removeProbesArray);
                        out.println(numTrue+" probes have a relative affinity above the threshold of "+threshold+".");

                        ///////////////////////////////////////////////////////////////////////////////////
                        // refit a new PSAM using removed weights and removed probes
                        ///////////////////////////////////////////////////////////////////////////////////
                        double[][] removedWeights = getRemovedWeights(startPositions, positionalWeights);
                        positionalWeights = removedWeights;

                        int numRemoveProbeIterations = (int)initFileCache.getInt("REGRESSION", "NumRemoveProbeIterations");

                        for (int psamIter=1; psamIter <= numRemoveProbeIterations; psamIter++) {

                            hamming1PSAM = getPsam(
                                psamAlgorithm,
                                fsam.getPosStrandPWM(),
                                positionalWeights,
                                //removedWeights,
                                minAffinityWeight,
                                backgroundDist,
                                effectiveStrand,
                                strand,
                                calc,
                                eToMu,
                                refSymList,
                                referenceKa,
                                unifiedColumn,
                                alpha,
                                beta,
                                nonSpecKa);

                            fsam.putPWMs(hamming1PSAM, complementTable);

                        }

                        // Get the fit for the new hamming1PSAM with no positional Weights
                        regressionFitData = getRegressionFit(
                            hamming1PSAM,
                            null,
                            "fit with the extended PSAM",
                            "lm",
                            effectiveStrand,
                            strand,
                            calc,
                            unifiedColumn);

                        noPosWeightsRSquared = regressionFitData.rSquared;

                        // Get the fit for the new hamming1PSAM with the positional Weights
                        regressionFitData = getRegressionFit(
                            hamming1PSAM,
                            positionalWeights,
                            "fit with the extended PSAM with weights",
                            "lm",
                            effectiveStrand,
                            strand,
                            calc,
                            unifiedColumn);

                        double posWeightsRSquared = regressionFitData.rSquared;

                        double revCompSimilarity = WeightMatrixTools.revCompSimilarity(hamming1PSAM);

                        nonSpecKa.setValue(getRelativeKaNonSpec(
                                regressionFitData.getIntercept(),
                                regressionFitData.getCoefficient(hamming1PSAM),
                                positionalWeights,
                                eToMu,
                                probeSeqLengths,
                                motifLength,
                                strand,
                                calc,
                                revCompSimilarity));

                        if (displayMotifs) {
                            BioJavaChart chart = null;
                            if (((String)initFileCache.get("AFFINITY MODELS", "GetPositionalBias")).equalsIgnoreCase("Yes")) {
//                             String[] datasetLabels = {"Raw Positive Strand", "Raw Negative Strand", "Cleaned Positive Strand", "Cleaned Negative Strand", "Smoothed Positive Strand", "Smoothed Negative Strand"};
//                             double[][] values = {ArrayTools.toDoubleArray(startPositions), rawWeights[0], rawWeights[1], cleanedWeights[0], cleanedWeights[1], smoothedWeights[0], smoothedWeights[1]};

//                                 String[] datasetLabels = {"Smoothed Positive Strand", "Smoothed Negative Strand", "Removal Positive Strand", "Removal Negative Strand", "Removed Positive Strand", "Removed Negative Strand"};
//                                 double[][] values = {ArrayTools.toDoubleArray(startPositions), smoothedWeights[0], smoothedWeights[1], removalWeights[0], removalWeights[1], removedWeights[0], removedWeights[1]};

                                String[] datasetLabels = {"Raw Positive Strand", "Raw Negative Strand", "Cleaned Positive Strand", "Cleaned Negative Strand", "Smoothed Positive Strand", "Smoothed Negative Strand", "Removal Positive Strand", "Removal Negative Strand", "Removed Positive Strand", "Removed Negative Strand"};
                                double[][] values = {ArrayTools.toDoubleArray(startPositions), rawWeights[0], rawWeights[1], cleanedWeights[0], cleanedWeights[1], smoothedWeights[0], smoothedWeights[1], removalWeights[0], removalWeights[1], removedWeights[0], removedWeights[1]};

                                chart = new BioJavaChart("lines", null, "Start Position", "Normalized Coefficients", null, datasetLabels, true, values, -1);
                            }
                            addMotif(strandLogoFrame, hamming1PSAM, "Hamming PSAM", (extLeftIndex+0)+" + "+seedSymList.seqString()+" + "+(extLeftIndex+0), effectiveStrand.toString(), noPosWeightsRSquared, chart, posWeightsRSquared);
                        }

                    } // end - RemoveProbes


                } // end - Perform regression to find any positional bias


                ////////////////////////////////////////////////////////////////////////////////////////////////
                // Remove all features from the current FSAM!!
                // Set Non-specific Ka to 0
                ////////////////////////////////////////////////////////////////////////////////////////////////
                fsam.removeFeatures();
                fsam.setNonSpecKa(0);
                hamming1PSAM = fsam.getPosStrandPWM();

                ////////////////////////////////////////////////////////////////////////////////////////////////
                // Perform regression for current fsam!!
                ////////////////////////////////////////////////////////////////////////////////////////////////
                regressionFitData = getRegressionFit(
                    fsam,
                    positionalWeights,
                    "fit with the current FSAM",
                    "lm",
                    //featureDiscoveryAlgorithm,
                    effectiveStrand,
                    strand,
                    calc,
                    unifiedColumn);

                rSquared = regressionFitData.rSquared;
                this.fsamUnivarCoeff = regressionFitData.getCoefficient(fsam);
                this.fsamUnivarInter = regressionFitData.getIntercept();

                nonSpecKa.setValue(getRelativeKaNonSpec(
                    regressionFitData.getIntercept(),
                    regressionFitData.getCoefficient(fsam),
                    positionalWeights,
                    eToMu,
                    probeSeqLengths,
                    motifLength,
                    strand,
                    calc,
                    fsam.getRevCompSimilarity()));

                ////////////////////////////////////////////////////////////////////////////////
                // Perform regression to find the PSAM and Kmer fits for the different Hamming Spheres
                ////////////////////////////////////////////////////////////////////////////////
                if (((String)initFileCache.get("AFFINITY MODELS", "GetHammingSpheresAffinities")).equalsIgnoreCase("Yes")) {

                    ////////////////////////////////////////////////////////////////////////////////
                    // Get the Hamming-spheres
                    ////////////////////////////////////////////////////////////////////////////////
                    // 0 = seed alone
                    // 1 = Hamming-1 sphere
                    // 2 = hamming-2 sphere .....

                    // seedSymList needs to change if we extended the original PSAM
                    if (extensionLength > 0) {
                        seedSymList = WeightMatrixTools.getMostLikelySymList(hamming1PSAM);
                    }

                    // hammingDist0Set contains the seed symbolList alone
                    LinkedHashSet<SymbolList> hammingDist0Set = new LinkedHashSet<SymbolList>(1);
                    hammingDist0Set.add(seedSymList);
                    hammingSpheres.add(hammingDist0Set);

                    // add the hammingDist1Set already built
                    if (extensionLength == 0) {
                        if (includeRevComps && (hammingDist1WithRevCompsSet != null)) {
                            hammingSpheres.add(hammingDist1WithRevCompsSet);
                        }
                        else if (!includeRevComps && (hammingDist1WithoutRevCompsSet != null)) {
                            hammingSpheres.add(hammingDist1WithoutRevCompsSet);
                        }
                        else {
                            hammingSpheres.add(SymbolListTools.getHammingSphere(seedSymList, 1, includeRevComps, complementTable));
                        }
                    }
                    else {
                        hammingSpheres.add(SymbolListTools.getHammingSphere(seedSymList, 1, includeRevComps, complementTable));
                    }

                    // build the larger hamming spheres
                    for (int j=2; j <= 5; j++) {
                        out.print("\n\nGenerating the Hamming-"+j+" Sphere about the seed... ");
                        LinkedHashSet<SymbolList> smallerSphere = hammingSpheres.get(j-1);
                        LinkedHashSet<SymbolList> largerSphere = new LinkedHashSet<SymbolList>(smallerSphere);
                        SymbolListTools.getHammingSphere(largerSphere, 1, includeRevComps, complementTable);
                        hammingSpheres.add(largerSphere);
                        out.println(" Done.");
                    }

                    RegressionFitData regressionFitData2 = null;
                    ArrayList<String> outputStrings = new ArrayList<String>();
                    outputStrings.add("# PSAM\tKmers");

                    out.println("\n*******************************************************************************");
                    for (int k=1; k <= 3; k++) {

                        out.println("\nHamming-"+k+" Sphere contains "+hammingSpheres.get(k).size()+" Kmers.");

                        // Fit the PSAM
                        motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();
                        motifToHammingSphereMap.put(hamming1PSAM, hammingSpheres.get(k));

                        regressionFitData = fitModel(
                            motifToHammingSphereMap,
                            null, //startPositions
                            positionalWeights,
                            null, //aMotifsToModsMap
                            null, //aMotifsToMandatoriesMap
                            "fit of the PSAM using only the Hamming-"+k+" Sphere",
                            linearAlgorithm,
                            effectiveStrand,
                            strand,
                            calc,
                            null, //startParams
                            null, //refSymList
                            Double.NaN, //refKaReal
                            unifiedColumn);

                        out.println("\nThe R-squared for the PSAM using only the Hamming-"+k+" Sphere is "+regressionFitData.rSquared+".");

                        // Fit the Kmers
                        if (k <= 2) {
                            regressionFitData2 = fitModel(
                                hammingSpheres.get(k),
                                null, //startPositions
                                positionalWeights,
                                "fit of the Kmers in the Hamming-"+k+" Sphere",
                                linearAlgorithm,
                                effectiveStrand,
                                strand,
                                calc,
                                null, //startParams
                                null, //refSymList
                                Double.NaN, //refKaReal
                                unifiedColumn);

                            out.println("\nThe R-squared for the Kmers in the Hamming-"+k+" Sphere is "+regressionFitData2.rSquared+".");

                            outputStrings.add(regressionFitData.rSquared +"\t"+ regressionFitData2.rSquared);
                        }
                        else {
                            outputStrings.add(regressionFitData.rSquared +"\t"+"0.0");
                        }
                    }
                    out.println("\n*******************************************************************************");

                    FileTools.write((String[])outputStrings.toArray(new String[0]), "psam.kmers.fits.byHammingDist", false);

                } // end - Perform regression to find the PSAM and Kmer fits for the different Hamming Spheres


                // start caching counts arrays
                if (cacheCounts) {
                    reduceData.setCacheCountsArraysFlag(true);
                }

                ////////////////////////////////////////////////////////////////////////////////
                // Prepare to perform regressions for features
                ////////////////////////////////////////////////////////////////////////////////

                ////////////////////////////////////////////////////////////////////////////////
                // Perform regression to find context-independent, inserted nucleotide-independent insertion features
                ////////////////////////////////////////////////////////////////////////////////

                if (((String)initFileCache.get("AFFINITY MODELS", "Get_CINI_InsertionAffinities")).equalsIgnoreCase("Yes")) {

                    motifsToModsMap = new LinkedHashMap<Object, Symbol[]>();
                    motifsToMandatoriesMap = new LinkedHashMap<Object, boolean[][]>();
                    motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();

                    if (useHammingDistForApprox) {
                        motifToHammingSphereMap.put(hamming1PSAM, hammingSpheres.get(hammingDistForApprox));
                    }
                    else {
                        motifToHammingSphereMap.put(hamming1PSAM, null);
                    }

                    //Distribution allOnesDist = HMMTools.createUniformDistribution(alphabet,true);
                    Distribution allOnesDist = new SimpleDistribution(alphabet);
                    DistributionTools.allWeightsTo(allOnesDist, 1.0);

                    // Add all possible single-point insertions to the motifToHammingSphereMap
                    // DO NOT insert at locations 0-1, last-(last+1)
                    // 0-based indexing
                    for (int insert=1+extensionLength; insert < motifLength-extensionLength; insert++) {
                        //for (int insert=1; insert < motifLength; insert+=motifLength-2) {

                        // create uniform dist insertion in the PSAM
                        WeightMatrix singleInsertionPSAM = WeightMatrixTools.insert(hamming1PSAM, insert, 1, allOnesDist);

                        // name has 1-based indexing
                        //singleInsertionPSAM.setName("[I("+ insert +"-"+(insert+1)+","+1+","+"N"+")]");
                        singleInsertionPSAM.setName("insert_"+ insert +"_"+(insert+1)+"_"+1+"_"+"N"+"");
                        Symbol[] modsArray = new Symbol[(motifLength*2)-1];
                        modsArray[(insert*2)-1] = DNATools.n();

                        boolean[][] mandatoryColumns = new boolean[2][motifLength+1];
                        mandatoryColumns[0][insert] = true;
                        mandatoryColumns[0][insert-1] = true;
                        mandatoryColumns[0][insert+1] = true;
                        mandatoryColumns[1] = ArrayTools.reverse(mandatoryColumns[0]);

                        LinkedHashSet<SymbolList> hammingSphere = null;

                        if (useHammingDistForApprox) {

                            // hammingSphere = hammingSpheres.get(hammingDistForApprox);

                            // create the insertion in the seedSymList
                            // 1-based indexing
                            // insert an A, C, G, T
                            hammingSphere = new LinkedHashSet<SymbolList>(4);

                            // iterate over every possible insert symbol
                            for (Iterator dnaIter2 = ((FiniteAlphabet)alphabet).iterator(); dnaIter2.hasNext(); ) {
                                Symbol insertSymbol = (Symbol)dnaIter2.next();
                                Symbol[] symbolArray = {insertSymbol};

                                SimpleSymbolList editedSeed = new SimpleSymbolList(seedSymList);
                                Edit edit = new Edit(insert+1, 0, new SimpleSymbolList(symbolArray, 1, alphabet));
                                editedSeed.edit(edit);
                                hammingSphere.add(editedSeed);
                            }

                            // create the hamming sphere around the 4 inserts
                            //LinkedHashSet<SymbolList> hammingSphere = null;
                            SymbolListTools.getHammingSphere(hammingSphere, hammingDistForApprox, includeRevComps, complementTable);
                            out.println((byte)0, "hammingSphere = "+symListCollectionDump(hammingSphere));
                        }

                        //out.println((byte)0, "WeightMatrix: "+StringTools.toString(singleInsertionPSAM));
                        out.println((byte)0, "WeightMatrix: "+WeightMatrixTools.toString(singleInsertionPSAM));
                        out.println((byte)0, "**********************************************************");

                        motifToHammingSphereMap.put(singleInsertionPSAM, hammingSphere);
                        motifsToModsMap.put(singleInsertionPSAM, modsArray);
                        motifsToMandatoriesMap.put(singleInsertionPSAM, mandatoryColumns);
                    }

                    // perform penalized regression to find significant features
                    regressionFitData = fitModel(
                        motifToHammingSphereMap,
                        null, //startPositions
                        positionalWeights,
                        motifsToModsMap,
                        motifsToMandatoriesMap,
                        "context-independent, nucleotide-independent insertion features",
                        featureDiscoveryAlgorithm,
                        effectiveStrand,
                        strand,
                        calc,
                        null, //startParams
                        null, //refSymList
                        Double.NaN, //refKaReal);
                        unifiedColumn);

                    // coefficients are as returned in R
                    // First coefficient is the intercept, then the rest are in independent-var order
                    out.print(regressionFitData.toString());

                    fsam.addFeatures(
                        !regressionFitData.isCoefficientZero(hamming1PSAM),
                        regressionFitData.getFeatures(hamming1PSAM, featuresScalar));
                    //fsam.setIntercept(regressionFitData.getIntercept(hamming1PSAM));
                } // end - Perform regression to find context-independent, inserted nucleotide-independent insertion features


                ////////////////////////////////////////////////////////////////////////////////
                // Perform regression to find context-independent, inserted nucleotide-dependent insertion features
                ////////////////////////////////////////////////////////////////////////////////
                if (((String)initFileCache.get("AFFINITY MODELS", "Get_CIND_InsertionAffinities")).equalsIgnoreCase("Yes")) {

                    motifsToModsMap = new LinkedHashMap<Object, Symbol[]>();
                    motifsToMandatoriesMap = new LinkedHashMap<Object, boolean[][]>();
                    motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();

                    if (useHammingDistForApprox) {
                        motifToHammingSphereMap.put(hamming1PSAM, hammingSpheres.get(hammingDistForApprox));
                    }
                    else {
                        motifToHammingSphereMap.put(hamming1PSAM, null);
                    }

                    // Add all possible single-point insertions to the motifToHammingSphereMap
                    // DO NOT insert at locations 0-1, last-(last+1)
                    // 0-based indexing
                    for (int insert=1+extensionLength; insert < motifLength-extensionLength; insert++) {
                        //for (int insert=1; insert < motifLength; insert+=motifLength-2) {

                        boolean[][] mandatoryColumns = new boolean[2][motifLength+1];
                        mandatoryColumns[0][insert] = true;
                        mandatoryColumns[0][insert-1] = true;
                        mandatoryColumns[0][insert+1] = true;
                        mandatoryColumns[1] = ArrayTools.reverse(mandatoryColumns[0]);

                        // iterate over every possible insert symbol
                        for (Iterator dnaIter2 = ((FiniteAlphabet)alphabet).iterator(); dnaIter2.hasNext(); ) {
                            Symbol insertSymbol = (Symbol)dnaIter2.next();

                            // create insertion in the PSAM
                            WeightMatrix singleInsertionPSAM = WeightMatrixTools.insert(hamming1PSAM, insert, insertSymbol);

                            // name has 1-based indexing
                            //singleInsertionPSAM.setName("[I("+ insert +"-"+(insert+1)+","+1+","+insertSymbol.getName().substring(0,1).toUpperCase()+")]");
                            singleInsertionPSAM.setName("insert_"+ insert +"_"+(insert+1)+"_"+1+"_"+insertSymbol.getName().substring(0,1).toUpperCase()+"");
                            Symbol[] modsArray = new Symbol[(motifLength*2)-1];
                            modsArray[(insert*2)-1] = insertSymbol;

                            LinkedHashSet<SymbolList> hammingSphere = null;

                            if (useHammingDistForApprox) {

                                // hammingSphere = hammingSpheres.get(hammingDistForApprox);

                                // create the insertion in the seedSymList
                                // 1-based indexing
                                Symbol[] symbolArray = {insertSymbol};
                                SimpleSymbolList editedSeed = new SimpleSymbolList(seedSymList);
                                Edit edit = new Edit(insert+1, 0, new SimpleSymbolList(symbolArray, 1, alphabet));
                                editedSeed.edit(edit);

                                // create the hamming sphere around the insertion
                                hammingSphere = SymbolListTools.getHammingSphere(editedSeed, hammingDistForApprox, includeRevComps, complementTable);
                                out.println((byte)0, "hammingSphere = "+symListCollectionDump(hammingSphere));
                            }

                            out.println((byte)0, "WeightMatrix: "+WeightMatrixTools.toString(singleInsertionPSAM));
                            out.println((byte)0, "**********************************************************");

                            motifToHammingSphereMap.put(singleInsertionPSAM, hammingSphere);
                            motifsToModsMap.put(singleInsertionPSAM, modsArray);
                            motifsToMandatoriesMap.put(singleInsertionPSAM, mandatoryColumns);

                        }
                    }

                    // perform penalized regression to find significant features
                    regressionFitData = fitModel(
                        motifToHammingSphereMap,
                        null, //startPositions
                        positionalWeights,
                        motifsToModsMap,
                        motifsToMandatoriesMap,
                        "context-independent, nucleotide-dependent insertion features",
                        featureDiscoveryAlgorithm,
                        effectiveStrand,
                        strand,
                        calc,
                        null, //startParams
                        null, //refSymList
                        Double.NaN, //refKaReal);
                        unifiedColumn);

                    fsam.addFeatures(
                        !regressionFitData.isCoefficientZero(hamming1PSAM),
                        regressionFitData.getFeatures(hamming1PSAM, featuresScalar));
                    //fsam.setIntercept(regressionFitData.getIntercept(hamming1PSAM));

                    // coefficients are as returned in R
                    // First coefficient is the intercept, then the rest are in independent-var order
                    out.print(regressionFitData.toString());

                } // end - Perform regression to find context-independent, inserted nucleotide-dependent insertion features


                ////////////////////////////////////////////////////////////////////////////////
                // Perform regression to find context-dependent, inserted nucleotide-dependent insertion features
                ////////////////////////////////////////////////////////////////////////////////
                if (((String)initFileCache.get("AFFINITY MODELS", "Get_CDND_InsertionAffinities")).equalsIgnoreCase("Yes")) {

                    motifsToModsMap = new LinkedHashMap<Object, Symbol[]>();
                    motifsToMandatoriesMap = new LinkedHashMap<Object, boolean[][]>();
                    motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();

                    if (useHammingDistForApprox) {
                        motifToHammingSphereMap.put(hamming1PSAM, hammingSpheres.get(hammingDistForApprox));
                    }
                    else {
                        motifToHammingSphereMap.put(hamming1PSAM, null);
                    }

                    // Add all possible single-point insertions to the motifToHammingSphereMap
                    // DO NOT insert at locations 1-2, last-(last+1)
                    // 0-based indexing
                    for (int insert=1+extensionLength; insert < motifLength-extensionLength; insert++) {
                        //for (int insert=1; insert < motifLength; insert+=motifLength-2) {

                        boolean[][] mandatoryColumns = new boolean[2][motifLength+1];
                        mandatoryColumns[0][insert] = true;
                        mandatoryColumns[0][insert-1] = true;
                        mandatoryColumns[0][insert+1] = true;
                        mandatoryColumns[1] = ArrayTools.reverse(mandatoryColumns[0]);

                        // iterate over every possible insert symbol
                        for (Iterator dnaIter2 = ((FiniteAlphabet)alphabet).iterator(); dnaIter2.hasNext(); ) {
                            Symbol insertSymbol = (Symbol)dnaIter2.next();

                            // iterate over every symbol at the mut1 position (to the left of the insert)
                            for (Iterator dnaIter3 = ((FiniteAlphabet)alphabet).iterator(); dnaIter3.hasNext(); ) {
                                Symbol mut1Symbol = (Symbol)dnaIter3.next();

                                // iterate over every symbol at the mut2 position (to the right of the insert)
                                for (Iterator dnaIter4 = ((FiniteAlphabet)alphabet).iterator(); dnaIter4.hasNext(); ) {
                                    Symbol mut2Symbol = (Symbol)dnaIter4.next();

                                    // create insertion in the PSAM
                                    WeightMatrix contextInsertionPSAM = WeightMatrixTools.insert(hamming1PSAM, insert, insertSymbol);
                                    contextInsertionPSAM = WeightMatrixTools.mutate(contextInsertionPSAM, insert-1, mut1Symbol, false);
                                    contextInsertionPSAM = WeightMatrixTools.mutate(contextInsertionPSAM, insert+1, mut2Symbol, false);

                                    // name has 1-based indexing
                                    contextInsertionPSAM.setName(
                                        //                                         "["+(insert)+"("+mut1Symbol.getName().substring(0,1).toUpperCase()
                                        //                                         +"), I("+ insert +"-"+(insert+1)+","+1+","+insertSymbol.getName().substring(0,1).toUpperCase()
                                        //                                         +"), "+(insert+1)+"("+mut2Symbol.getName().substring(0,1).toUpperCase()+")]");
                                        "mutate_"+insert+"_"+mut1Symbol.getName().substring(0,1).toUpperCase()
                                        +"__"+"insert_"+ insert +"_"+(insert+1)+"_"+1+"_"+insertSymbol.getName().substring(0,1).toUpperCase()
                                        +"__"+"mutate_"+(insert+1)+"_"+mut2Symbol.getName().substring(0,1).toUpperCase()+"");

                                    Symbol[] modsArray = new Symbol[(motifLength*2)-1];
                                    modsArray[(insert*2)-2] = mut1Symbol;
                                    modsArray[(insert*2)-1] = insertSymbol;
                                    modsArray[(insert*2)] = mut2Symbol;

                                    LinkedHashSet<SymbolList> hammingSphere = null;

                                    if (useHammingDistForApprox) {

                                        //hammingSphere = hammingSpheres.get(hammingDistForApprox);

                                        // create the insertion in the seedSymList
                                        // 1-based indexing
                                        Symbol[] symbolArray = {insertSymbol};
                                        SimpleSymbolList editedSeed = new SimpleSymbolList(seedSymList);
                                        Edit edit = new Edit(insert+1, 0, new SimpleSymbolList(symbolArray, 1, alphabet));
                                        editedSeed.edit(edit);
                                        edit = new Edit((insert+1)-1, alphabet, mut1Symbol);
                                        editedSeed.edit(edit);
                                        edit = new Edit((insert+1)+1, alphabet, mut2Symbol);
                                        editedSeed.edit(edit);

                                        // create the hamming sphere around the mutation
                                        hammingSphere = SymbolListTools.getHammingSphere(editedSeed, hammingDistForApprox, includeRevComps, complementTable);
                                        out.println((byte)0, "hammingSphere = "+symListCollectionDump(hammingSphere));
                                    }

                                    out.println((byte)0, "WeightMatrix: "+WeightMatrixTools.toString(contextInsertionPSAM));
                                    out.println((byte)0, "**********************************************************");

                                    motifToHammingSphereMap.put(contextInsertionPSAM, hammingSphere);
                                    motifsToModsMap.put(contextInsertionPSAM, modsArray);
                                    motifsToMandatoriesMap.put(contextInsertionPSAM, mandatoryColumns);

                                }
                            }
                        }
                    }

                    // perform penalized regression to find significant features
                    regressionFitData = fitModel(
                        motifToHammingSphereMap,
                        null, //startPositions
                        positionalWeights,
                        motifsToModsMap,
                        motifsToMandatoriesMap,
                        "context-dependent, nucleotide-dependent insertion features",
                        featureDiscoveryAlgorithm,
                        effectiveStrand,
                        strand,
                        calc,
                        null, //startParams
                        null, //refSymList
                        Double.NaN, //refKaReal);
                        unifiedColumn);

                    fsam.addFeatures(
                        !regressionFitData.isCoefficientZero(hamming1PSAM),
                        regressionFitData.getFeatures(hamming1PSAM, featuresScalar));
                    //fsam.setIntercept(regressionFitData.getIntercept(hamming1PSAM));

                    // coefficients are as returned in R
                    // First coefficient is the intercept, then the rest are in independent-var order
                    out.print(regressionFitData.toString());

                } // end - Perform regression to find context-dependent, inserted nucleotide-dependent insertion features


                ////////////////////////////////////////////////////////////////////////////////
                // Perform regression to find context-independent deletion features
                ////////////////////////////////////////////////////////////////////////////////
                if (((String)initFileCache.get("AFFINITY MODELS", "Get_CI_DeletionAffinities")).equalsIgnoreCase("Yes")) {

                    motifsToModsMap = new LinkedHashMap<Object, Symbol[]>();
                    motifsToMandatoriesMap = new LinkedHashMap<Object, boolean[][]>();
                    motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();

                    if (useHammingDistForApprox) {
                        motifToHammingSphereMap.put(hamming1PSAM, hammingSpheres.get(hammingDistForApprox));
                    }
                    else {
                        motifToHammingSphereMap.put(hamming1PSAM, null);
                    }

                    // Add all possible single-point deletions to the motifToHammingSphereMap
                    // DO NOT delete another column that is in extension
                    // 0-based indexing
                    for (int del=0+extensionLength; del < motifLength-extensionLength; del++) {
                        //for (int del=0; del < motifLength; del+=motifLength-1) {

                        boolean[][] mandatoryColumns = new boolean[2][motifLength-1];
                        mandatoryColumns[0][del] = true;
                        mandatoryColumns[0][del-1] = true;
                        mandatoryColumns[0][del+1] = true;
                        mandatoryColumns[1] = ArrayTools.reverse(mandatoryColumns[0]);

                        // create deletion in the PSAM
                        WeightMatrix singleDeletionPSAM = WeightMatrixTools.delete(hamming1PSAM, del);
                        //singleDeletionPSAM.setName("[D("+ (del+1) +")]");
                        singleDeletionPSAM.setName("delete_"+ (del+1) +"");
                        Symbol[] modsArray = new Symbol[(motifLength*2)-1];
                        modsArray[(del*2)] = DNATools.n();

                        // create the hamming sphere around the mutation
                        LinkedHashSet<SymbolList> hammingSphere = null;

                        if (useHammingDistForApprox) {

                            //hammingSphere = hammingSpheres.get(hammingDistForApprox);

                            // create the deletion in the seedSymList
                            // 1-based indexing
                            SimpleSymbolList editedSeed = new SimpleSymbolList(seedSymList);
                            Edit edit = new Edit(del+1, 1, SymbolList.EMPTY_LIST);
                            editedSeed.edit(edit);

                            hammingSphere = SymbolListTools.getHammingSphere(editedSeed, hammingDistForApprox, includeRevComps, complementTable);
                            out.println((byte)0, "hammingSphere = "+symListCollectionDump(hammingSphere));
                        }

                        out.println((byte)0, "WeightMatrix: "+WeightMatrixTools.toString(singleDeletionPSAM));
                        out.println((byte)0, "**********************************************************");

                        motifToHammingSphereMap.put(singleDeletionPSAM, hammingSphere);
                        motifsToModsMap.put(singleDeletionPSAM, modsArray);
                        motifsToMandatoriesMap.put(singleDeletionPSAM, mandatoryColumns);

                    }

                    // perform penalized regression to find significant features
                    regressionFitData = fitModel(
                        motifToHammingSphereMap,
                        null, //startParams
                        positionalWeights,
                        motifsToModsMap,
                        motifsToMandatoriesMap,
                        "context-independent nucleotide deletion features",
                        featureDiscoveryAlgorithm,
                        effectiveStrand,
                        strand,
                        calc,
                        null, //startParams
                        null, //refSymList
                        Double.NaN, //refKaReal);
                        unifiedColumn);

                    fsam.addFeatures(
                        !regressionFitData.isCoefficientZero(hamming1PSAM),
                        regressionFitData.getFeatures(hamming1PSAM, featuresScalar));
                    //fsam.setIntercept(regressionFitData.getIntercept(hamming1PSAM));

                    // coefficients are as returned in R
                    // First coefficient is the intercept, then the rest are in independent-var order
                    out.print(regressionFitData.toString());

                } // end - Perform regression to find context-independent deletion features


                ////////////////////////////////////////////////////////////////////////////////
                // Perform regression to find context-dependent deletion features
                ////////////////////////////////////////////////////////////////////////////////
                if (((String)initFileCache.get("AFFINITY MODELS", "Get_CD_DeletionAffinities")).equalsIgnoreCase("Yes")) {

                    motifsToModsMap = new LinkedHashMap<Object, Symbol[]>();
                    motifsToMandatoriesMap = new LinkedHashMap<Object, boolean[][]>();
                    motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();

                    if (useHammingDistForApprox) {
                        motifToHammingSphereMap.put(hamming1PSAM, hammingSpheres.get(hammingDistForApprox));
                    }
                    else {
                        motifToHammingSphereMap.put(hamming1PSAM, null);
                    }

                    // Add all possible single-point deletions to the motifToHammingSphereMap
                    // DO delete first and last position
                    // 0-based indexing
                    for (int del=1+extensionLength; del < motifLength-1-extensionLength; del++) {
                        //for (int del=1; del < motifLength-1; del += (motifLength-3)) {

                        //for (int del=0; del < motifLength; del++) {

                        boolean[][] mandatoryColumns = new boolean[2][motifLength-1];
                        mandatoryColumns[0][del] = true;
                        mandatoryColumns[0][del-1] = true;
                        mandatoryColumns[0][del+1] = true;
                        mandatoryColumns[1] = ArrayTools.reverse(mandatoryColumns[0]);

                        // iterate over every symbol at the mut1 position (to the left of del)
                        for (Iterator dnaIter2 = ((FiniteAlphabet)alphabet).iterator(); dnaIter2.hasNext(); ) {
                            Symbol mut1Symbol = (Symbol)dnaIter2.next();

                            // iterate over every symbol at the mut2 position (to the right of del)
                            for (Iterator dnaIter3 = ((FiniteAlphabet)alphabet).iterator(); dnaIter3.hasNext(); ) {
                                Symbol mut2Symbol = (Symbol)dnaIter3.next();

                                // create deletion and the 2 mutations in the PSAM
                                WeightMatrix contextDeletionPSAM = WeightMatrixTools.delete(hamming1PSAM, del);
                                contextDeletionPSAM = WeightMatrixTools.mutate(contextDeletionPSAM, del-1, mut1Symbol, false);
                                contextDeletionPSAM = WeightMatrixTools.mutate(contextDeletionPSAM, del, mut2Symbol, false);

                                // name has 1-based indexing
                                contextDeletionPSAM.setName(
                                    //                                     "["+(del)+"("+mut1Symbol.getName().substring(0,1).toUpperCase()
                                    //                                     +"), D("+ (del+1)
                                    //                                     +"), "+(del+2)+"("+mut2Symbol.getName().substring(0,1).toUpperCase()+")]");
                                    "mutate_"+(del)+"_"+mut1Symbol.getName().substring(0,1).toUpperCase()
                                    +"__"+"delete_"+ (del+1)
                                    +"__"+"mutate_"+(del+2)+"_"+mut2Symbol.getName().substring(0,1).toUpperCase()+"");

                                Symbol[] modsArray = new Symbol[(motifLength*2)-1];
                                modsArray[(del*2)-2] = mut1Symbol;
                                modsArray[(del*2)] = DNATools.n();
                                modsArray[(del*2)+2] = mut2Symbol;

                                LinkedHashSet<SymbolList> hammingSphere = null;

                                if (useHammingDistForApprox) {

                                    //hammingSphere = hammingSpheres.get(hammingDistForApprox);

                                    // create the deletion in the seedSymList
                                    // 1-based indexing
                                    SimpleSymbolList editedSeed = new SimpleSymbolList(seedSymList);
                                    Edit edit = new Edit(del+1, 1, SymbolList.EMPTY_LIST);
                                    editedSeed.edit(edit);
                                    edit = new Edit((del+1)-1, alphabet, mut1Symbol);
                                    editedSeed.edit(edit);
                                    edit = new Edit((del+1), alphabet, mut2Symbol);
                                    editedSeed.edit(edit);

                                    // create the hamming sphere around the mutation
                                    hammingSphere = SymbolListTools.getHammingSphere(editedSeed, hammingDistForApprox, includeRevComps, complementTable);
                                    out.println((byte)0, "hammingSphere = "+symListCollectionDump(hammingSphere));
                                }
                                out.println((byte)0, "WeightMatrix : "+WeightMatrixTools.toString(contextDeletionPSAM));
                                out.println((byte)0, "**********************************************************");

                                motifsToModsMap.put(contextDeletionPSAM, modsArray);
                                motifToHammingSphereMap.put(contextDeletionPSAM, hammingSphere);
                                motifsToMandatoriesMap.put(contextDeletionPSAM, mandatoryColumns);
                            }
                        }
                    }

                    // perform penalized regression to find significant features
                    regressionFitData = fitModel(
                        motifToHammingSphereMap,
                        null, //startParams
                        positionalWeights,
                        motifsToModsMap,
                        motifsToMandatoriesMap,
                        "context-dependent nucleotide deletion features",
                        featureDiscoveryAlgorithm,
                        effectiveStrand,
                        strand,
                        calc,
                        null, //startParams
                        null, //refSymList
                        Double.NaN, //refKaReal);
                        unifiedColumn);

                    fsam.addFeatures(
                        !regressionFitData.isCoefficientZero(hamming1PSAM),
                        regressionFitData.getFeatures(hamming1PSAM, featuresScalar));
                    //fsam.setIntercept(regressionFitData.getIntercept(hamming1PSAM));

                    // coefficients are as returned in R
                    // First coefficient is the intercept, then the rest are in independent-var order
                    out.print(regressionFitData.toString());

                } // end - Perform regression to find context-dependent deletion features



                ////////////////////////////////////////////////////////////////////////////////////////////////
                // Perform regression for current fsam!!
                ////////////////////////////////////////////////////////////////////////////////////////////////
                //                 rSquared = getRsquared(
                //                     fsam,
                //                     positionalWeights,
                //                     "fit with the current FSAM",
                //                     "lm",
                //                     strand,
                //                     strand,
                //                     calc,
                //                     null);



                ////////////////////////////////////////////////////////////////////////////////
                // Perform regression to find Nearest-Neighbor dinucleotide dependencies
                ////////////////////////////////////////////////////////////////////////////////

                double[] intensities = reduceData.getUsedIntensities(1, unifiedColumn);
                double[] residuals = intensities;

                out.println("\nThe normalized L1 Distance between the smoothed and raw positional bias weights is "+posBiasWeightsL1);

                double posBiasWeightsL1Threshold = initFileCache.getDouble("AFFINITY MODELS", "PosBiasWeightsL1Threshold");

                if (posBiasWeightsL1 > posBiasWeightsL1Threshold) {
                    out.println("\nThe normalized L1 Distance "+posBiasWeightsL1+" between the smoothed and raw positional bias weights is greater than the threshold "+posBiasWeightsL1Threshold+", so dinucleotide dependencies will NOT be attempted.");
                }
                else {

                    if (((String)initFileCache.get("AFFINITY MODELS", "Get_NNDD_AdditiveAffinities")).equalsIgnoreCase("Yes")) {

                        WeightMatrix hamming1RevComp = WeightMatrixTools.reverseComplement(hamming1PSAM, complementTable);

                        motifsToModsMap = new LinkedHashMap<Object, Symbol[]>();
                        motifsToMandatoriesMap = new LinkedHashMap<Object, boolean[][]>();
                        motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();

                        if (!featureDiscoveryAlgorithm.equalsIgnoreCase("lars")) {
                            if (useHammingDistForApprox) {
                                motifToHammingSphereMap.put(hamming1PSAM, hammingSpheres.get(hammingDistForApprox));
                            }
                            else {
                                motifToHammingSphereMap.put(hamming1PSAM, null);
                            }
                        }

                        //                     java.util.List avoidDimerResults = WeightMatrixTools.getMostLikelySymList(hamming1PSAM, 2);
                        //                     SymbolList avoidDimer =(SymbolList) avoidDimerResults.get(0);
                        //                     int avoidDimerStartPos = ((Integer)avoidDimerResults.get(1)).intValue();
                        //                     out.println("\nMost Likely Dimer is "+avoidDimer.seqString().toUpperCase()+" at start position "+(avoidDimerStartPos+1)+".");

                        // Add all double point mutations to the motifToHammingSphereMap
                        // 0-based indexing
                        for (int mut1=0; mut1 < motifLength-1; mut1++) { // from 0 to (L-1)

                            //for (int mut1=0; mut1 < 2; mut1++) {
                            //for (int mut1=0; mut1 < motifLength-1; mut1+=motifLength-2) {
                            //for (int mut1=0; mut1 < 1; mut1++) {

                            for (int mut2 = mut1 + 1; mut2 <= mut1 + 1; mut2++) { // Nearest-Neighbor
                                //for (int mut2 = mut1 + 1; mut2 < motifLength; mut2++) { // All

                                boolean[][] mandatoryColumns = new boolean[2][motifLength];
                                mandatoryColumns[0][mut1] = true;
                                mandatoryColumns[0][mut2] = true;
                                mandatoryColumns[1] = ArrayTools.reverse(mandatoryColumns[0]);

                                // iterate over every symbol at the mut1 position
                                for (Iterator dnaIter2 = ((FiniteAlphabet)alphabet).iterator(); dnaIter2.hasNext(); ) {
                                    Symbol mut1Symbol = (Symbol)dnaIter2.next();

                                    // iterate over every symbol at the mut2 position
                                    for (Iterator dnaIter3 = ((FiniteAlphabet)alphabet).iterator(); dnaIter3.hasNext(); ) {
                                        Symbol mut2Symbol = (Symbol)dnaIter3.next();

                                        // if the double-point mutation is different from the most likely nucleotides in the PSAM
                                        // then add them to the dinucleotide dependencies regression model

                                        /////////////////////////////////////////////////
                                        Symbol avoidMut1Symbol = null;
                                        Symbol avoidMut2Symbol = null;

                                        // if (isRevCompPalindrome) {
                                        //     avoidMut1Symbol = DistributionTools.getLeastLikelySymbol(hamming1PSAM.getColumn(mut1));
                                        //     avoidMut2Symbol = DistributionTools.getLeastLikelySymbol(hamming1PSAM.getColumn(mut2));
                                        // }
                                        // else {
                                        //     avoidMut1Symbol = DistributionTools.getMostLikelySymbol(hamming1RevComp.getColumn(mut1));
                                        //     avoidMut2Symbol = DistributionTools.getMostLikelySymbol(hamming1RevComp.getColumn(mut2));
                                        // }

                                        //Symbol avoidMut1Symbol = DistributionTools.getMostLikelySymbol(hamming1PSAM.getColumn(mut1));
                                        //Symbol avoidMut2Symbol = DistributionTools.getMostLikelySymbol(hamming1PSAM.getColumn(mut2));

                                        avoidMut1Symbol = DistributionTools.getLeastLikelySymbol(hamming1PSAM.getColumn(mut1));
                                        avoidMut2Symbol = DistributionTools.getLeastLikelySymbol(hamming1PSAM.getColumn(mut2));
                                        /////////////////////////////////////////////////

                                        //out.println((byte)100, "avoidMut1Symbol="+avoidMut1Symbol.getName().substring(0,1).toUpperCase()+"\tavoidMut2Symbol="+avoidMut2Symbol.getName().substring(0,1).toUpperCase());

                                        //                                     double avoidMut1Prob = hamming1PSAM.getColumn(mut1).getWeight(avoidMut1Symbol);
                                        //                                     double avoidMut2Prob = hamming1PSAM.getColumn(mut2).getWeight(avoidMut2Symbol);
                                        //                                     double dualProb = avoidMut1Prob * avoidMut2Prob;

                                        //                                     if (   ((mut1 == 2) && (mut1Symbol == DNATools.c()) && (mut2Symbol == DNATools.a()))
                                        //                                         || ((mut1 == 3) && (mut1Symbol == DNATools.a()) && (mut2Symbol == DNATools.c()))
                                        //                                         || ((mut1 == 4) && (mut1Symbol == DNATools.c()) && (mut2Symbol == DNATools.g()))
                                        //                                         || ((mut1 == 5) && (mut1Symbol == DNATools.g()) && (mut2Symbol == DNATools.t()))
                                        //                                         || ((mut1 == 6) && (mut1Symbol == DNATools.t()) && (mut2Symbol == DNATools.g()))
                                        //                                            ) {

                                        //if (false) {
                                        if ( (mut1Symbol == avoidMut1Symbol) && (mut2Symbol == avoidMut2Symbol) ) {

                                            //if ( (mut1 == avoidDimerStartPos) && (mut1Symbol == avoidDimer.symbolAt(1)) && (mut2Symbol == avoidDimer.symbolAt(2)) ) {
                                            //if ( (mut1Symbol == avoidMut1Symbol) && (mut2Symbol == avoidMut2Symbol) && (dualProb > .2) ) {

                                            /////////////////////////////////////////////////
                                            //out.println((byte)100, "\nNot including high affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression.");

                                            out.println((byte)100, "\nNot including low affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression.");

                                            // if (isRevCompPalindrome) {
                                            //     out.println((byte)100, "\nNot including low affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression.");
                                            // }
                                            // else {
                                            //     out.println((byte)100, "\nNot including the high affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) found in the reverse complement in regression.");
                                            // }
                                            /////////////////////////////////////////////////

                                            //out.println((byte)100, "\nNot including high affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression (dualProb = "+dualProb+", threshold = "+0.2+").");
                                        }
                                        //                                     else if ( (iterNum > 1) && (mut1Symbol == avoidMut1bSymbol) && (mut2Symbol == avoidMut2bSymbol) ) {

                                        //                                         //if ( (mut1 == avoidDimerStartPos) && (mut1Symbol == avoidDimer.symbolAt(1)) && (mut2Symbol == avoidDimer.symbolAt(2)) ) {
                                        //                                         //if ( (mut1Symbol == avoidMut1Symbol) && (mut2Symbol == avoidMut2Symbol) && (dualProb > .2) ) {

                                        //                                         //out.println((byte)100, "\nNot including high affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression.");
                                        //                                         out.println((byte)100, "\nNot including low affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression.");

                                        //                                         //out.println((byte)100, "\nNot including high affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression (dualProb = "+dualProb+", threshold = "+0.2+").");
                                        //                                     }
                                        else {

                                            WeightMatrix doubleMutationPSAM = WeightMatrixTools.mutate(hamming1PSAM, mut1, mut1Symbol, true);
                                            doubleMutationPSAM = WeightMatrixTools.mutate(doubleMutationPSAM, mut2, mut2Symbol, false);

                                            // name has 1-based indexing
                                            doubleMutationPSAM.setName(
                                                // "["+ (mut1+1) +"("+ mut1Symbol.getName().substring(0,1).toUpperCase()
                                                // + "), "+ (mut2+1) +"("+ mut2Symbol.getName().substring(0,1).toUpperCase() +")]");
                                                "mutate_"+ (mut1+1) +"_"+ mut1Symbol.getName().substring(0,1).toUpperCase()
                                                +"__"+ "mutate_"+ (mut2+1) +"_"+ mut2Symbol.getName().substring(0,1).toUpperCase() +"");

                                            // modsArray includes 1bp insertions
                                            Symbol[] modsArray = new Symbol[(motifLength*2)-1];
                                            modsArray[(mut1*2)] = mut1Symbol;
                                            modsArray[(mut2*2)] = mut2Symbol;

                                            LinkedHashSet<SymbolList> hammingSphere = null;

                                            if (useHammingDistForApprox) {

                                                hammingSphere = hammingSpheres.get(hammingDistForApprox);

                                                //                                         // create the 2 point-mutations in the seed
                                                //                                         // 1-based indexing
                                                //                                         SimpleSymbolList editedSeed = new SimpleSymbolList(seedSymList);
                                                //                                         Edit edit = new Edit(mut1+1, alphabet, mut1Symbol);
                                                //                                         editedSeed.edit(edit);
                                                //                                         edit = new Edit(mut2+1, alphabet, mut2Symbol);
                                                //                                         editedSeed.edit(edit);

                                                //                                         // create the hamming sphere around the double point-mutation
                                                //                                         hammingSphere = SymbolListTools.getHammingSphere(editedSeed, 1, includeRevComps, complementTable);
                                                //                                         out.println((byte)0, "hammingSphere = "+symListCollectionDump(hammingSphere));
                                            }

                                            //                                     out.println((byte)0, "WeightMatrix: "+WeightMatrixTools.toString(doubleMutationPSAM));
                                            //                                     out.println((byte)0, "**********************************************************");

                                            motifsToModsMap.put(doubleMutationPSAM, modsArray);
                                            motifToHammingSphereMap.put(doubleMutationPSAM, hammingSphere);
                                            motifsToMandatoriesMap.put(doubleMutationPSAM, mandatoryColumns);
                                        }
                                    }
                                }
                            }
                        }

                        //////////////////////////////////////////
                        // Do LARS regression
                        //////////////////////////////////////////
                        if (featureDiscoveryAlgorithm.equalsIgnoreCase("lars")) {

                            // perform penalized regression to find significant features
                            regressionFitData = fitModel(
                                motifToHammingSphereMap,
                                null, //startPositions
                                positionalWeights,
                                motifsToModsMap,
                                motifsToMandatoriesMap,
                                "nearest-neighbor additive dinucleotide dependencies",
                                featureDiscoveryAlgorithm,
                                //"lm",
                                effectiveStrand,
                                strand,
                                //WeightMatrixTools.BindingStrand.BOTH,
                                calc,
                                null, //startParams
                                null, //refSymList
                                Double.NaN, //refKaReal);
                                unifiedColumn);

                            //regressionFitData = fitModel(motifToHammingSphereMap, null, positionalWeights, motifsToModsMap, "nearest-neighbor dinucleotide dependencies", "lm", strand);

                            // coefficients are as returned in R
                            // First coefficient is the intercept, then the rest are in independent-var order
                            out.print(regressionFitData.toString());

                            //yeast
                            motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 5.0);
                            //motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 3.0);

                            //human
                            //motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 2.0);
                            //motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 1.8);

                            regressionFitData = fitModel(
                                motifToHammingSphereMap,
                                null, //startPositions
                                positionalWeights,
                                motifsToModsMap,
                                motifsToMandatoriesMap,
                                "nearest-neighbor additive dinucleotide dependencies",
                                //featureDiscoveryAlgorithm,
                                //"glm.poisson",
                                //"rlm",
                                "lm",
                                effectiveStrand,
                                strand,
                                //WeightMatrixTools.BindingStrand.BOTH,
                                calc,
                                null, //startParams
                                null, //refSymList
                                Double.NaN, //refKaReal);
                                unifiedColumn);

                        }

                        ////////////////////////////////////////////////////////////////////////////////////
                        // Do "featureDiscoveryAlgorithm" regression with psamOffset
                        ////////////////////////////////////////////////////////////////////////////////////
                        else {

                            //numRegressions = 10;

                            // Perform univariate regression for current fsam!!
                            regressionFitData = getRegressionFit(
                                //hamming1PSAM,
                                fsam,
                                positionalWeights,
                                "fit with the current FSAM",
                                //"lm",
                                featureDiscoveryAlgorithm,
                                effectiveStrand,
                                strand,
                                calc,
                                unifiedColumn);

                            this.fsamUnivarCoeff = regressionFitData.getCoefficient(fsam);
                            this.fsamUnivarInter = regressionFitData.getIntercept();

                            psamOffset = true;

                            regressionFitData = fitModel(
                                motifToHammingSphereMap,
                                null, //startPositions
                                positionalWeights,
                                motifsToModsMap,
                                motifsToMandatoriesMap,
                                "nearest-neighbor additive dinucleotide dependencies",
                                featureDiscoveryAlgorithm,
                                //"lm",
                                effectiveStrand,
                                strand,
                                //WeightMatrixTools.BindingStrand.BOTH,
                                calc,
                                null, //startParams
                                null, //refSymList
                                Double.NaN, //refKaReal);
                                unifiedColumn);

                            psamOffset = false;


                            numRegressions = 1;


                            //////////////////////////////////////////
                            // If necessary, Do "rlm.default" regression
                            //////////////////////////////////////////
                            if (regressionFitData == null) {
                                psamOffset = true;

                                regressionFitData = fitModel(
                                    motifToHammingSphereMap,
                                    null, //startPositions
                                    positionalWeights,
                                    motifsToModsMap,
                                    motifsToMandatoriesMap,
                                    "nearest-neighbor additive dinucleotide dependencies",
                                    //featureDiscoveryAlgorithm,
                                    "rlm.default",
                                    //"lm",
                                    effectiveStrand,
                                    strand,
                                    //WeightMatrixTools.BindingStrand.BOTH,
                                    calc,
                                    null, //startParams
                                    null, //refSymList
                                    Double.NaN, //refKaReal);
                                    unifiedColumn);

                                psamOffset = false;
                            }

                        }

                        if (regressionFitData == null) {

                            out.println("WARNING - RegressionFitData is null! Could not fit dinucleotide dependency terms.");

                        }
                        else { // add the dinuc corrections to the model

                            out.print(regressionFitData.toString());

                            fsam.addFeatures(
                                !regressionFitData.isCoefficientZero(hamming1PSAM),
                                regressionFitData.getFeatures(hamming1PSAM, featuresScalar, true),
                                true); // these are additive features

                            //fsam.setIntercept(regressionFitData.getIntercept(hamming1PSAM));
                            out.print("\n"+fsam.toString());
				//add boolean additive true in the above statement - pa2399
                            regressionFitData = getRegressionFit(
                                fsam,
                                positionalWeights,
                                "fit with FSAM before equilibrate",
                                "lm",
                                effectiveStrand,
                                strand,
                                calc,
                                unifiedColumn);
                            System.out.println("R2 with FSAM before equilibrate is "+regressionFitData.rSquared);

                            /////////////////////////////////////////////////////////////////////////////////////
                            // equilibrate - creates a new PSAM and new multiplicative features
                            /////////////////////////////////////////////////////////////////////////////////////
                            fsam.equilibrate();

                            regressionFitData = getRegressionFit(
                                fsam,
                                positionalWeights,
                                "fit with FSAM after equilibrate",
                                "lm",
                                effectiveStrand,
                                strand,
                                calc,
                                unifiedColumn);
                            System.out.println("R2 with FSAM after equilibrate is "+regressionFitData.rSquared);

                            // coefficients are as returned in R
                            // First coefficient is the intercept, then the rest are in independent-var order
                            out.print(regressionFitData.toString());

                            out.print("\n"+fsam.toString());

                            /////////////////////////////////////////////////////////////////////////////////////
                            // set nonSpecKa
                            /////////////////////////////////////////////////////////////////////////////////////
                            nonSpecKa.setValue(getRelativeKaNonSpec(
                                    regressionFitData.getIntercept(),
                                    regressionFitData.getCoefficient(fsam),
                                    positionalWeights,
                                    eToMu,
                                    probeSeqLengths,
                                    motifLength,
                                    strand,
                                    calc,
                                    fsam.getRevCompSimilarity()));

                            fsam.setNonSpecKa(nonSpecKa.doubleValue());

                            regressionFitData = getRegressionFit(
                                fsam,
                                positionalWeights,
                                "fit with FSAM after setting nonSpecKa",
                                "lm",
                                effectiveStrand,
                                strand,
                                calc,
                                unifiedColumn);
                            System.out.println("R2 with FSAM after setting nonSpecKa is "+regressionFitData.rSquared);


                            // set the scaler and intercept for the FSAM
                            //this.intercept = regressionFitData.getIntercept();
                            //this.intensitiesScaler = regressionFitData.getCoefficient(fsam);
                            fsam.setIntercept(regressionFitData.getIntercept());
                            fsam.setScaler(regressionFitData.getCoefficient(fsam));


                            //residuals = MathTools.subtract(residuals, regressionFitData.getResiduals());

                            //residuals = regressionFitData.getResiduals();

                        } // end - if (regressionFitData != null)

                    } // end - Perform regression to find Additive Nearest-Neighbor dinucleotide dependencies

                    ////////////////////////////////////////////////////////
                    // Set hamming1PSAM after re-equilibration
                    hamming1PSAM = fsam.getPosStrandPWM();
                    ////////////////////////////////////////////////////////

                    ////////////////////////////////////////////////////////////////////////////////
                    // Perform regression to find Nearest-Neighbor dinucleotide dependencies
                    ////////////////////////////////////////////////////////////////////////////////
                    if (((String)initFileCache.get("AFFINITY MODELS", "Get_NN_DinucleotideAffinities")).equalsIgnoreCase("Yes")) {

                        // get all the symbol lists, and the hamming Sphere around the highest affinity seq
                        System.out.print("\nGenerating all sequences of length "+hamming1PSAM.columns()+"...");
                        LinkedHashSet<SymbolList> allSymLists = SymbolListTools.getAllSymbolLists(this.alphabet, hamming1PSAM.columns());
                        System.out.println("Done.");

                        int numNNDinucleotideInterations = (int)initFileCache.getInt("REGRESSION", "NumNNDinucleotideInterations");

                        // int numLoops = 10;
                        // int numLoops = 1;
                        // for (int loopCounter=0; loopCounter < numLoops; loopCounter++) {
                        for (int loopCounter=0; loopCounter < numNNDinucleotideInterations; loopCounter++) {
                            // Add all double point mutations to the motifToHammingSphereMap
                            // 0-based indexing
                            for (int mut1=0; mut1 < motifLength-1; mut1++) { // from 0 to (L-1)
                                //for (int mut1=0; mut1 < 5; mut1++) { // from 0 to (L-1)
                                //for (int mut1=0; mut1 < motifLength-1; mut1 +=2) { // from 0 to (L-1)

                                //for (int mut1=0; mut1 < 2; mut1++) {
                                //for (int mut1=0; mut1 < motifLength-1; mut1+=motifLength-2) {
                                //for (int mut1=0; mut1 < 1; mut1++) {

                                for (int mut2 = mut1 + 1; mut2 <= mut1 + 1; mut2++) { // Nearest-Neighbor
                                    //for (int mut2 = mut1 + 1; mut2 < motifLength; mut2++) { // All

                                    motifsToModsMap = new LinkedHashMap<Object, Symbol[]>();
                                    motifsToMandatoriesMap = new LinkedHashMap<Object, boolean[][]>();
                                    motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();
                                    aFeatureRemovalList = new ArrayList<FeatureKey>();

                                    boolean[][] mandatoryColumns = new boolean[2][motifLength];
                                    mandatoryColumns[0][mut1] = true;
                                    mandatoryColumns[0][mut2] = true;
                                    mandatoryColumns[1] = ArrayTools.reverse(mandatoryColumns[0]);

                                    // create an fsam clone that removes all the features at contain only mut1 and mut2 positions
                                    // FeaturedWeightMatrix fsamMinusMuts = fsam.clone(0.0, mandatoryColumns[0]);

                                    if (!featureDiscoveryAlgorithm.equalsIgnoreCase("lars")) {
                                        if (useHammingDistForApprox) {
                                            motifToHammingSphereMap.put(hamming1PSAM, hammingSpheres.get(hammingDistForApprox));
                                        }
                                        else {
                                            motifToHammingSphereMap.put(hamming1PSAM, null);
                                        }
                                    }
                                    //                     java.util.List avoidDimerResults = WeightMatrixTools.getMostLikelySymList(hamming1PSAM, 2);
                                    //                     SymbolList avoidDimer =(SymbolList) avoidDimerResults.get(0);
                                    //                     int avoidDimerStartPos = ((Integer)avoidDimerResults.get(1)).intValue();
                                    //                     out.println("\nMost Likely Dimer is "+avoidDimer.seqString().toUpperCase()+" at start position "+(avoidDimerStartPos+1)+".");


                                    // iterate over every symbol at the mut1 position
                                    for (Iterator dnaIter2 = ((FiniteAlphabet)alphabet).iterator(); dnaIter2.hasNext(); ) {
                                        Symbol mut1Symbol = (Symbol)dnaIter2.next();

                                        // iterate over every symbol at the mut2 position
                                        for (Iterator dnaIter3 = ((FiniteAlphabet)alphabet).iterator(); dnaIter3.hasNext(); ) {
                                            Symbol mut2Symbol = (Symbol)dnaIter3.next();

                                            // if the double-point mutation is different from the most likely nucleotides in the PSAM
                                            // then add them to the dinucleotide dependencies regression model

                                            Symbol avoidMut1Symbol = DistributionTools.getMostLikelySymbol(hamming1PSAM.getColumn(mut1));
                                            Symbol avoidMut2Symbol = DistributionTools.getMostLikelySymbol(hamming1PSAM.getColumn(mut2));


                                            //Symbol avoidMut1Symbol = DistributionTools.getLeastLikelySymbol(hamming1PSAM.getColumn(mut1));
                                            //Symbol avoidMut2Symbol = DistributionTools.getLeastLikelySymbol(hamming1PSAM.getColumn(mut2));
                                            //out.println((byte)100, "avoidMut1Symbol="+avoidMut1Symbol.getName().substring(0,1).toUpperCase()+"\tavoidMut2Symbol="+avoidMut2Symbol.getName().substring(0,1).toUpperCase());

                                            //                                     double avoidMut1Prob = hamming1PSAM.getColumn(mut1).getWeight(avoidMut1Symbol);
                                            //                                     double avoidMut2Prob = hamming1PSAM.getColumn(mut2).getWeight(avoidMut2Symbol);
                                            //                                     double dualProb = avoidMut1Prob * avoidMut2Prob;

                                            //                                     if (   ((mut1 == 2) && (mut1Symbol == DNATools.c()) && (mut2Symbol == DNATools.a()))
                                            //                                         || ((mut1 == 3) && (mut1Symbol == DNATools.a()) && (mut2Symbol == DNATools.c()))
                                            //                                         || ((mut1 == 4) && (mut1Symbol == DNATools.c()) && (mut2Symbol == DNATools.g()))
                                            //                                         || ((mut1 == 5) && (mut1Symbol == DNATools.g()) && (mut2Symbol == DNATools.t()))
                                            //                                         || ((mut1 == 6) && (mut1Symbol == DNATools.t()) && (mut2Symbol == DNATools.g()))
                                            //                                            ) {

                                            if (false) {
                                                //if ( (mut1Symbol == avoidMut1Symbol) && (mut2Symbol == avoidMut2Symbol) ) {

                                                //if ( (mut1 == avoidDimerStartPos) && (mut1Symbol == avoidDimer.symbolAt(1)) && (mut2Symbol == avoidDimer.symbolAt(2)) ) {
                                                //if ( (mut1Symbol == avoidMut1Symbol) && (mut2Symbol == avoidMut2Symbol) && (dualProb > .2) ) {

                                                out.println((byte)100, "\nNot including high affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression.");
                                                //out.println((byte)100, "\nNot including low affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression.");
                                                //out.println((byte)100, "\nNot including high affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression (dualProb = "+dualProb+", threshold = "+0.2+").");
                                            }
                                            else {

                                                WeightMatrix doubleMutationPSAM = WeightMatrixTools.mutate(hamming1PSAM, mut1, mut1Symbol, true);
                                                doubleMutationPSAM = WeightMatrixTools.mutate(doubleMutationPSAM, mut2, mut2Symbol, false);

                                                // name has 1-based indexing
                                                doubleMutationPSAM.setName(
                                                    // "["+ (mut1+1) +"("+ mut1Symbol.getName().substring(0,1).toUpperCase()
                                                    // + "), "+ (mut2+1) +"("+ mut2Symbol.getName().substring(0,1).toUpperCase() +")]");
                                                    "mutate_"+ (mut1+1) +"_"+ mut1Symbol.getName().substring(0,1).toUpperCase()
                                                    +"__"+ "mutate_"+ (mut2+1) +"_"+ mut2Symbol.getName().substring(0,1).toUpperCase() +"");

                                                Symbol[] modsArray = new Symbol[(motifLength*2)-1];
                                                modsArray[(mut1*2)] = mut1Symbol;
                                                modsArray[(mut2*2)] = mut2Symbol;

                                                LinkedHashSet<SymbolList> hammingSphere = null;

                                                if (useHammingDistForApprox) {

                                                    hammingSphere = hammingSpheres.get(hammingDistForApprox);

                                                    //                                         // create the 2 point-mutations in the seed
                                                    //                                         // 1-based indexing
                                                    //                                         SimpleSymbolList editedSeed = new SimpleSymbolList(seedSymList);
                                                    //                                         Edit edit = new Edit(mut1+1, alphabet, mut1Symbol);
                                                    //                                         editedSeed.edit(edit);
                                                    //                                         edit = new Edit(mut2+1, alphabet, mut2Symbol);
                                                    //                                         editedSeed.edit(edit);

                                                    //                                         // create the hamming sphere around the double point-mutation
                                                    //                                         hammingSphere = SymbolListTools.getHammingSphere(editedSeed, 1, includeRevComps, complementTable);
                                                    //                                         out.println((byte)0, "hammingSphere = "+symListCollectionDump(hammingSphere));
                                                }

                                                //                                     out.println((byte)0, "WeightMatrix: "+WeightMatrixTools.toString(doubleMutationPSAM));
                                                //                                     out.println((byte)0, "**********************************************************");

                                                motifsToModsMap.put(doubleMutationPSAM, modsArray);
                                                motifToHammingSphereMap.put(doubleMutationPSAM, hammingSphere);
                                                motifsToMandatoriesMap.put(doubleMutationPSAM, mandatoryColumns);
                                                aFeatureRemovalList.add(new FeatureKey(modsArray));
                                            }
                                        }
                                    }

                                    if (featureDiscoveryAlgorithm.equalsIgnoreCase("lars")) {

                                        // perform penalized regression to find significant features
                                        regressionFitData = fitModel(
                                            motifToHammingSphereMap,
                                            null, //startPositions
                                            positionalWeights,
                                            motifsToModsMap,
                                            motifsToMandatoriesMap,
                                            "nearest-neighbor dinucleotide dependencies",
                                            featureDiscoveryAlgorithm,
                                            //"lm",
                                            effectiveStrand,
                                            strand,
                                            //WeightMatrixTools.BindingStrand.BOTH,
                                            calc,
                                            null, //startParams
                                            null, //refSymList
                                            Double.NaN, //refKaReal);
                                            unifiedColumn);

                                        //regressionFitData = fitModel(motifToHammingSphereMap, null, positionalWeights, motifsToModsMap, "nearest-neighbor dinucleotide dependencies", "lm", strand);

                                        // coefficients are as returned in R
                                        // First coefficient is the intercept, then the rest are in independent-var order
                                        out.print(regressionFitData.toString());

                                        //yeast
                                        motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 5.0);
                                        //motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 3.0);

                                        //human
                                        //motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 2.0);
                                        //motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 1.8);

                                        regressionFitData = fitModel(
                                            motifToHammingSphereMap,
                                            null, //startPositions
                                            positionalWeights,
                                            motifsToModsMap,
                                            motifsToMandatoriesMap,
                                            "nearest-neighbor dinucleotide dependencies",
                                            //featureDiscoveryAlgorithm,
                                            //"glm.poisson",
                                            //"rlm",
                                            "lm",
                                            effectiveStrand,
                                            strand,
                                            //WeightMatrixTools.BindingStrand.BOTH,
                                            calc,
                                            null, //startParams
                                            null, //refSymList
                                            Double.NaN, //refKaReal);
                                            unifiedColumn);

                                    }

                                    ////////////////////////////////////////////////////////////////////////////////////
                                    // Do "featureDiscoveryAlgorithm" regression with psamOffset
                                    ////////////////////////////////////////////////////////////////////////////////////
                                    else {

                                        //numRegressions = 10;

                                        // Perform univariate regression for current fsam!!
                                        regressionFitData = getRegressionFit(
                                            fsam,
                                            positionalWeights,
                                            "fit with the current FSAM",
                                            //"lm",
                                            featureDiscoveryAlgorithm,
                                            effectiveStrand,
                                            strand,
                                            calc,
                                            unifiedColumn);

                                        this.fsamUnivarCoeff = regressionFitData.getCoefficient(fsam);
                                        this.fsamUnivarInter = regressionFitData.getIntercept();

                                        motifToFwmMap = fsam.clones(
                                            motifToHammingSphereMap,
                                            motifsToMandatoriesMap,
                                            aFeatureRemovalList);

                                        psamOffset = true;

                                        regressionFitData = fitModel(
                                            motifToHammingSphereMap,
                                            null, //startPositions
                                            positionalWeights,
                                            motifsToModsMap,
                                            motifsToMandatoriesMap,
                                            "nearest-neighbor additive dinucleotide dependencies",
                                            featureDiscoveryAlgorithm,
                                            //"lm",
                                            effectiveStrand,
                                            strand,
                                            //WeightMatrixTools.BindingStrand.BOTH,
                                            calc,
                                            null, //startParams
                                            null, //refSymList
                                            Double.NaN, //refKaReal);
                                            unifiedColumn);

                                        psamOffset = false;

                                        motifToFwmMap = null;

                                        numRegressions = 1;


                                        //////////////////////////////////////////
                                        // If necessary, Do "rlm.default" regression
                                        //////////////////////////////////////////
                                        if (regressionFitData == null) {

                                            motifToFwmMap = fsam.clones(
                                                motifToHammingSphereMap,
                                                motifsToMandatoriesMap,
                                                aFeatureRemovalList);

                                            psamOffset = true;

                                            regressionFitData = fitModel(
                                                motifToHammingSphereMap,
                                                null, //startPositions
                                                positionalWeights,
                                                motifsToModsMap,
                                                motifsToMandatoriesMap,
                                                "nearest-neighbor additive dinucleotide dependencies",
                                                //featureDiscoveryAlgorithm,
                                                "rlm.default",
                                                //"lm",
                                                effectiveStrand,
                                                strand,
                                                //WeightMatrixTools.BindingStrand.BOTH,
                                                calc,
                                                null, //startParams
                                                null, //refSymList
                                                Double.NaN, //refKaReal);
                                                unifiedColumn);

                                            psamOffset = false;

                                            motifToFwmMap = null;

                                        }

                                    }

                                    if (regressionFitData == null) {

                                        out.println("WARNING - RegressionFitData is null! Could not fit dinucleotide dependency terms.");

                                    }
                                    else { // add the dinuc corrections to the model

                                        out.print(regressionFitData.toString());

                                        // get the new features from this regression
                                        ArrayList<WeightMatrixFeature> newFeatures = regressionFitData.getFeatures(
                                            hamming1PSAM,
                                            featuresScalar,
                                            false); // get multiplicative features

                                        fsam.avgFeatures(
                                            !regressionFitData.isCoefficientZero(hamming1PSAM),
                                            newFeatures,
                                            false, // these are multiplicative features
                                            .3); // newFeature weight

                                        //fsam.setIntercept(regressionFitData.getIntercept(hamming1PSAM));
                                        out.print("\n"+fsam.toString(true));

                                        regressionFitData = getRegressionFit(
                                            fsam,
                                            positionalWeights,
                                            //"fit with FSAM before equilibrate",
                                            "fit with FSAM after setting dinuc features",
                                            "lm",
                                            effectiveStrand,
                                            strand,
                                            calc,
                                            unifiedColumn);

                                        System.out.println("R2 with FSAM after setting dinuc features is "+regressionFitData.rSquared);
                                        //System.out.println("R2 with FSAM before equilibrate is "+regressionFitData.rSquared);

                                        /////////////////////////////////////////////////////////////////////////////////////
                                        // equilibrate - creates a new PSAM and new multiplicative features
                                        /////////////////////////////////////////////////////////////////////////////////////
                                        // fsam.equilibrate();

                                        // regressionFitData = getRegressionFit(
                                        //     fsam,
                                        //     positionalWeights,
                                        //     "fit with FSAM after equilibrate",
                                        //     "lm",
                                        //     effectiveStrand,
                                        //     strand,
                                        //     calc,
                                        //     unifiedColumn);
                                        // System.out.println("R2 with FSAM after equilibrate is "+regressionFitData.rSquared);

                                        // // coefficients are as returned in R
                                        // // First coefficient is the intercept, then the rest are in independent-var order
                                        // out.print(regressionFitData.toString());

                                        // out.print("\n"+fsam.toString());

                                        /////////////////////////////////////////////////////////////////////////////////////
                                        // set nonSpecKa
                                        /////////////////////////////////////////////////////////////////////////////////////
                                        nonSpecKa.setValue(getRelativeKaNonSpec(
                                                regressionFitData.getIntercept(),
                                                regressionFitData.getCoefficient(fsam),
                                                positionalWeights,
                                                eToMu,
                                                probeSeqLengths,
                                                motifLength,
                                                strand,
                                                calc,
                                                fsam.getRevCompSimilarity()));

                                        fsam.setNonSpecKa(nonSpecKa.doubleValue());

                                        regressionFitData = getRegressionFit(
                                            fsam,
                                            positionalWeights,
                                            "fit with FSAM after setting nonSpecKa",
                                            "lm",
                                            effectiveStrand,
                                            strand,
                                            calc,
                                            unifiedColumn);
                                        System.out.println("R2 with FSAM after setting nonSpecKa is "+regressionFitData.rSquared);


                                        // set the scaler and intercept for the FSAM
                                        //this.intercept = regressionFitData.getIntercept();
                                        //this.intensitiesScaler = regressionFitData.getCoefficient(fsam);
                                        fsam.setIntercept(regressionFitData.getIntercept());
                                        fsam.setScaler(regressionFitData.getCoefficient(fsam));


                                        //residuals = MathTools.subtract(residuals, regressionFitData.getResiduals());

                                        //residuals = regressionFitData.getResiduals();

                                    } // end - if (regressionFitData != null)

                                } // mut2 = postions 1 to length-1

                            } // mut1 = positions 0 to length-2

                            out.print("\n"+fsam.toString());

                            //if (loopCounter < (numLoops-1)) {
                            if (loopCounter < 4) {
                                //fsam.equilibrateMultModel(allSymLists);

                                out.print("\n"+fsam.toString());
                            }

                        } // loopCounter = 0 to N-1

                    } // end - Perform regression to find Nearest-Neighbor dinucleotide dependencies


                } // (posBiasWeightsL1 > posBiasWeightsL1Threshold)


                ////////////////////////////////////////////////////////////////////////////////
                // Create the All-Kmer Model
                ////////////////////////////////////////////////////////////////////////////////
                if (initFileCache.getBoolean("AFFINITY MODELS", "GetAllKmerModel")) {

                    // kmer length must <= 8!!!!!
                    //char kmerLength = (char)initFileCache.getInt("AFFINITY MODELS", "KmerModelLength");
                    boolean[] nonNegativeFlags = initFileCache.getBooleans("AFFINITY MODELS", "NonNegativeRegression");
                    boolean[] includeRevCompsArray = initFileCache.getBooleans("AFFINITY MODELS", "IncludeRevComps");

                    //double[] fractionDrops = {.20, .20, .20, .20, .20, .20, .20};

                    double[] fractionDrops = initFileCache.getDoubles("AFFINITY MODELS", "KmerModelFractionDrops");
                    //double[] fractionDrops = {.30, .30, .30, .30, .30, .30, .30};

                    //boolean[] nonNegativeFlags = {true, false, false, false, false};

                    double[] percentSteps = initFileCache.getDoubles("AFFINITY MODELS", "KmerModelPercentSteps");
                    //double[] percentSteps = {.10, .10, .10, .10, .10, .10, .10};

                    //double[] percentSteps = {.15, .15, .15, .15, .15, .15, .15};

                    int[] itersArray = initFileCache.getInts("AFFINITY MODELS", "KmerModelIterations");
                    //int[] itersArray = {100, 100, 100, 100, 100, 100, 100};

                    //String[] initValFlags = {"avg", "0"};
                    String[] initValFlags = {"avg", "avg", "avg", "avg", "avg", "avg", "avg", "avg"};

                    double[] initTopDrops = initFileCache.getDoubles("AFFINITY MODELS", "KmerModelInitTopDrops");
                    double[] initBottomDrops = initFileCache.getDoubles("AFFINITY MODELS", "KmerModelInitBottomDrops");
                    //double[] initTopDrops = {.15, .15, .15, .15, .15, .15, .15};
                    //double[] initBottomDrops = {.15, .15, .15, .15, .15, .15, .15};

                    // double[] initTopDrops = {.50, .50, .50, .50, .50, .50, .50};
                    // double[] initBottomDrops = {.10, .10, .10, .10, .10, .10, .10};

                    double[] probesPerKmerMultFactors = initFileCache.getDoubles("AFFINITY MODELS", "ProbesPerKmerMultFactors");

                    //boolean[] includeRevCompsArray = {true, false, false, false, false};

                    int[] maxKmers = initFileCache.getInts("AFFINITY MODELS", "MaxKmersPerModel");
                    //int[] maxKmers = {4095, 4095, 4095, 4095, 4095, 4095, 4095};

                    boolean[] useSparseRegression = initFileCache.getBooleans("AFFINITY MODELS", "UseSparseKmerRegression");

                    fsam.kmerStartPositions = new int[kmerModelLengths.length][];
                    fsam.revCompMatrix = new char[kmerModelLengths.length][];
                    fsam.kmerToAffinityMatrix = new double[kmerModelLengths.length][];

                    //System.out.println("Hello 2");

                    for (int z=0; z < kmerModelLengths.length; z++) {

                        /////////////////////////////////////////
                        // Create the start positions for z
                        /////////////////////////////////////////
                        int[] psamStartPositions = getStartPositions(probeSeqLengths, fsam.kmerLengths[z]);
                        fsam.kmerStartPositions[z] = noOverHangStartPos(psamStartPositions);


                        /////////////////////////////////////////
                        // Set the revComp Martrix for z
                        /////////////////////////////////////////
                        KmerMatrix kmerMatrix = reduceData.getKmerMatrix(fsam.kmerLengths[z], 1, unifiedColumn);
                        fsam.revCompMatrix[z] = kmerMatrix.getRevCompMatrix();

                        /////////////////////////////////////////
                        // Set the includeRevComps for z
                        /////////////////////////////////////////
                        fsam.includeRevComps = includeRevCompsArray;

                        /////////////////////////////////////////
                        // Create the positional weights for z
                        /////////////////////////////////////////
                        double[][] noOverHangWeights = null;
//                         if (isLearnedPosBias[z]) {
//                         }
                        if (isUniformPosBias[z]) {
                            // create uniform 1.0 weights
                            noOverHangWeights = new double[2][fsam.kmerStartPositions[z].length];
                            Arrays.fill(noOverHangWeights[0], 1.0);
                            Arrays.fill(noOverHangWeights[1], 1.0);
                        }
                        else {
                            // create no-over-hang weights
                            noOverHangWeights = noOverHangWeights(psamStartPositions, fsam.kmerPositionalWeights[z]);
                        }

                        fsam.kmerPositionalWeights[z] = noOverHangWeights;
                        //fsam.symListToWordMap     = kmerMatrix.getSymListToWordMap();

                        /////////////////////////////////////////
                        // Make the words Array
                        /////////////////////////////////////////
                        char[] wordsArray = null;
                        //TCharHashSet wordsHS = null;
                        if (fsam.includeRevComps[z]) {
                            wordsArray = kmerMatrix.makeRevCompEquivKmers(fsam.kmerLengths[z]);
                            //wordsHS = makeRevCompEquivKmers(kmerLength);
                            //wordsArray = wordsHS.toArray();
                        }
                        else {
                            wordsArray = kmerMatrix.makeAllKmers(fsam.kmerLengths[z]);
                        }

                        /////////////////////////////////////////
                        // Create the Kmer to affinity matrix
                        /////////////////////////////////////////
                        out.println("\nCreating the \"All K-mer Affinity Model\" using "+((int)fsam.kmerLengths[z])+"-mers and the following positional-bias profile:");
                        out.println(""+StringTools.toString(fsam.kmerStartPositions[z], "\t"));
                        out.println(""+StringTools.toString(fsam.kmerPositionalWeights[z])+"....");

                        System.out.print("Creating Kmer to Probes lookup table...");

                        kmerMatrix.wordToProbesMatrix = kmerMatrix.makeKmerToProbesMatrix(
                            fsam.kmerLengths[z],
                            wordsArray,
                            fsam.includeRevComps[z]);

                        System.out.println("Done.");

                        if (useSparseRegression[z]) {
                        //if (fsam.kmerLengths[z] <= 4) {
                        //if (fsam.kmerLengths[z] <= 7) {
                        //if (!nonNegativeFlags[z]) {
                        //if (true) {
                        //if (false) {

                            fsam.kmerToAffinityMatrix[z] = fitSparseKmerModel(
                                fsam.kmerLengths[z],
                                fsam.kmerPositionalWeights[z],
                                wordsArray, // need to add to KmerMatrix object
                                kmerMatrix.wordToProbesMatrix,
                                residuals,
                                kmerMatrix,
                                //4000,
                                maxKmers[z],
                                fractionDrops[z],
                                fsam.includeRevComps[z],
                                true); // include intercept

                        }
                        else {

                            fsam.kmerToAffinityMatrix[z] = fitGradDescentKmerModel(
                                fsam.kmerLengths[z],
                                fsam.kmerPositionalWeights[z],
                                wordsArray, // need to add to KmerMatrix object
                                kmerMatrix.wordToProbesMatrix,
                                residuals,
                                kmerMatrix,
                                fractionDrops[z],
                                fsam.includeRevComps[z],
                                initValFlags[z],
                                itersArray[z],
                                initTopDrops[z],
                                initBottomDrops[z],
                                percentSteps[z],
                                probesPerKmerMultFactors[z],
                                nonNegativeFlags[z],
                                true); // include intercept
                        }

                        // DONE TRAINING!!
                        if (fsam.kmerToAffinityMatrix[z] != null) {
                            out.print("\n\tFinal Affinities = "
                                +"  "+format(fsam.kmerToAffinityMatrix[z][0], 3, 3)
                                +"  "+format(fsam.kmerToAffinityMatrix[z][1], 3, 3)
                                +"  "+format(fsam.kmerToAffinityMatrix[z][2], 3, 3)
                                +"  "+format(fsam.kmerToAffinityMatrix[z][3], 3, 3)
                                +"  "+format(fsam.kmerToAffinityMatrix[z][4], 3, 3)
                                +"  "+format(fsam.kmerToAffinityMatrix[z][5], 3, 3)
                                +"  "+format(fsam.kmerToAffinityMatrix[z][6], 3, 3)
                                +"  "+format(fsam.kmerToAffinityMatrix[z][7], 3, 3)
                                +"...");

                        }

                        // done training, now calculate residuals after applying
                        // this model
                        double[] affinityCorrections = kmerMatrix.getProbeAffinities(
                            fsam.kmerLengths[z],
                            fsam.kmerToAffinityMatrix[z],
                            fsam.kmerPositionalWeights[z],
                            fsam.includeRevComps[z]);

                        residuals = MathTools.subtract(residuals, affinityCorrections);

                        out.println("Done.");

                    }
                }


                ////////////////////////////////////////////////////////////////////////////////////////////////
                // Perform regression for current fsam!!
                ////////////////////////////////////////////////////////////////////////////////////////////////
                //                 rSquared = getRsquared(
                //                     fsam,
                //                     positionalWeights,
                //                     "fit with the current FSAM",
                //                     "lm",
                //                     strand,
                //                     strand,
                //                     calc,
                //                     null);

                ////////////////////////////////////////////////////////////////////////////////
                // Perform regression to find Non-Nearest-Neighbor dinucleotide dependencies
                ////////////////////////////////////////////////////////////////////////////////
                if (((String)initFileCache.get("AFFINITY MODELS", "Get_NNN_DinucleotideAffinities")).equalsIgnoreCase("Yes")) {

                    // Add all double point mutations to the motifToHammingSphereMap
                    // 0-based indexing
                    //for (int mut1=0; mut1 < motifLength-2; mut1++) { // from 0 to (L-1)
                    for (int mut1=0; mut1 < 5; mut1++) { // from 0 to (L-1)
                        //for (int mut1=0; mut1 < motifLength-1; mut1 +=2) { // from 0 to (L-1)

                        //for (int mut1=0; mut1 < 2; mut1++) {
                        //for (int mut1=0; mut1 < motifLength-1; mut1+=motifLength-2) {
                        //for (int mut1=0; mut1 < 1; mut1++) {

                        for (int mut2 = mut1 + 2; mut2 < motifLength; mut2++) { // Non-Nearest Neighbor
                            //for (int mut2 = mut1 + 1; mut2 <= mut1 + 1; mut2++) { // Nearest-Neighbor
                            //for (int mut2 = mut1 + 1; mut2 < motifLength; mut2++) { // All

                            motifsToModsMap = new LinkedHashMap<Object, Symbol[]>();
                            motifsToMandatoriesMap = new LinkedHashMap<Object, boolean[][]>();
                            motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();
                            aFeatureRemovalList = new ArrayList<FeatureKey>();


                            if (!featureDiscoveryAlgorithm.equalsIgnoreCase("lars")) {
                                if (useHammingDistForApprox) {
                                    motifToHammingSphereMap.put(hamming1PSAM, hammingSpheres.get(hammingDistForApprox));
                                }
                                else {
                                    motifToHammingSphereMap.put(hamming1PSAM, null);
                                }
                            }
                            //                     java.util.List avoidDimerResults = WeightMatrixTools.getMostLikelySymList(hamming1PSAM, 2);
                            //                     SymbolList avoidDimer =(SymbolList) avoidDimerResults.get(0);
                            //                     int avoidDimerStartPos = ((Integer)avoidDimerResults.get(1)).intValue();
                            //                     out.println("\nMost Likely Dimer is "+avoidDimer.seqString().toUpperCase()+" at start position "+(avoidDimerStartPos+1)+".");


                            boolean[][] mandatoryColumns = new boolean[2][motifLength];
                            mandatoryColumns[0][mut1] = true;
                            mandatoryColumns[0][mut2] = true;
                            mandatoryColumns[1] = ArrayTools.reverse(mandatoryColumns[0]);

                            // iterate over every symbol at the mut1 position
                            for (Iterator dnaIter2 = ((FiniteAlphabet)alphabet).iterator(); dnaIter2.hasNext(); ) {
                                Symbol mut1Symbol = (Symbol)dnaIter2.next();

                                // iterate over every symbol at the mut2 position
                                for (Iterator dnaIter3 = ((FiniteAlphabet)alphabet).iterator(); dnaIter3.hasNext(); ) {
                                    Symbol mut2Symbol = (Symbol)dnaIter3.next();

                                    // if the double-point mutation is different from the most likely nucleotides in the PSAM
                                    // then add them to the dinucleotide dependencies regression model

                                    Symbol avoidMut1Symbol = DistributionTools.getMostLikelySymbol(hamming1PSAM.getColumn(mut1));
                                    Symbol avoidMut2Symbol = DistributionTools.getMostLikelySymbol(hamming1PSAM.getColumn(mut2));


                                    //Symbol avoidMut1Symbol = DistributionTools.getLeastLikelySymbol(hamming1PSAM.getColumn(mut1));
                                    //Symbol avoidMut2Symbol = DistributionTools.getLeastLikelySymbol(hamming1PSAM.getColumn(mut2));
                                    //out.println((byte)100, "avoidMut1Symbol="+avoidMut1Symbol.getName().substring(0,1).toUpperCase()+"\tavoidMut2Symbol="+avoidMut2Symbol.getName().substring(0,1).toUpperCase());

                                    //                                     double avoidMut1Prob = hamming1PSAM.getColumn(mut1).getWeight(avoidMut1Symbol);
                                    //                                     double avoidMut2Prob = hamming1PSAM.getColumn(mut2).getWeight(avoidMut2Symbol);
                                    //                                     double dualProb = avoidMut1Prob * avoidMut2Prob;

                                    //                                     if (   ((mut1 == 2) && (mut1Symbol == DNATools.c()) && (mut2Symbol == DNATools.a()))
                                    //                                         || ((mut1 == 3) && (mut1Symbol == DNATools.a()) && (mut2Symbol == DNATools.c()))
                                    //                                         || ((mut1 == 4) && (mut1Symbol == DNATools.c()) && (mut2Symbol == DNATools.g()))
                                    //                                         || ((mut1 == 5) && (mut1Symbol == DNATools.g()) && (mut2Symbol == DNATools.t()))
                                    //                                         || ((mut1 == 6) && (mut1Symbol == DNATools.t()) && (mut2Symbol == DNATools.g()))
                                    //                                            ) {

                                    if (false) {
                                    //if ( (mut1Symbol == avoidMut1Symbol) && (mut2Symbol == avoidMut2Symbol) ) {

                                        //if ( (mut1 == avoidDimerStartPos) && (mut1Symbol == avoidDimer.symbolAt(1)) && (mut2Symbol == avoidDimer.symbolAt(2)) ) {
                                        //if ( (mut1Symbol == avoidMut1Symbol) && (mut2Symbol == avoidMut2Symbol) && (dualProb > .2) ) {

                                        out.println((byte)100, "\nNot including high affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression.");
                                        //out.println((byte)100, "\nNot including low affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression.");
                                        //out.println((byte)100, "\nNot including high affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression (dualProb = "+dualProb+", threshold = "+0.2+").");
                                    }
                                    else {

                                        WeightMatrix doubleMutationPSAM = WeightMatrixTools.mutate(hamming1PSAM, mut1, mut1Symbol, true);
                                        doubleMutationPSAM = WeightMatrixTools.mutate(doubleMutationPSAM, mut2, mut2Symbol, false);

                                        // name has 1-based indexing
                                        doubleMutationPSAM.setName(
                                            // "["+ (mut1+1) +"("+ mut1Symbol.getName().substring(0,1).toUpperCase()
                                            // + "), "+ (mut2+1) +"("+ mut2Symbol.getName().substring(0,1).toUpperCase() +")]");
                                            "mutate_"+ (mut1+1) +"_"+ mut1Symbol.getName().substring(0,1).toUpperCase()
                                            +"__"+ "mutate_"+ (mut2+1) +"_"+ mut2Symbol.getName().substring(0,1).toUpperCase() +"");

                                        Symbol[] modsArray = new Symbol[(motifLength*2)-1];
                                        modsArray[(mut1*2)] = mut1Symbol;
                                        modsArray[(mut2*2)] = mut2Symbol;

                                        LinkedHashSet<SymbolList> hammingSphere = null;

                                        if (useHammingDistForApprox) {

                                            hammingSphere = hammingSpheres.get(hammingDistForApprox);

                                            //                                         // create the 2 point-mutations in the seed
                                            //                                         // 1-based indexing
                                            //                                         SimpleSymbolList editedSeed = new SimpleSymbolList(seedSymList);
                                            //                                         Edit edit = new Edit(mut1+1, alphabet, mut1Symbol);
                                            //                                         editedSeed.edit(edit);
                                            //                                         edit = new Edit(mut2+1, alphabet, mut2Symbol);
                                            //                                         editedSeed.edit(edit);

                                            //                                         // create the hamming sphere around the double point-mutation
                                            //                                         hammingSphere = SymbolListTools.getHammingSphere(editedSeed, 1, includeRevComps, complementTable);
                                            //                                         out.println((byte)0, "hammingSphere = "+symListCollectionDump(hammingSphere));
                                        }

                                        //                                     out.println((byte)0, "WeightMatrix: "+WeightMatrixTools.toString(doubleMutationPSAM));
                                        //                                     out.println((byte)0, "**********************************************************");

                                        motifsToModsMap.put(doubleMutationPSAM, modsArray);
                                        motifToHammingSphereMap.put(doubleMutationPSAM, hammingSphere);
                                        motifsToMandatoriesMap.put(doubleMutationPSAM, mandatoryColumns);
                                        aFeatureRemovalList.add(new FeatureKey(modsArray));
                                    }
                                }
                            }

                            if (featureDiscoveryAlgorithm.equalsIgnoreCase("lars")) {

                                // perform penalized regression to find significant features
                                regressionFitData = fitModel(
                                    motifToHammingSphereMap,
                                    null, //startPositions
                                    positionalWeights,
                                    motifsToModsMap,
                                    motifsToMandatoriesMap,
                                    "non-nearest-neighbor dinucleotide dependencies",
                                    featureDiscoveryAlgorithm,
                                    //"lm",
                                    effectiveStrand,
                                    strand,
                                    //WeightMatrixTools.BindingStrand.BOTH,
                                    calc,
                                    null, //startParams
                                    null, //refSymList
                                    Double.NaN, //refKaReal);
                                    unifiedColumn);

                                //regressionFitData = fitModel(motifToHammingSphereMap, null, positionalWeights, motifsToModsMap, "nearest-neighbor dinucleotide dependencies", "lm", strand);

                                // coefficients are as returned in R
                                // First coefficient is the intercept, then the rest are in independent-var order
                                out.print(regressionFitData.toString());

                                //yeast
                                motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 5.0);
                                //motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 3.0);

                                //human
                                //motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 2.0);
                                //motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 1.8);

                                regressionFitData = fitModel(
                                    motifToHammingSphereMap,
                                    null, //startPositions
                                    positionalWeights,
                                    motifsToModsMap,
                                    motifsToMandatoriesMap,
                                    "non-nearest-neighbor dinucleotide dependencies",
                                    //featureDiscoveryAlgorithm,
                                    //"glm.poisson",
                                    //"rlm",
                                    "lm",
                                    effectiveStrand,
                                    strand,
                                    //WeightMatrixTools.BindingStrand.BOTH,
                                    calc,
                                    null, //startParams
                                    null, //refSymList
                                    Double.NaN, //refKaReal);
                                    unifiedColumn);

                            }

                            else {

                                ////////////////////////////////////////////////////////////////////////////////////////////////
                                // Perform univariate regression for current fsam!!
                                ////////////////////////////////////////////////////////////////////////////////////////////////
                                regressionFitData = getRegressionFit(
                                    fsam,
                                    positionalWeights,
                                    "fit with the current FSAM",
                                    "lm",
                                    //featureDiscoveryAlgorithm,
                                    effectiveStrand,
                                    strand,
                                    calc,
                                    unifiedColumn);

                                rSquared = regressionFitData.rSquared;
                                this.fsamUnivarCoeff = regressionFitData.getCoefficient(fsam);
                                this.fsamUnivarInter = regressionFitData.getIntercept();

                                nonSpecKa.setValue(getRelativeKaNonSpec(
                                        regressionFitData.getIntercept(),
                                        regressionFitData.getCoefficient(fsam),
                                        positionalWeights,
                                        eToMu,
                                        probeSeqLengths,
                                        motifLength,
                                        strand,
                                        calc,
                                        fsam.getRevCompSimilarity()));


                                // perform multivariate regression to fit with features

                                motifToFwmMap = fsam.clones(
                                    motifToHammingSphereMap,
                                    motifsToMandatoriesMap,
                                    aFeatureRemovalList);

                                psamOffset = true;
                                regressionFitData = fitModel(
                                    motifToHammingSphereMap,
                                    null, //startPositions
                                    positionalWeights,
                                    motifsToModsMap,
                                    motifsToMandatoriesMap,
                                    "non-nearest-neighbor dinucleotide dependencies",
                                    featureDiscoveryAlgorithm,
                                    //"lm",
                                    effectiveStrand,
                                    strand,
                                    //WeightMatrixTools.BindingStrand.BOTH,
                                    calc,
                                    null, //startParams
                                    null, //refSymList
                                    Double.NaN, //refKaReal);
                                    unifiedColumn);
                                psamOffset = false;

                            }

                            motifToFwmMap = null;

                            fsam.addFeatures(
                                !regressionFitData.isCoefficientZero(hamming1PSAM),
                                regressionFitData.getFeatures(hamming1PSAM, featuresScalar));
                            //fsam.setIntercept(regressionFitData.getIntercept(hamming1PSAM));

                            // coefficients are as returned in R
                            // First coefficient is the intercept, then the rest are in independent-var order
                            out.print(regressionFitData.toString());

                        }
                    }

                    out.print("\n"+fsam.toString());

                } // end - Perform regression to find Nearest-Neighbor dinucleotide dependencies


                ////////////////////////////////////////////////////////////////////////////////
                // Perform regression to find Nearest-Neighbor trinucleotide dependencies
                ////////////////////////////////////////////////////////////////////////////////
                if (((String)initFileCache.get("AFFINITY MODELS", "Get_NN_TrinucleotideAffinities")).equalsIgnoreCase("Yes")) {

                    // Add all triple point mutations to the motifToHammingSphereMap
                    // 0-based indexing
                    //for (int mut1=0; mut1 < motifLength-2; mut1++) { // from 0 to (L-1)
                    for (int mut1=0; mut1 < 5; mut1++) { // from 0 to (L-1)
                        //for (int mut1=0; mut1 < motifLength-1; mut1 +=2) { // from 0 to (L-1)

                        //for (int mut1=0; mut1 < 2; mut1++) {
                        //for (int mut1=0; mut1 < motifLength-1; mut1+=motifLength-2) {
                        //for (int mut1=0; mut1 < 1; mut1++) {

                        for (int mut2 = mut1 + 1; mut2 <= mut1 + 1; mut2++) { // Nearest-Neighbor
                            //for (int mut2 = mut1 + 1; mut2 < motifLength; mut2++) { // All

                            //for (int mut3 = mut2 + 1; mut3 < motifLength; mut3++) { // All
                            for (int mut3 = mut2 + 1; mut3 <= mut2 + 1; mut3++) { // Nearest-Neighbor

                                motifsToModsMap = new LinkedHashMap<Object, Symbol[]>();
                                motifsToMandatoriesMap = new LinkedHashMap<Object, boolean[][]>();
                                motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();
                                aFeatureRemovalList = new ArrayList<FeatureKey>();


                                if (!featureDiscoveryAlgorithm.equalsIgnoreCase("lars")) {
                                    if (useHammingDistForApprox) {
                                        motifToHammingSphereMap.put(hamming1PSAM, hammingSpheres.get(hammingDistForApprox));
                                    }
                                    else {
                                        motifToHammingSphereMap.put(hamming1PSAM, null);
                                    }
                                }
                                //                     java.util.List avoidDimerResults = WeightMatrixTools.getMostLikelySymList(hamming1PSAM, 2);
                                //                     SymbolList avoidDimer =(SymbolList) avoidDimerResults.get(0);
                                //                     int avoidDimerStartPos = ((Integer)avoidDimerResults.get(1)).intValue();
                                //                     out.println("\nMost Likely Dimer is "+avoidDimer.seqString().toUpperCase()+" at start position "+(avoidDimerStartPos+1)+".");


                                boolean[][] mandatoryColumns = new boolean[2][motifLength];
                                mandatoryColumns[0][mut1] = true;
                                mandatoryColumns[0][mut2] = true;
                                mandatoryColumns[0][mut3] = true;
                                mandatoryColumns[1] = ArrayTools.reverse(mandatoryColumns[0]);

                                // iterate over every symbol at the mut1 position
                                for (Iterator dnaIter2 = ((FiniteAlphabet)alphabet).iterator(); dnaIter2.hasNext(); ) {
                                    Symbol mut1Symbol = (Symbol)dnaIter2.next();

                                    // iterate over every symbol at the mut2 position
                                    for (Iterator dnaIter3 = ((FiniteAlphabet)alphabet).iterator(); dnaIter3.hasNext(); ) {
                                        Symbol mut2Symbol = (Symbol)dnaIter3.next();

                                        // iterate over every symbol at the mut3 position
                                        for (Iterator dnaIter4 = ((FiniteAlphabet)alphabet).iterator(); dnaIter4.hasNext(); ) {
                                            Symbol mut3Symbol = (Symbol)dnaIter4.next();

                                            // if the triple-point mutation is different from the most likely nucleotides in the PSAM
                                            // then add them to the dinucleotide dependencies regression model

                                            Symbol avoidMut1Symbol = DistributionTools.getMostLikelySymbol(hamming1PSAM.getColumn(mut1));
                                            Symbol avoidMut2Symbol = DistributionTools.getMostLikelySymbol(hamming1PSAM.getColumn(mut2));
                                            Symbol avoidMut3Symbol = DistributionTools.getMostLikelySymbol(hamming1PSAM.getColumn(mut3));


                                            //Symbol avoidMut1Symbol = DistributionTools.getLeastLikelySymbol(hamming1PSAM.getColumn(mut1));
                                            //Symbol avoidMut2Symbol = DistributionTools.getLeastLikelySymbol(hamming1PSAM.getColumn(mut2));
                                            //out.println((byte)100, "avoidMut1Symbol="+avoidMut1Symbol.getName().substring(0,1).toUpperCase()+"\tavoidMut2Symbol="+avoidMut2Symbol.getName().substring(0,1).toUpperCase());

                                            //                                     double avoidMut1Prob = hamming1PSAM.getColumn(mut1).getWeight(avoidMut1Symbol);
                                            //                                     double avoidMut2Prob = hamming1PSAM.getColumn(mut2).getWeight(avoidMut2Symbol);
                                            //                                     double dualProb = avoidMut1Prob * avoidMut2Prob;

                                            //                                     if (   ((mut1 == 2) && (mut1Symbol == DNATools.c()) && (mut2Symbol == DNATools.a()))
                                            //                                         || ((mut1 == 3) && (mut1Symbol == DNATools.a()) && (mut2Symbol == DNATools.c()))
                                            //                                         || ((mut1 == 4) && (mut1Symbol == DNATools.c()) && (mut2Symbol == DNATools.g()))
                                            //                                         || ((mut1 == 5) && (mut1Symbol == DNATools.g()) && (mut2Symbol == DNATools.t()))
                                            //                                         || ((mut1 == 6) && (mut1Symbol == DNATools.t()) && (mut2Symbol == DNATools.g()))
                                            //                                            ) {

                                            if (false) {
                                                //if ( !((mut1Symbol == avoidMut1Symbol) && (mut2Symbol == avoidMut2Symbol)&& (mut3Symbol == avoidMut3Symbol)) ) {
                                                //if ( (mut1Symbol == avoidMut1Symbol) && (mut2Symbol == avoidMut2Symbol) ) {

                                                //if ( (mut1 == avoidDimerStartPos) && (mut1Symbol == avoidDimer.symbolAt(1)) && (mut2Symbol == avoidDimer.symbolAt(2)) ) {
                                                //if ( (mut1Symbol == avoidMut1Symbol) && (mut2Symbol == avoidMut2Symbol) && (dualProb > .2) ) {

                                                out.println((byte)100, "\nNot including high affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression.");
                                                //out.println((byte)100, "\nNot including low affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression.");
                                                //out.println((byte)100, "\nNot including high affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression (dualProb = "+dualProb+", threshold = "+0.2+").");
                                            }
                                            else {

                                                WeightMatrix tripleMutationPSAM = WeightMatrixTools.mutate(hamming1PSAM, mut1, mut1Symbol, true);
                                                tripleMutationPSAM = WeightMatrixTools.mutate(tripleMutationPSAM, mut2, mut2Symbol, false);
                                                tripleMutationPSAM = WeightMatrixTools.mutate(tripleMutationPSAM, mut3, mut3Symbol, false);

                                                // name has 1-based indexing
                                                tripleMutationPSAM.setName(
                                                    //                                                 "["+ (mut1+1) +"("+ mut1Symbol.getName().substring(0,1).toUpperCase()
                                                    //                                                 + "), "+ (mut2+1) +"("+ mut2Symbol.getName().substring(0,1).toUpperCase()
                                                    //                                                 + "), "+ (mut3+1) +"("+ mut3Symbol.getName().substring(0,1).toUpperCase() +")]");
                                                    "mutate_"+ (mut1+1) +"_"+ mut1Symbol.getName().substring(0,1).toUpperCase()
                                                    +"__"+ "mutate_"+ (mut2+1) +"_"+ mut2Symbol.getName().substring(0,1).toUpperCase()
                                                    +"__"+ "mutate_"+ (mut3+1) +"_"+ mut3Symbol.getName().substring(0,1).toUpperCase() +"");

                                                Symbol[] modsArray = new Symbol[(motifLength*2)-1];
                                                modsArray[(mut1*2)] = mut1Symbol;
                                                modsArray[(mut2*2)] = mut2Symbol;
                                                modsArray[(mut3*2)] = mut3Symbol;

                                                LinkedHashSet<SymbolList> hammingSphere = null;

                                                if (useHammingDistForApprox) {

                                                    hammingSphere = hammingSpheres.get(hammingDistForApprox);

                                                    //                                         // create the 2 point-mutations in the seed
                                                    //                                         // 1-based indexing
                                                    //                                         SimpleSymbolList editedSeed = new SimpleSymbolList(seedSymList);
                                                    //                                         Edit edit = new Edit(mut1+1, alphabet, mut1Symbol);
                                                    //                                         editedSeed.edit(edit);
                                                    //                                         edit = new Edit(mut2+1, alphabet, mut2Symbol);
                                                    //                                         editedSeed.edit(edit);

                                                    //                                         // create the hamming sphere around the triple point-mutation
                                                    //                                         hammingSphere = SymbolListTools.getHammingSphere(editedSeed, 1, includeRevComps, complementTable);
                                                    //                                         out.println((byte)0, "hammingSphere = "+symListCollectionDump(hammingSphere));
                                                }

                                                //                                     out.println((byte)0, "WeightMatrix: "+WeightMatrixTools.toString(tripleMutationPSAM));
                                                //                                     out.println((byte)0, "**********************************************************");

                                                motifsToModsMap.put(tripleMutationPSAM, modsArray);
                                                motifToHammingSphereMap.put(tripleMutationPSAM, hammingSphere);
                                                motifsToMandatoriesMap.put(tripleMutationPSAM, mandatoryColumns);
                                                aFeatureRemovalList.add(new FeatureKey(modsArray));
                                            }
                                        }
                                    }
                                }

                                if (featureDiscoveryAlgorithm.equalsIgnoreCase("lars")) {

                                    // perform penalized regression to find significant features
                                    regressionFitData = fitModel(
                                        motifToHammingSphereMap,
                                        null, //startPositions
                                        positionalWeights,
                                        motifsToModsMap,
                                        motifsToMandatoriesMap,
                                        "nearest-neighbor trinucleotide dependencies",
                                        featureDiscoveryAlgorithm,
                                        //"lm",
                                        effectiveStrand,
                                        strand,
                                        //WeightMatrixTools.BindingStrand.BOTH,
                                        calc,
                                        null, //startParams
                                        null, //refSymList
                                        Double.NaN, //refKaReal);
                                        unifiedColumn);

                                    //regressionFitData = fitModel(motifToHammingSphereMap, null, positionalWeights, motifsToModsMap, "nearest-neighbor dinucleotide dependencies", "lm", strand);

                                    // coefficients are as returned in R
                                    // First coefficient is the intercept, then the rest are in independent-var order
                                    out.print(regressionFitData.toString());

                                    //yeast
                                    motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 5.0);
                                    //motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 3.0);

                                    //human
                                    //motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 2.0);
                                    //motifToHammingSphereMap = regressionFitData.significantMotifs(hamming1PSAM, 1.8);

                                    regressionFitData = fitModel(
                                        motifToHammingSphereMap,
                                        null, //startPositions
                                        positionalWeights,
                                        motifsToModsMap,
                                        motifsToMandatoriesMap,
                                        "nearest-neighbor trinucleotide dependencies",
                                        //featureDiscoveryAlgorithm,
                                        //"glm.poisson",
                                        //"rlm",
                                        "lm",
                                        effectiveStrand,
                                        strand,
                                        //WeightMatrixTools.BindingStrand.BOTH,
                                        calc,
                                        null, //startParams
                                        null, //refSymList
                                        Double.NaN, //refKaReal);
                                        unifiedColumn);

                                }

                                else {

                                    ////////////////////////////////////////////////////////////////////////////////////////////////
                                    // Perform univariate regression for current fsam!!
                                    ////////////////////////////////////////////////////////////////////////////////////////////////
                                    regressionFitData = getRegressionFit(
                                        fsam,
                                        positionalWeights,
                                        "fit with the current FSAM",
                                        "lm",
                                        //featureDiscoveryAlgorithm,
                                        effectiveStrand,
                                        strand,
                                        calc,
                                        unifiedColumn);

                                    rSquared = regressionFitData.rSquared;
                                    this.fsamUnivarCoeff = regressionFitData.getCoefficient(fsam);
                                    this.fsamUnivarInter = regressionFitData.getIntercept();

                                    nonSpecKa.setValue(getRelativeKaNonSpec(
                                            regressionFitData.getIntercept(),
                                            regressionFitData.getCoefficient(fsam),
                                            positionalWeights,
                                            eToMu,
                                            probeSeqLengths,
                                            motifLength,
                                            strand,
                                            calc,
                                            fsam.getRevCompSimilarity()));


                                    // perform multivariate regression to fit with features

                                    motifToFwmMap = fsam.clones(
                                        motifToHammingSphereMap,
                                        motifsToMandatoriesMap,
                                        aFeatureRemovalList);

                                    psamOffset = true;
                                    regressionFitData = fitModel(
                                        motifToHammingSphereMap,
                                        null, //startPositions
                                        positionalWeights,
                                        motifsToModsMap,
                                        motifsToMandatoriesMap,
                                        "nearest-neighbor trinucleotide dependencies",
                                        featureDiscoveryAlgorithm,
                                        //"lm",
                                        effectiveStrand,
                                        strand,
                                        //WeightMatrixTools.BindingStrand.BOTH,
                                        calc,
                                        null, //startParams
                                        null, //refSymList
                                        Double.NaN, //refKaReal);
                                        unifiedColumn);
                                    psamOffset = false;

                                }

                                motifToFwmMap = null;

                                fsam.addFeatures(
                                    !regressionFitData.isCoefficientZero(hamming1PSAM),
                                    regressionFitData.getFeatures(hamming1PSAM, featuresScalar));
                                //fsam.setIntercept(regressionFitData.getIntercept(hamming1PSAM));

                                // coefficients are as returned in R
                                // First coefficient is the intercept, then the rest are in independent-var order
                                out.print(regressionFitData.toString());

                            }
                        }
                    }

                    out.print("\n"+fsam.toString());

                } // end - Perform regression to find Nearest-Neighbor dinucleotide dependencies


                ////////////////////////////////////////////////////////////////////////////////
                // Perform regression to find spacer-length features
                ////////////////////////////////////////////////////////////////////////////////
                if (((String)initFileCache.get("AFFINITY MODELS", "GetSpacerLengthAffinities")).equalsIgnoreCase("Yes")) {

                    String spacerLenRegressionAlgorithm = (String)initFileCache.get("REGRESSION", "spacerLenRegressionAlgorithm");

                    motifsToMandatoriesMap = new LinkedHashMap<Object, boolean[][]>();
                    motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();

                    //                     if (useHammingDistForApprox) {
                    //                         motifToHammingSphereMap.put(hamming1PSAM, hammingSpheres.get(hammingDistForApprox));
                    //                     }
                    //                     else {
                    //                         motifToHammingSphereMap.put(hamming1PSAM, null);
                    //                     }

                    // Don't send in a hamming sphere ever
                    motifToHammingSphereMap.put(hamming1PSAM, null);

                    //Distribution allOnesDist = HMMTools.createUniformDistribution(alphabet,true);
                    Distribution allOnesDist = new SimpleDistribution(alphabet);
                    DistributionTools.allWeightsTo(allOnesDist, 1.0);

                    // Add all spacer-lengths between Min and Max inclusive
                    // 0-based indexing
                    int spacerLenMin = (int)initFileCache.getInt("AFFINITY MODELS", "SpacerLengthMin");
                    int spacerLenMax = (int)initFileCache.getInt("AFFINITY MODELS", "SpacerLengthMax");
                    int spacerPos = (int)Math.round(Math.floor((float) motifLength / 2));

                    for (int spacerLen = spacerLenMin; spacerLen <= spacerLenMax; spacerLen++) {

                        boolean[][] mandatoryColumns = new boolean[2][motifLength+spacerLen];
                        // include the positions neighboring the spacer as mandatory
                        Arrays.fill(mandatoryColumns[0], spacerPos-1, spacerPos+spacerLen+2, true);
                        mandatoryColumns[1] = ArrayTools.reverse(mandatoryColumns[0]);

                        // create uniform-dist spacer in the PSAM
                        WeightMatrix addedSpacerPSAM = WeightMatrixTools.insert(hamming1PSAM, spacerPos, spacerLen, allOnesDist);

                        // name has 1-based indexing
                        //addedSpacerPSAM.setName("[I("+ spacerPos +"-"+(spacerPos+1)+","+spacerLen+","+"N"+")]");
                        addedSpacerPSAM.setName("insert_"+ spacerPos +"_"+(spacerPos+1)+"_"+spacerLen+"_"+"N"+"");

                        //                         out.println((byte)0, "**********************************************************");
                        //                         out.println((byte)0, "WeightMatrix: "+WeightMatrixTools.toString(addedSpacerPSAM));
                        //                         out.println((byte)0, "**********************************************************");

                        // Don't send in a hamming sphere ever
                        motifsToMandatoriesMap.put(addedSpacerPSAM, mandatoryColumns);
                        motifToHammingSphereMap.put(addedSpacerPSAM, null);

                    }

                    // Add Positions to the RHS
                    this.addWeightsToLeft = false;
                    this.subtractOtherAffinities = true;

                    // perform linear regression to find coefficients for each spacer-length
                    regressionFitData = fitModel(
                        motifToHammingSphereMap,
                        null, //startPositions
                        positionalWeights,
                        null, //motifsToModsMap
                        motifsToMandatoriesMap,
                        "nucleotide-independent spacer-length features",
                        //linearAlgorithm,
                        spacerLenRegressionAlgorithm,
                        effectiveStrand,
                        strand,
                        calc,
                        null, //startParams
                        null, //refSymList
                        Double.NaN, //refKaReal);
                        unifiedColumn);

                    this.addWeightsToLeft = true;
                    this.subtractOtherAffinities = false;

                    // coefficients are as returned in R
                    // First coefficient is the intercept, then the rest are in independent-var order
                    out.print(regressionFitData.toString());

                    double[] coefs = regressionFitData.getCoefficients();
                    double[] coefsWithoutIntercept = Arrays.copyOfRange(coefs, 1, coefs.length);
                    double maxCoef = MathTools.max(coefsWithoutIntercept);
                    double[] relAffinities = MathTools.divide(coefsWithoutIntercept, maxCoef);

                    out.println("\nSpacer Length Relative Affiniities are: "+StringTools.toString(relAffinities, ", "));

                } //end - Perform regression to find spacer-length features



            } // end - for (int iterNum = 1; iterNum <= 3; iterNum++) {

            ////////////////////////////////////////////////////////////////////////////////////////////////
            // Perform regression for final fsam!!
            ////////////////////////////////////////////////////////////////////////////////////////////////
            rSquared = -1;
            noPosWeightsRSquared = -1;

            boolean performFinalFsamFit = initFileCache.getBoolean("REGRESSION", "PerformFinalFsamFit", true);
            if (performFinalFsamFit) {

                regressionFitData = getRegressionFit(
                    fsam,
                    positionalWeights,
                    "fit with the final FSAM",
                    "lm",
                    effectiveStrand,
                    strand,
                    calc,
                    unifiedColumn);

                rSquared = regressionFitData.rSquared;

//                 this.intercept = regressionFitData.getIntercept();
//                 this.intensitiesScaler = regressionFitData.getCoefficient(fsam);
//                 fsam.setScaler(regressionFitData.getCoefficient(fsam));

//                 nonSpecKa.setValue(getRelativeKaNonSpec(
//                         regressionFitData.getIntercept(),
//                         regressionFitData.getCoefficient(fsam),
//                         positionalWeights,
//                         eToMu,
//                         probeSeqLengths,
//                         motifLength,
//                         strand,
//                         calc,
//                         fsam.getRevCompSimilarity()));

//                 fsam.setNonSpecKa(nonSpecKa.doubleValue());

                regressionFitData = getRegressionFit(
                    fsam,
                    null,
                    "fit with the final FSAM and no posWeights",
                    "lm",
                    effectiveStrand,
                    strand,
                    calc,
                    unifiedColumn);

                noPosWeightsRSquared = regressionFitData.rSquared;
            }


            // end the R session
            rengine.end();

            if (saveAll) {
                WeightMatrixTools.writeToXML( hamming1PSAM, resultsDir + File.separator + proteinLabel+".init.psam."+motifLength+"nt.xml");
                FileTools.writeSerializedFile(hamming1PSAM, resultsDir + File.separator + proteinLabel+".init.psam."+motifLength+"nt.ser");

                WeightMatrix fsamPSAM = fsam.getPosStrandPWM();
                WeightMatrixTools.writeToXML( fsamPSAM, resultsDir + File.separator + proteinLabel+".final.psam."+motifLength+"nt.xml");
                FileTools.writeSerializedFile(fsamPSAM, resultsDir + File.separator + proteinLabel+".final.psam."+motifLength+"nt.ser");

                FileTools.writeSerializedFile(fsam, resultsDir + File.separator + proteinLabel+".fsam."+motifLength+"nt.ser");

                if (displayMotifs) {

                    GraphicsTools.writeImage(strandLogoFrame.getPanel(), resultsDir + File.separator + proteinLabel+".results."+motifLength+"nt.png", false);
                    GraphicsTools.writeImage(strandLogoFrame.getPanel(), resultsDir + File.separator + proteinLabel+".results."+motifLength+"nt.pdf", false);

                    //                 GridBagLayoutFrame finalLogoFrame = new GridBagLayoutFrame(proteinLabel+" Final PSAM", true);
                    //                 WeightMatrixLogo finalLogo = new WeightMatrixLogo(fsamPSAM, true, true, true, 0, WeightMatrixLogo.PositionLoc.BOTTOM, null);
                    //                 finalLogoFrame.add(finalLogo, 10, 10, false, true);
                    GridBagLayoutFrame finalLogoFrame = display(proteinLabel+" Final PSAM", fsamPSAM);
                    GraphicsTools.writeImage(finalLogoFrame.getPanel(), resultsDir + File.separator + proteinLabel+".psam."+motifLength+"nt.png", false);
                    GraphicsTools.writeImage(finalLogoFrame.getPanel(), resultsDir + File.separator + proteinLabel+".psam."+motifLength+"nt.pdf", false);

                    GridBagLayoutFrame fsamLogo = display(proteinLabel+" FSAM Logo", fsam, featureThresh);
                    GraphicsTools.writeImage(fsamLogo.getPanel(), resultsDir + File.separator + proteinLabel+".fsam."+motifLength+"nt.png", false);
                    GraphicsTools.writeImage(fsamLogo.getPanel(), resultsDir + File.separator + proteinLabel+".fsam."+motifLength+"nt.pdf", false);

                }

                FileTools.write(proteinLabel +"\t"+ seedSymList.seqString() +"\t"+ noPosWeightsRSquared +"\t"+ rSquared +"\n", resultsDir + File.separator + "unifiedMotifs."+motifLength+"nt.fsam.info", true);

                FileTools.write(                    "seedMotif" +"\t"+ "noPosWeightsRSquared" +"\t"+ "posWeightsRSquared" +"\n", resultsDir + File.separator + proteinLabel+"."+motifLength+"nt.fsam.info", false);
                FileTools.write(                    seedSymList.seqString() +"\t"+ noPosWeightsRSquared +"\t"+ rSquared +"\n", resultsDir + File.separator + proteinLabel+"."+motifLength+"nt.fsam.info", true);

                FileTools.write(fsam.toString(), resultsDir + File.separator + proteinLabel+".fsam."+motifLength+"nt.txt", false);

                if (positionalWeights != null) {
                    save(positionalWeights,
                        startPositions,
                        resultsDir + File.separator + proteinLabel+".psam."+motifLength+"nt.positionalBias.table");

                    // String[] rowLabels = {"posStrandRelBias", "negStrandRelBias"};
                    // String[] columnLabels = StringTools.toStrings(startPositions);
                    // String outputTable = StringTools.toString(positionalWeights, rowLabels, columnLabels);
                    // FileTools.write(
                    //     outputTable,
                    //     resultsDir + File.separator + proteinLabel+".psam."+motifLength+"nt.positionalBias.table",
                    //     false);
                }

                if ((fsam.kmerToAffinityMatrix != null)
                    && (fsam.kmerPositionalWeights != null)) {
                    for(int z=0; z < fsam.kmerPositionalWeights.length; z++) {
                        String[] rowLabels = {"posStrandRelBias", "negStrandRelBias"};
                        String[] columnLabels = StringTools.toStrings(fsam.kmerStartPositions[z]);
                        String outputTable = StringTools.toString(fsam.kmerPositionalWeights[z], rowLabels, columnLabels);
                        FileTools.write(
                            outputTable,
                            resultsDir + File.separator + proteinLabel+".kmers."+((int)fsam.kmerLengths[z])+"nt.positionalBias.table",
                            false);

                        if (displayMotifs) {
                            String[] datasetLabels = {"Smoothed Positive Strand", "Smoothed Negative Strand"};
                            double[][] values = {ArrayTools.toDoubleArray(fsam.kmerStartPositions[z]), fsam.kmerPositionalWeights[z][0], fsam.kmerPositionalWeights[z][1]};
                            BioJavaChart chart = new BioJavaChart("lines", null, "Start Position", "Normalized Coefficients", null, datasetLabels, true, values, -1);

                            GridBagLayoutFrame gridBagFrame = display("All Kmer Positional Bias Profile", chart);
                            GraphicsTools.writeImage(gridBagFrame.getPanel(), resultsDir + File.separator + proteinLabel+".kmers."+((int)fsam.kmerLengths[z])+"nt.positionalBias.png");
                            GraphicsTools.writeImage(gridBagFrame.getPanel(), resultsDir + File.separator + proteinLabel+".kmers."+((int)fsam.kmerLengths[z])+"nt.positionalBias.pdf");
                        }

                    }
                }

                if (fsam.kmerToAffinityMatrix != null) {
                    for(int z=0; z < fsam.kmerToAffinityMatrix.length; z++) {
                        String sortedKmerAffinitiesTable = KmerMatrix.getSortedTableString(fsam.kmerLengths[z], fsam.kmerToAffinityMatrix[z]);
                        FileTools.write(
                            sortedKmerAffinitiesTable,
                            resultsDir + File.separator + proteinLabel+".kmers."+((int)fsam.kmerLengths[z])+"nt.affinities.table",
                            false);
                    }
                }

                //FileTools.write(""+rSquared, resultsDir + File.separator + proteinLabel+".final.rSquared", false);
                //GraphicsTools.writeImage(psamsFrame.getPanel(), resultsDir + File.separator + proteinLabel+".psams.png");
            }

            return(fsam);
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(null);
    }

    public Table fitModel(
        String proteinAlignmentFileString,
        String affinityModelListFileString) {

        // protein sequences in the alignment must all be the same length

        return(null);
    }

    public RegressionFitData fitModel(
        LinkedHashSet motifs,
        int[] startPositions,
        double[][] positionalWeights,
        String motifsFamily,
        String method,
        WeightMatrixTools.BindingStrand effectiveStrand,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double[] startParams,
        SymbolList refSymList,
        double refKaReal,
        int unifiedColumn)
    {
        LinkedHashMap<Object, LinkedHashSet<SymbolList>> motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>(motifs.size());

        // set all the motifs to hash to null
        for (Object motif : motifs) {
            motifToHammingSphereMap.put(motif, null);
        }

        return(fitModel(
                motifToHammingSphereMap,
                startPositions,
                positionalWeights,
                null,
                null,
                motifsFamily,
                method,
                effectiveStrand,
                strand,
                calc,
                startParams,
                refSymList,
                refKaReal,
                unifiedColumn));
    }

    //
    // motifToHammingSphereMap has 3 different kind of entries
    // 1. SymbolList -> null
    // 2. WeightMatrix -> null
    // 3. WeightMatrix -> LinkedHashSet<SymbolList>
    //

    // Perform a fit for each experiment and get the average of the coefficients
    public RegressionFitData fitModel(
        LinkedHashMap<Object, LinkedHashSet<SymbolList>> motifToHammingSphereMap,
        int[] startPositions,
        double[][] positionalWeights,
        LinkedHashMap<Object, Symbol[]> aMotifsToModsMap,
        LinkedHashMap<Object, boolean[][]> aMotifsToMandatoriesMap,
        String motifsFamily,
        String method,
        WeightMatrixTools.BindingStrand effectiveStrand,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double[] startParams,
        SymbolList refSymList,
        double refKaReal,
        int unifiedColumn)
    {
        try {
            ArrayList<RegressionFitData> regressionFitDataArrayList = new ArrayList<RegressionFitData>();
            int numExperiments = reduceData.getNumExperiments();

            for (int experNum = 1; experNum <= numExperiments; experNum++) {

                RegressionFitData fitData = fitModel(
                    motifToHammingSphereMap,
                    startPositions,
                    positionalWeights,
                    aMotifsToModsMap,
                    aMotifsToMandatoriesMap,
                    motifsFamily,
                    method,
                    effectiveStrand,
                    strand,
                    calc,
                    startParams,
                    refSymList,
                    refKaReal,
                    unifiedColumn,
                    experNum);

                regressionFitDataArrayList.add(fitData);
            }

            // if there is just one experiment then just return that RegressionFitData
            if (regressionFitDataArrayList.size() == 1) {
                return(regressionFitDataArrayList.get(0));
            }
            else { // else average the coefficients
                return(RegressionFitData.average(regressionFitDataArrayList));
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(null);
    }

    public RegressionFitData fitModel(
        LinkedHashMap<Object, LinkedHashSet<SymbolList>> motifToHammingSphereMap,
        int[] startPositions,
        double[][] positionalWeights,
        LinkedHashMap<Object, Symbol[]> aMotifsToModsMap,
        LinkedHashMap<Object, boolean[][]> aMotifsToMandatoriesMap,
        String motifsFamily,
        String method,
        WeightMatrixTools.BindingStrand effectiveStrand,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double[] startParams,
        SymbolList refSymList,
        double refKaReal,
        int unifiedColumn,
        int experimentNumber)
    {
        //String linearAlgorithm = (String)initFileCache.get("REGRESSION", "LinearAlgorithm");
        //String featureDiscoveryAlgorithm = (String)initFileCache.get("REGRESSION", "FeatureDiscoveryAlgorithm");

        ArrayList<RegressionFitData> regressionFitDataArrayList = new ArrayList<RegressionFitData>();
        WeightMatrixTools.BindingStrand[] strandsArray = WeightMatrixTools.BindingStrand.values();
        RegressionFitData fitData;
        int maxStrandIndex;
        int startStrandIndex;

        // if motif is RC-palindrome then use only the positive strand
        if (strand == WeightMatrixTools.BindingStrand.POS) {
            startStrandIndex = 0;
            maxStrandIndex = 0;
        }
        else if (strand == WeightMatrixTools.BindingStrand.NEG) {
            startStrandIndex = 1;
            maxStrandIndex = 1;
        }
//         else if (isRevCompPalindrome) { // BOTH & SelfRevComp
//             startStrandIndex = 0;
//             maxStrandIndex = 0;
//         }
//         else if (method.equals(featureDiscoveryAlgorithm)) { // BOTH and featureDiscoveryAlgorithm
//             startStrandIndex = 2;
//             maxStrandIndex = 2;
//         }
        else { // BOTH and featureDiscoveryAlgorithm
//             startStrandIndex = 0;
//             maxStrandIndex = 1;
            startStrandIndex = 2;
            maxStrandIndex = 2;
        }

        // loop for the pos and maybe neg strand
        for (int strandIndex = startStrandIndex; strandIndex <= maxStrandIndex; strandIndex++) {
            WeightMatrixTools.BindingStrand currentStrand = strandsArray[strandIndex];

            fitData = fitModel(
                motifToHammingSphereMap,
                startPositions,
                positionalWeights,
                aMotifsToModsMap,
                aMotifsToMandatoriesMap,
                motifsFamily,
                method,
                currentStrand,
                calc,
                startParams,
                refSymList,
                refKaReal,
                unifiedColumn,
                experimentNumber);

            regressionFitDataArrayList.add(fitData);
        }

        if (regressionFitDataArrayList.size() == 1) {
            return(regressionFitDataArrayList.get(0));
        }
        else { // else average the coefficients
            RegressionFitData averagedRegressionFitData = RegressionFitData.average(regressionFitDataArrayList);
            if (strand == WeightMatrixTools.BindingStrand.BOTH) {
                averagedRegressionFitData.averageRevComps();
            }
            return(averagedRegressionFitData);
        }

    }

    // Perform the fit for a particular experiment
    //
    // If startPositions == null, then there are more than 1 motif
    // If startPositions != null, then there is just 1 motif at many different positions
    public RegressionFitData fitModel(
        LinkedHashMap<Object, LinkedHashSet<SymbolList>> motifToHammingSphereMap,
        int[] startPositions,
        double[][] positionalWeights,
        LinkedHashMap<Object, Symbol[]> aMotifsToModsMap,
        LinkedHashMap<Object, boolean[][]> aMotifsToMandatoriesMap,
        String motifsFamily,
        String method,
        //WeightMatrixTools.BindingStrand effectiveStrand,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double[] startParams,
        SymbolList refSymList,
        double refKaReal,
        int unifiedColumn,
        int experimentNumber)
    {

        StringBuffer regressionModel = new StringBuffer();
        StringBuffer cbindColumns = new StringBuffer();
        REXP rexp = null;

        try {

            final String frHome = System.getenv("FR_HOME");

            String rScriptsDir = frHome + File.separator + (String)initFileCache.get("REGRESSION", "RScriptsDir");
            String rNonlinearScript = (String)initFileCache.get("REGRESSION", "RNonlinearScript");
            String rPoissonScript = (String)initFileCache.get("REGRESSION", "RPoissonScript");
            String rBgSubtractScript = (String)initFileCache.get("REGRESSION", "RBgSubtractScript");

            String initPoolExpectedFreqsFile = (String) initFileCache.get("INTENSITIES", "InitPoolExpectedFreqsFile");

            double  minProbeAffinity                = initFileCache.getDouble("REGRESSION", "MinProbeAffinity");
            boolean backgroundSubtractIntensities   = (boolean)initFileCache.getBoolean("REGRESSION", "BackgroundSubtractIntensities");
            boolean removeBackgroundIntensityProbes = (boolean)initFileCache.getBoolean("REGRESSION", "RemoveBackgroundIntensityProbes");


            ArrayList<String> independentVarArrayList = new ArrayList<String>();
            ArrayList<double[]> motifCountsArrayList = new ArrayList<double[]>();

            String strandUsageString = null;
            switch (strand) {
            case POS:
            case NEG:
                strandUsageString = "the "+strand.toString()+" strand";
                break;
            case BOTH:
                strandUsageString = "BOTH strands with "+calc.toString()+" calc";
                break;
            }

            out.println("\n*******************************************************************************");
            out.print("Perparing to perform "+method+" regression in R to find the "+motifsFamily+" for experiment "+experimentNumber+".");
            out.println("\n*******************************************************************************");

//             LinkedHashSet motifsInRegression = new LinkedHashSet(motifToHammingSphereMap.size());
//             LinkedHashSet motifsNotInRegression = new LinkedHashSet(motifToHammingSphereMap.size());
            LinkedHashSet motifsInRegression = new LinkedHashSet();
            LinkedHashSet motifsNotInRegression = new LinkedHashSet();

            int probeSeqLengths = reduceData.getSeqLength();
            boolean isRevCompPalindrome = false;
            int numKmerWindows = -1;
            ArrayList<double[]> xTranspose = null;

            // call the multi-variate R routine with
            // 1. The probe intensity values
            // 2. Motif counts for each probe-associated sequence and each motif (weight matrix)

            // Get the intensities and set them to y
            double[] intensities = null;

            if (backgroundSubtractIntensities) {
                if (this.bgSubtractedIs == null) {
                    intensities = reduceData.getUsedIntensities(experimentNumber, unifiedColumn);
                    rengine.assign("intensities", intensities);
                    double[] minProbeAffinityArray = {minProbeAffinity};
                    rengine.assign("minProbeAffinity", minProbeAffinityArray);
                    //FileTools.write(intensities, "intensities.txt", false);

                    // load scripts
                    ArrayList<String> regressionStatements = new ArrayList<String>();

                    java.util.List<String> functions = getRStatements(rScriptsDir + File.separator + rBgSubtractScript);
                    regressionStatements.addAll(functions);
                    regressionStatements.add("bgSubtractedIs = getBgSubIntensities(intensities, minProbeAffinity)");
                    rengine.assign("FeatureReduceObject", rengine.createRJavaRef(this));

                    boolean displayReturnData = (boolean)initFileCache.getBoolean("REGRESSION", "DisplayReturnData");
                    out.println("\nEvaluating the following regression statements in order in R:");
                    for (String regressionStatement : regressionStatements) {
                        out.println(regressionStatement);
                        rexp = rengine.eval(regressionStatement);
                        if (displayReturnData) {
                            out.println("\t\treturn XT type is "+rexp.getType()+" : "+rexp.toString());
                        }
                    }
                    out.println("Done.");

                    this.bgSubtractedIs = rexp.asDoubleArray();

                    if (removeBackgroundIntensityProbes) {
                        removeProbesArray = MathTools.isEqualTo(this.bgSubtractedIs, minProbeAffinity);
                        int numTrue = MathTools.sum(true, removeProbesArray);
                        out.println(numTrue+" background subtracted intensities have intensity "+minProbeAffinity+".");
                    }

                    intensities = this.bgSubtractedIs;
                }
                else {
                    intensities = this.bgSubtractedIs;
                }
            }
            // No background subtraction
            else {
                intensities = reduceData.getUsedIntensities(experimentNumber, unifiedColumn);
            }

            // FileTools.write(intensities, "intensities.txt", false);

            // set intensities to residuals if needed
            if ((residualsFSAMs != null) && (!residualsFSAMs.isEmpty()) && (residualsSemaphore == false)) {
                residualsSemaphore = true;

                if (this.residuals == null) {
                    // FIX ME!!!!!! This only works for one FSAM, doesn't loop correctly
                    for (FeaturedWeightMatrix anFSAM : residualsFSAMs) {
                        RegressionFitData fitData = getRegressionFit(anFSAM, method, unifiedColumn);
                        out.println("Calculating the residuals from the fit of the fsam "+anFSAM.getName()+".");
                        //this.residuals = MathTools.lowerBound(fitData.getResiduals(), 0);
                        this.residuals = fitData.getResiduals();
                        out.println("Done. Will now fit to these new residuals.");
                        //out.println("Hello!!!");
                        // FileTools.write(this.residuals, "residuals.txt", false);
                    }
                }
                intensities = this.residuals;

                residualsSemaphore = false;
            }


            // remove probes if we are truncating probes at the loose end of the probes
            if (removeProbesArray != null) {
                int numPurge = MathTools.sum(true, removeProbesArray);
                out.print("\nRemoving "+numPurge+" probes from the intensities variable...");
                intensities = MathTools.purge(true, removeProbesArray, intensities);
                out.println(" Done.");
            }

            if (method.equalsIgnoreCase("nls.levMar")) {
                this.globalY = intensities;
                this.globalX = new double[4][][];
            }
            else {
                out.print("\nCreating the intensities variable in R...");
                rengine.assign("y", intensities);
                regressionModel.append("y ~");
                //rexp = rengine.eval("y");
                out.println(" Done.");
            }

            //////////////////////////////////////////////////////////////////////////////////////
            // Go through the motifs and get their counts and put them in X1, X2, X3,....
            //////////////////////////////////////////////////////////////////////////////////////
            int numMotifsInRegression = 0;
            int motifCount = 0;
            for (Object motif : motifToHammingSphereMap.keySet()) {

                //System.out.println("Here 1.");
                Object regressionMotif = motif;
                if (motifToFwmMap != null) {
                    regressionMotif = motifToFwmMap.get(motif);
                }

                String independentVar = null;
                int motifLength = -1;
                Object revCompMotif = null;

                if (motif instanceof WeightMatrix) {
                    independentVar = ((WeightMatrix)motif).getName();
                    motifLength = ((WeightMatrix)motif).columns();
                    revCompMotif = WeightMatrixTools.reverseComplement((WeightMatrix)motif, complementTable);
                }
                else if (motif instanceof FeaturedWeightMatrix) {
                    independentVar = ((FeaturedWeightMatrix)motif).getName();
                    motifLength = ((FeaturedWeightMatrix)motif).getPosStrandPWM().columns();
                    revCompMotif = null;
                }
                else { // it's a SymbolList
                    //////////////////////////////////////////////
                    // set method tp lm!!!!!!!!!!
                    //////////////////////////////////////////////
                    //method = "lm";

                    independentVar = ((SymbolList)motif).seqString();
                    motifLength = ((SymbolList)motif).length();
                    revCompMotif = DNATools.reverseComplement((SymbolList)motif);
                }
                numKmerWindows = probeSeqLengths + motifLength - 1;

                if (motif instanceof WeightMatrix) {
                    out.println((byte)0, "\nWeightMatrix: "+WeightMatrixTools.toString((WeightMatrix)motif));
                    out.println((byte)0, "RevComp: "+WeightMatrixTools.toString((WeightMatrix)revCompMotif)+"\n");
                }

                if (motifsInRegression.contains(motif)) {
                    out.println("An equivalent motif to the "+independentVar+" independent variable has already been added to the regression statement, not including "+independentVar+" in regression.");
                    motifsNotInRegression.add(motif);
                }
                //else if ((effectiveStrand == WeightMatrixTools.BindingStrand.BOTH) && (revCompMotif != null) && motifsInRegression.contains(revCompMotif)) {
                else if ((strand == WeightMatrixTools.BindingStrand.BOTH) && (revCompMotif != null) && motifsInRegression.contains(revCompMotif)) {
                    out.println("The reverse complement of the "+independentVar+" independent variable has already been added to the regression statement for both strands, not including "+independentVar+" in regression.");
                    motifsNotInRegression.add(motif);
                }

                ////////////////////////////////////////////////////////////////
                // add one variable with no specific startPos
                // we will use the positional weights and aMotifsToMandatoriesMap
                ////////////////////////////////////////////////////////////////
                else if (startPositions == null) {
                    //System.out.println("Here 2.");

                    LinkedHashSet<SymbolList> hammingSphere = null;
                    if (motifToHammingSphereMap != null) {
                        hammingSphere = motifToHammingSphereMap.get(motif);
                    }

                    boolean[][] mandatoryColumns = null;
                    if (aMotifsToMandatoriesMap != null) {
                        mandatoryColumns = aMotifsToMandatoriesMap.get(motif);
                    }

                    if (psamsPosWeights != null) {
                        positionalWeights = psamsPosWeights[motifCount];
                    }

                    double[][] modPositionalWeights = positionalWeights;

                    if (positionalWeights != null) {

                        // May need to remove posWeights from posBiasProfile if motifLength has gotten smaller
                        if (numKmerWindows < positionalWeights[0].length) {
                            System.out.println("modPositionalWeights!");
                            modPositionalWeights = new double[2][];

                            // Remove posWeights from the LHS of the posBiasProfile
                            if (this.addWeightsToLeft) {
                                modPositionalWeights[0] = Arrays.copyOfRange(positionalWeights[0], positionalWeights[0].length - numKmerWindows, positionalWeights[0].length);
                                modPositionalWeights[1] = Arrays.copyOfRange(positionalWeights[1], positionalWeights[1].length - numKmerWindows, positionalWeights[1].length);
                            }
                            // Remove posWeights from the RHS of the posBiasProfile
                            else {
                                modPositionalWeights[0] = Arrays.copyOfRange(positionalWeights[0], 0, numKmerWindows);
                                modPositionalWeights[1] = Arrays.copyOfRange(positionalWeights[1], 0, numKmerWindows);
                            }
                        }

                        // May need to pad 0s to the posBiasProfile if motifLength has gotten bigger
                        if (numKmerWindows > positionalWeights[0].length) {
                            System.out.println("modPositionalWeights!");
                            modPositionalWeights = new double[2][numKmerWindows]; //start with all 0s

                            // Pad 0s to the LHS of the posBiasProfile
                            if (this.addWeightsToLeft) {
                                modPositionalWeights[0] = ArrayTools.replaceRange(positionalWeights[0], 0, positionalWeights[0].length, modPositionalWeights[0], numKmerWindows - positionalWeights[0].length);
                                modPositionalWeights[1] = ArrayTools.replaceRange(positionalWeights[1], 0, positionalWeights[1].length, modPositionalWeights[1], numKmerWindows - positionalWeights[1].length);
                            }
                            // Pad 0s to the RHS of the posBiasProfile
                            else {
                                modPositionalWeights[0] = ArrayTools.replaceRange(positionalWeights[0], 0, positionalWeights[0].length, modPositionalWeights[0], 0);
                                modPositionalWeights[1] = ArrayTools.replaceRange(positionalWeights[1], 0, positionalWeights[1].length, modPositionalWeights[1], 0);
                            }
                        }
                        //out.println("numKmerWindows="+numKmerWindows+"; positionalWeights[0].length="+positionalWeights[0].length);
                    }

                    out.print("Retrieving the counts for the "+independentVar+" independent variable using "+strandUsageString+" ...");

                    double[] motifCounts =  reduceData.getCountsArray(
                        //motif,
                        regressionMotif,
                        hammingSphere,
                        null, //startPosition
                        mandatoryColumns,
                        //positionalWeights,
                        modPositionalWeights,
                        strand,
                        calc,
                        experimentNumber,
                        unifiedColumn);

                    out.println(" Done. The sum is "+MathTools.sum(motifCounts)+".");

                    if (MathTools.isAllZeros(motifCounts) || MathTools.hasInfiniteOrNaN(motifCounts)) {
                        out.println("The counts for the "+independentVar+" independent variable are all zeros or has Inf or NaN, trying again without positional weights.");
                        out.print("Retrieving the counts for the "+independentVar+" independent variable...");

                        motifCounts =  reduceData.getCountsArray(
                            //motif,
                            regressionMotif,
                            hammingSphere,
                            null, //startPosition
                            mandatoryColumns,
                            null, //positionalWeights,
                            strand,
                            calc,
                            experimentNumber,
                            unifiedColumn);

                        out.println(" Done. The sum is "+MathTools.sum(motifCounts)+".");
                    }

                    if (removeProbesArray != null) {
                        int numPurge = MathTools.sum(true, removeProbesArray);
                        out.print("\nRemoving "+numPurge+" probes from the "+independentVar+" motif counts...");
                        motifCounts = MathTools.purge(true, removeProbesArray, motifCounts);
                        out.println(" Done.");
                    }

                    // perform normalizing
//                     if ((refSymList != null) && (!Double.isNaN(refKaReal))) {
//                         double refKaFitted = getKaRef(motif, refSymList, strand, calc);
//                         //out.println("Before normalization, motifCounts[0] = "+motifCounts[0]);
//                         motifCounts = normalizeKaArray(motifCounts, refKaFitted, refKaReal);
//                         //out.println("After normalization, motifCounts[0] = "+motifCounts[0]);
//                     }


                    // int finalRoundTotalCounts = -1;
                    // if (method.equalsIgnoreCase("glm.poisson")) {
                    //     finalRoundTotalCounts = initFileCache.getInt("INTENSITIES", "FinalRoundTotalCounts");
                    // }

                    // check that (number of counts) = (number of intensities)
                    if ( motifCounts.length != intensities.length) {
                        out.println("\nError: number of counts for motif ("+motifCounts.length+") does not match number of intensities ("+intensities.length+").\n");
                        return(null);
                    }

                    if (MathTools.hasInfiniteOrNaN(motifCounts)) {
                        out.println("\nError: the counts for motif ("+motifCounts.length+") contains an infinite value or an NaN!\n");
                        return(null);
                    }

                    if (MathTools.isAllZeros(motifCounts)) {
                        out.println("The counts for the "+independentVar+" independent variable are all zeros, not including in regression.");
                        motifsNotInRegression.add(motif);
                    }
                    else {
                        //out.println("The counts sum for the "+independentVar+" independent variable is "+MathTools.sum(motifCounts)+".");

                        motifsInRegression.add(motif);

                        if (method.equalsIgnoreCase("nls.levMar")) {
                            xTranspose.add(motifCounts);
                        }
                        else {

                            if (method.equalsIgnoreCase("glm.poisson")) {
                                if (globalInitPoolFreqs == null) {
                                    globalInitPoolFreqs = reduceData.getOrderedInitPoolFreqs(experimentNumber, unifiedColumn);
                                    // globalExpectedNoEnrichCounts = MathTools.multiply(finalRoundTotalCounts, globalInitPoolFreqs);
                                }
                                //motifCounts = MathTools.multiply(motifCounts, globalExpectedNoEnrichCounts);
                                motifCounts = MathTools.multiply(motifCounts, globalInitPoolFreqs);
//                                 for (int i=0; i < motifCounts.length; i++) {
//                                     motifCounts[i] = motifCounts[i] * finalRoundTotalCounts * globalInitPoolFreqs[i];
//                                 }
                            }

                            if ((this.psamOffset) && (independentVar.equals("PSAM"))) {
                                motifCounts = MathTools.multiply(fsamUnivarCoeff, motifCounts);
                                //motifCounts = MathTools.multiply(this.fsamUnivarCoeff, motifCounts) + this.fsamUnivarInter;
                            }

                            //out.print("\nCreating the "+independentVar+" independent variable in R...");

                            rengine.assign(independentVar, motifCounts);

                            if (this.subtractOtherAffinities) {
                                independentVarArrayList.add(independentVar);
                                motifCountsArrayList.add(motifCounts);
                            }

                            //rexp = rengine.eval(independentVar);
                            //out.println(" Done.");

                            if (numMotifsInRegression == 0) {
                                if (this.psamOffset) {
                                    regressionModel.append(" offset(" + independentVar +")");
                                    cbindColumns.append(independentVar);
                                }
                                else {
                                    regressionModel.append(" " + independentVar);
                                    cbindColumns.append(independentVar);
                                }
                                //                                 regressionModel.append(" " + independentVar);
                                //                                 cbindColumns.append(independentVar);
                            }
                            else {
                                regressionModel.append(" + " + independentVar);
                                cbindColumns.append(", "+independentVar);
                            }

                            numMotifsInRegression++;

                        }
                    }

//                     // POS = 0, NEG = 1, BOTH = 2
//                     String origIndependentVar = independentVar;
//                     WeightMatrixTools.BindingStrand[] strandsArray = WeightMatrixTools.BindingStrand.values();

//                     boolean isRevCompPalindrome = isSelfReverseComplement(motif, complementTable);

//                     // look at both strands as default
//                     int maxStrandIndex;
//                     int startStrandIndex;

//                     // if motif is RC-palindrome then use only the positive strand
//                     if (strand == WeightMatrixTools.BindingStrand.POS) {
//                         startStrandIndex = 0;
//                         maxStrandIndex = 1;
//                     }
//                     else if (strand == WeightMatrixTools.BindingStrand.NEG) {
//                         startStrandIndex = 1;
//                         maxStrandIndex = 2;
//                     }
//                     else if (isRevCompPalindrome) { // BOTH & SelfRevComp
//                         startStrandIndex = 0;
//                         maxStrandIndex = 1;
//                     }
//                     else { // BOTH & !SelfRevComp
//                         startStrandIndex = 0;
//                         maxStrandIndex = 2;
//                     }

//                     // loop for the pos and neg strands as required
//                     for (int strandIndex = startStrandIndex; strandIndex < maxStrandIndex; strandIndex++) {

//                         WeightMatrixTools.BindingStrand currentStrand = strandsArray[strandIndex];
//                         independentVar = origIndependentVar.concat("_"+currentStrand.toString()+"strand");

//                         LinkedHashSet<SymbolList> hammingSphere = null;
//                         if (motifToHammingSphereMap != null) {
//                             hammingSphere = motifToHammingSphereMap.get(motif);
//                         }

//                         boolean[][] mandatoryColumns = null;
//                         if (aMotifsToMandatoriesMap != null) {
//                             mandatoryColumns = aMotifsToMandatoriesMap.get(motif);
//                         }

//                         out.print("Retrieving the counts for the "+independentVar+" independent variable using the "+currentStrand.toString()+" strand ...");

//                         double[] motifCounts =  reduceData.getCountsArray(
//                             motif,
//                             hammingSphere,
//                             null,
//                             mandatoryColumns,
//                             positionalWeights,
//                             currentStrand,
//                             calc,
//                             experimentNumber);

//                         out.println(" Done.");

//                         // check that (number of counts) = (number of intensities)
//                         if ( motifCounts.length != intensities.length) {
//                             out.println("\nError: number of counts for motif ("+motifCounts.length+") does not match number of intensities ("+intensities.length+").\n");
//                             return(null);
//                         }

//                         if (MathTools.hasInfiniteOrNaN(motifCounts)) {
//                             out.println("\nError: the counts for motif ("+motifCounts.length+") contains an infinite value or an NaN!\n");
//                             return(null);
//                         }

//                         if (MathTools.isAllZeros(motifCounts)) {
//                             out.println("The counts for the "+independentVar+" independent variable are all zeros, not including in regression.");
//                             motifsNotInRegression.add(motif);
//                         }
//                         else {
//                             motifsInRegression.add(motif);

//                             //out.print("\nCreating the "+independentVar+" independent variable in R...");
//                             rengine.assign(independentVar, motifCounts);
//                             //rexp = rengine.eval(independentVar);
//                             //out.println(" Done.");

//                             if (numMotifsInRegression == 0) {
//                                 regressionModel.append(" " + independentVar);
//                                 cbindColumns.append(independentVar);
//                             }
//                             else {
//                                 regressionModel.append(" + " + independentVar);
//                                 cbindColumns.append(", "+independentVar);
//                             }
//                             numMotifsInRegression++;
//                         }
//                     }

                }


                ////////////////////////////////////////////////////////////////
                // add a variable for each startPos using the same motif
                //
                // if we have non-null startPositions and >1 motif
                // then we are performing nonlinear column-reduce!!!
                //
                // we will not use the aMotifsToMandatoriesMap
                ////////////////////////////////////////////////////////////////
                else { // STARTPOSITIONS != NULL
                    //else if (motifToHammingSphereMap.size() == 1) {
                    // POS = 0, NEG = 1, BOTH = 2
                    String origIndependentVar = independentVar;
                    WeightMatrixTools.BindingStrand[] strandsArray = WeightMatrixTools.BindingStrand.values();

                    if (method.equalsIgnoreCase("nls.levMar")) {
                        xTranspose = new ArrayList<double[]>();
                        if (positionalWeights != null) {
                            this.globalConcGammas = ArrayTools.flatten(positionalWeights);
                        }
                    }

                    isRevCompPalindrome = isSelfReverseComplement(motif, complementTable);

                    // look at both strands as default
                    int maxStrandIndex;
                    int startStrandIndex;

                    // if motif is RC-palindrome then use only the positive strand
                    if (strand == WeightMatrixTools.BindingStrand.POS) {
                        startStrandIndex = 0;
                        maxStrandIndex = 1;
                    }
                    else if (strand == WeightMatrixTools.BindingStrand.NEG) {
                        startStrandIndex = 1;
                        maxStrandIndex = 2;
                    }
                    else if (isRevCompPalindrome) { // BOTH & SelfRevComp
                        startStrandIndex = 0;
                        maxStrandIndex = 1;
                    }
                    else { // BOTH & !SelfRevComp
                        startStrandIndex = 0;
                        maxStrandIndex = 2;
                    }

                    // loop for the pos and maybe neg strand
                    for (int strandIndex = startStrandIndex; strandIndex < maxStrandIndex; strandIndex++) {

                        // loop for each startPosition
                        for (int posIndex=0; posIndex < startPositions.length; posIndex++) {

                            WeightMatrixTools.BindingStrand currentStrand = strandsArray[strandIndex];

                            if (startPositions[posIndex] < 0) {
                                //out.println("\norigIndependentVar="+origIndependentVar);
                                //out.println("\ncurrentStrand="+currentStrand.toString());
                                independentVar = origIndependentVar.concat("_neg"+Math.abs(startPositions[posIndex])+"_"+currentStrand.toString()+"strand");
                            }
                            else {
                                independentVar = origIndependentVar.concat("_"+startPositions[posIndex]+"_"+currentStrand.toString()+"strand");
                            }

                            LinkedHashSet<SymbolList> hammingSphere = null;
                            if (motifToHammingSphereMap != null) {
                                hammingSphere = motifToHammingSphereMap.get(motif);
                            }

                            boolean[][] mandatoryColumns = null;
                            if (aMotifsToMandatoriesMap != null) {
                                mandatoryColumns = aMotifsToMandatoriesMap.get(motif);
                            }

                            out.print("Retrieving the counts for the "+independentVar+" independent variable...");

                            double[] motifCounts =  reduceData.getCountsArray(
                                motif,
                                hammingSphere,
                                startPositions[posIndex],
                                mandatoryColumns,
                                positionalWeights,
                                currentStrand,
                                calc,
                                experimentNumber,
                                unifiedColumn);

                            out.println(" Done. The sum is "+MathTools.sum(motifCounts)+".");

                            if (MathTools.isAllZeros(motifCounts) || MathTools.hasInfiniteOrNaN(motifCounts)) {
                                out.println("The counts for the "+independentVar+" independent variable are all zeros or has Inf or NaN, trying again without positional weights.");
                                out.print("Retrieving the counts for the "+independentVar+" independent variable...");

                                motifCounts =  reduceData.getCountsArray(
                                    motif,
                                    hammingSphere,
                                    startPositions[posIndex],
                                    mandatoryColumns,
                                    null, //positionalWeights,
                                    currentStrand,
                                    calc,
                                    experimentNumber,
                                    unifiedColumn);

                                out.println(" Done. The sum is "+MathTools.sum(motifCounts)+".");
                            }

                            if (removeProbesArray != null) {
                                int numPurge = MathTools.sum(true, removeProbesArray);
                                out.print("\nRemoving "+numPurge+" probes from the "+independentVar+" motif counts...");
                                motifCounts = MathTools.purge(true, removeProbesArray, motifCounts);
                                out.println(" Done.");
                            }

                            // perform normalizing
                            if ((refSymList != null) && (!Double.isNaN(refKaReal))) {
                                //out.println("Before normalization, motifCounts[0] = "+motifCounts[0]);
                                double refKaFitted = getKaRef(motif, refSymList, strand, calc);
                                motifCounts = normalizeKaArray(motifCounts, refKaFitted, refKaReal);
                                //out.println("After normalization, motifCounts[0] = "+motifCounts[0]);
                            }

                            // check that (number of counts) = (number of intensities)
                            if ( motifCounts.length != intensities.length) {
                                out.println("\nError: number of counts for motif ("+motifCounts.length+") does not match number of intensities ("+intensities.length+").\n");
                                return(null);
                            }

                            if (MathTools.hasInfiniteOrNaN(motifCounts)) {
                                out.println("\nError: the counts for motif ("+motifCounts.length+") contains an infinite value or an NaN!\n");
                                return(null);
                            }

                            if (MathTools.isAllZeros(motifCounts)) {
                                out.println("The counts for the "+independentVar+" independent variable are all zeros, not including in regression.");
                                motifsNotInRegression.add(motif);
                            }
                            else {
                                //out.println("The counts sum for the "+independentVar+" independent variable is "+MathTools.sum(motifCounts)+".");

                                motifsInRegression.add(motif);

                                if (method.equalsIgnoreCase("nls.levMar")) {
                                    xTranspose.add(motifCounts);
                                }
                                else {

                                    if (initPoolExpectedFreqsFile != null) {
                                    // if (method.equalsIgnoreCase("glm.poisson") || method.equalsIgnoreCase("nnls") || method.equalsIgnoreCase("rlm")) {
                                        if (globalInitPoolFreqs == null) {
                                            globalInitPoolFreqs = reduceData.getOrderedInitPoolFreqs(experimentNumber, unifiedColumn);
                                            // globalExpectedNoEnrichCounts = MathTools.multiply(finalRoundTotalCounts, globalInitPoolFreqs);
                                        }
                                        //motifCounts = MathTools.multiply(motifCounts, globalExpectedNoEnrichCounts);
                                        out.println("\nMultiplying motifCounts by globalInitPoolFreqs.");
                                        motifCounts = MathTools.multiply(motifCounts, globalInitPoolFreqs);
                                        //                                 for (int i=0; i < motifCounts.length; i++) {
                                        //                                     motifCounts[i] = motifCounts[i] * finalRoundTotalCounts * globalInitPoolFreqs[i];
                                        //                                 }
                                    }

                                    //out.print("\nCreating the "+independentVar+" independent variable in R...");
                                    rengine.assign(independentVar, motifCounts);
                                    //rexp = rengine.eval(independentVar);
                                    //out.println(" Done.");

                                    if (numMotifsInRegression == 0) {
//                                         if (psamOffset) {
//                                             regressionModel.append(" offset(" + independentVar +")");
//                                             cbindColumns.append(independentVar);
//                                         }
//                                         else {
//                                             regressionModel.append(" " + independentVar);
//                                             cbindColumns.append(independentVar);
//                                         }
                                        regressionModel.append(" " + independentVar);
                                        cbindColumns.append(independentVar);
                                    }
                                    else {
                                        regressionModel.append(" + " + independentVar);
                                        cbindColumns.append(", "+independentVar);
                                    }

                                    numMotifsInRegression++;

                                } //else !method.equalsIgnoreCase("nls.levMar")
                            } // else !MathTools.isAllZeros(motifCounts)

                        } // end for each posIndex
                    } //end for each strand

                    if (method.equalsIgnoreCase("nls.levMar")) {
                        double[][] xTransposeMatrix = xTranspose.toArray(new double[0][0]);
                        // FIX!!!!!!
                        //this.globalX[numMotifsInRegression] = MathTools.transpose(xTransposeMatrix);
                        this.globalX[0] = MathTools.transpose(xTransposeMatrix);
                    }

                } // end (add a var for each start positions

                motifCount++;
            } // end for(each motif)


            //////////////////////////////////////////////////////
            // Subtract affinities from each other
            //////////////////////////////////////////////////////
            if (this.subtractOtherAffinities) {
                for (int i=0; i < independentVarArrayList.size(); i++) {
                    String independentVar = independentVarArrayList.get(i);
                    double[] motifCounts = motifCountsArrayList.get(i);

                    //double trimmedMean = MathTools.trimmedMean(motifCounts, .0, .33); // drop the top 1/3
                    double minAffinity = MathTools.min(motifCounts);

                    for (int j=0; j < independentVarArrayList.size(); j++) {
                        // for all other independent vars perform subtraction
                        if (!independentVar.equals(independentVarArrayList.get(j))) {
                            motifCounts = MathTools.subtract(motifCounts, motifCountsArrayList.get(j));
                        }
                    }

                    // replace negative values with a positive constant
                    motifCounts = MathTools.lowerBound(motifCounts, minAffinity);

                    rengine.assign(independentVar, motifCounts);
                }
            }


            //////////////////////////////////////////////////////
            // Call the appropriate regression model
            // Methods = lm, rlm (robust), lqs (resistant), biglm (big models)
            //////////////////////////////////////////////////////
            boolean includeIntercept = (boolean)initFileCache.getBoolean("REGRESSION", "IncludeIntercept");

            ArrayList<String> regressionStatements = new ArrayList<String>();

            if (method.equalsIgnoreCase("lars") || method.equalsIgnoreCase("enet")) {
                boolean trace = (boolean)initFileCache.getBoolean("REGRESSION", "Trace");
                boolean normalize = (boolean)initFileCache.getBoolean("REGRESSION", "Normalize");
                int cvNum = (int)initFileCache.getInt("REGRESSION", "CvNum");
                this.globalLambda = (double)initFileCache.getDouble("REGRESSION", "Lambda");
                double epsScalar = (double)initFileCache.getDouble("REGRESSION", "EpsScalar");

                String traceString;
                String normalizeString;
                String interceptString;
                if (trace) {
                    traceString = "TRUE";
                }
                else {
                    traceString = "FALSE";
                }
                if (normalize) {
                    normalizeString = "TRUE";
                }
                else {
                    normalizeString = "FALSE";
                }
                if (includeIntercept) {
                    interceptString = "TRUE";
                }
                else {
                    interceptString = "FALSE";
                }

                //////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////
                // Testing!!!
                //interceptString = "FALSE";
                //////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////

                if (method.equalsIgnoreCase("lars")) {

                    regressionStatements.add("x = cbind("+cbindColumns.toString()+")");
                    regressionStatements.add("rm("+cbindColumns.toString()+")");
                    regressionStatements.add("model.cv = cv.lars(x, y, K="+cvNum+", normalize="+normalizeString+", trace="+traceString+", intercept="+interceptString+", use.Gram=TRUE, plot.it=FALSE, eps=.Machine$double.eps*"+epsScalar+")");
                    //regressionStatements.add("model.s = model.cv$fraction[order(model.cv$cv)[1]]");
                    regressionStatements.add("model.s = model.cv$fraction[which.min(model.cv$cv)]");
                    regressionStatements.add("model = lars(x, y, type=\"lasso\", normalize="+normalizeString+", trace="+traceString+", intercept="+interceptString+", use.Gram=TRUE, eps=.Machine$double.eps*"+epsScalar+")");
                    regressionStatements.add("coefsObject = predict.lars(model, s=model.s, type=\"coef\", mode=\"fraction\")");
                    regressionStatements.add("coefs = coefsObject$coefficients");
                }
                else { // enet!!!
                    regressionStatements.add("x = cbind("+cbindColumns.toString()+")");
                    regressionStatements.add("rm("+cbindColumns.toString()+")");
                    regressionStatements.add("model.cv = cv.enet(x, y, K="+cvNum+", lambda="+this.globalLambda+", s=seq(0,1,length=100), mode=\"fraction\", normalize="+normalizeString+", trace="+traceString+", intercept="+interceptString+", plot.it=FALSE, eps=.Machine$double.eps*"+epsScalar+")");
                    //regressionStatements.add("model.s = model.cv$fraction[order(model.cv$cv)[1]]");
                    //regressionStatements.add("model.cv$s");
                    //regressionStatements.add("model.cv$cv");
                    //regressionStatements.add("model.cv$cv.error");
                    regressionStatements.add("model.s = model.cv$s[which.min(model.cv$cv)]");

                    regressionStatements.add("model = enet(x, y, lambda="+this.globalLambda+", normalize="+normalizeString+", trace="+traceString+", intercept="+interceptString+", eps=.Machine$double.eps*"+epsScalar+")");

                    //regressionStatements.add("coefs = predict.enet(model, s=model.s, type=\"coef\", mode=\"fraction\")");
                    regressionStatements.add("coefsObject = predict.enet(model, s=model.s, type=\"coef\", mode=\"fraction\")");
                    regressionStatements.add("coefs = coefsObject$coefficients");
                }

            }
            else if (method.equalsIgnoreCase("nnls")) {
                if (includeIntercept) {
                    regressionStatements.add("x = cbind(1, "+cbindColumns.toString()+")");
                }
                else {
                    regressionStatements.add("x = cbind("+cbindColumns.toString()+")");
                }
                regressionStatements.add("rm("+cbindColumns.toString()+")");
                regressionStatements.add("model = nnls(x, y)");
                regressionStatements.add("coefs = coef(model)");
            }
            else if (method.equalsIgnoreCase("glm.poisson")) {
//                 java.util.List<String> functions = getRStatements(rScriptsDir + File.separator + rPoissonScript);
//                 regressionStatements.addAll(functions);

                if (globalInitPoolFreqs == null) {
                    // int finalRoundTotalCounts = initFileCache.getInt("INTENSITIES", "FinalRoundTotalCounts");

                    globalInitPoolFreqs = reduceData.getOrderedInitPoolFreqs(experimentNumber, unifiedColumn);
                    // globalExpectedNoEnrichCounts = MathTools.multiply(finalRoundTotalCounts, globalInitPoolFreqs);
                }

                if (includeIntercept) {

                    // double[] logExpectedNoEnrichCounts = MathTools.log(globalExpectedNoEnrichCounts);
                    // rengine.assign("expectedNoEnrichCounts", globalExpectedNoEnrichCounts);
                    // rengine.assign("logExpectedNoEnrichCounts", logExpectedNoEnrichCounts);
                    rengine.assign("initPoolFreqs", globalInitPoolFreqs);

                    ///////////////////////////////
                    // First multiplicative Model
                    ///////////////////////////////
//                     //regressionStatements.add("multModel = glm("+regressionModel+" + offset(logExpectedNoEnrichCounts), family=poisson(link=\"log\"), x=TRUE)");

//                     //includeIntercept = false;
//                     //regressionStatements.add("model = glm("+regressionModel+" + offset(logExpectedNoEnrichCounts) -1, family=poisson(link=\"log\"), x=TRUE)");

//                     //regressionStatements.add("model = glm("+regressionModel+" + log(initPoolFreqs), family=poisson(link=\"log\"), x=TRUE)");
//                     regressionStatements.add("multModel = glm("+regressionModel+" + log(initPoolFreqs), family=poisson(link=\"log\"), x=TRUE)");

//                     //regressionStatements.add("model = glm("+regressionModel+" + offset(log(initPoolFreqs)), family=poisson(link=\"log\"), x=TRUE)");

                    ///////////////////////////////
                    // Now Additive Model
                    ///////////////////////////////
                    // include the nonSpecIntercept at beginning of the model to model the intercept term
                    //double[] nonSpecIntercept = MathTools.multiply(numKmerWindows, globalExpectedNoEnrichCounts);


                    double[] nonSpecIntercept = MathTools.multiply(numKmerWindows, globalInitPoolFreqs);
                    rengine.assign("nonSpecIntercept", nonSpecIntercept);
                    regressionModel = regressionModel.replace(0, 3, "y ~ nonSpecIntercept +"); //replaces "y ~"



                    //regressionModel = regressionModel.append("+ nonSpecIntercept");

//                     //regressionStatements.add("newData = expectedNoEnrichCounts * multModel$x");
//                     regressionStatements.add("newData = data.frame(initPoolFreqs * multModel$x)");

//                     // drop the last term, then exponentiate the coeffs to go from the multiplicative model to additive
//                     //regressionStatements.add("multCoeffs = predict(multModel, type=\"response\")");
//                     regressionStatements.add("multCoeffs = multModel$coefficients");
//                     regressionStatements.add("startParams = exp(multCoeffs[-length(multCoeffs)])");
//                     //regressionStatements.add("startParams = exp(multModel$fitted[-length(multModel$fitted)])");

                    //regressionStatements.add("startParams = c(.000000001, 1.0, 1.0, 1.0, 1.0)");
                    //regressionStatements.add("startParams = c(1.0e7, 1.0e7, 1.0e7, 1.0e7)");
                    //regressionStatements.add("startParams = c(1.0, 1.0e7, 1.0e7, 1.0e7, 1.0e7)");
                    //regressionStatements.add("startParams = c(1.0e2, 1.0e7, 1.0e7, 1.0e7, 1.0e7)");
                    //regressionStatements.add("startParams = c(0.0, 1.0e7, 1.0e7, 1.0e7, 1.0e7, 1.0e2)");

                    if (startPositions == null) {
                        int repeatNum = numMotifsInRegression;
                        if (this.psamOffset) {
                            repeatNum--;
                        }
                        regressionStatements.add("startParams = c(1.0e4, rep.int(1.0e9, "+repeatNum+"))");
                        //regressionStatements.add("startParams = rep.int(1.0e9, "+repeatNum+")");
                    }
                    else {
                        regressionStatements.add("startParams = c(1.0e4, rep.int(1.0e9, "+2*startPositions.length+"))");
                        //regressionStatements.add("startParams = rep.int(1.0e9, "+2*startPositions.length+")");
                    }

                    //regressionStatements.add("model = glm("+regressionModel+" - 1, data=newData, family=poisson(link=\"identity\"), start=startParams)");
                    regressionStatements.add("model = glm("+regressionModel+" - 1, family=poisson(link=\"identity\"), start=startParams)");
                    //regressionStatements.add("model = glm("+regressionModel+", family=poisson(link=\"identity\"), start=startParams)");
                    //regressionStatements.add("model = glm("+regressionModel+" - 1, family=poisson(link=\"identity\"))");

//                     //regressionStatements.add("model = glm("+regressionModel+" - 1, family=poisson.additive, link=hyperbolic))");
//                     //regressionStatements.add("model = glm("+regressionModel+" - 1, family=poisson.additive(link=\"hyperbolic\"))");
//                     //regressionStatements.add("model = glm("+regressionModel+" - 1, family=poisson(link=\"identity\"))");
                }
                else { // don't include intercept
                    //regressionStatements.add("model = glm("+regressionModel+" - 1, family=poisson.additive, link=hyperbolic)");
                    //regressionStatements.add("model = glm("+regressionModel+" - 1, family=poisson.additive(link=\"hyperbolic\")");
                    //regressionStatements.add("model = glm("+regressionModel+" - 1, family=poisson(link=\"identity\"))");
                    //regressionStatements.add("model = glm("+regressionModel+" - 1, family=poisson(link=\"log\"))");
                    //regressionStatements.add("startParams = c(1.0e7, 1.0e7, 1.0e7, 1.0e7)");

                    if (startPositions == null) {
                        int repeatNum = numMotifsInRegression;
                        if (this.psamOffset) {
                            repeatNum--;
                        }
                        //regressionStatements.add("startParams = c(1.0e4, rep.int(1.0e9, "+repeatNum+"))");
                        regressionStatements.add("startParams =  rep.int(1.0e9, "+repeatNum+")");
                    }
                    else {
                        //regressionStatements.add("startParams = c(1.0e4, rep.int(1.0e9, "+2*startPositions.length+"))");
                        regressionStatements.add("startParams = rep.int(1.0e9, "+2*startPositions.length+")");
                    }
                    regressionStatements.add("model = glm("+regressionModel+" - 1, family=poisson(link=\"identity\"), start=startParams)");
                }

                regressionStatements.add("coefs = coef(model)");
                //regressionStatements.add("coefs = exp(coef(model))");

            }
            else if (method.equalsIgnoreCase("nls.levMar")) {
                //String functions = FileTools.readString(rScriptsDir + File.separator + "FeatureReduce.nonKa.levMar.r", false);
                java.util.List<String> functions = getRStatements(rScriptsDir + File.separator + rNonlinearScript);
                regressionStatements.addAll(functions);
                rengine.assign("FeatureReduceObject", rengine.createRJavaRef(this));

                double[] startParamsSqrt = MathTools.pow(startParams, .5);
                rengine.assign("startParams", startParamsSqrt);
                //rengine.assign("startParams", startParams);

                //regressionStatements.add("x = cbind("+cbindColumns.toString()+")");
                //regressionStatements.add("rm("+cbindColumns.toString()+")");

                this.globalNumWindows = numKmerWindows;
                this.globalIsSelfRevComp = isRevCompPalindrome;
                this.globalLambda = (double)initFileCache.getDouble("REGRESSION", "Lambda");

                if (startPositions != null) {
                    //regressionStatements.add("model = nls.lm(par=startParams, fn=getResidualGammaSqsL2, jac=getJacobianGammaSqsL2, control = nls.lm.control(nprint=1))");
                    //regressionStatements.add("model = nls.lm(par=startParams, fn=getResidualGammaSqsL1, jac=getJacobianGammaSqsL1, control = nls.lm.control(nprint=1, ptol=.01))");
                    //regressionStatements.add("model = nls.lm(par=startParams, fn=getResidualGammaSqsRg, jac=getJacobianGammaSqsRg, control = nls.lm.control(nprint=1))");
                    //regressionStatements.add("model = nls.lm(par=startParams, fn=getResidualGammaSqs, jac=getJacobianGammaSqs, control = nls.lm.control(nprint=1, ptol=.01))");
                    regressionStatements.add("model = nls.lm(par=startParams, fn=getResidualGammas, jac=getJacobianGammas, control = nls.lm.control(nprint=1))");
                }
                else {
                    regressionStatements.add("model = nls.lm(par=startParams, fn=getResidualKaSqs, jac=getJacobianKaSqs, control = nls.lm.control(nprint=1))");
                    //regressionStatements.add("model = nls.lm(par=startParams, fn=getResidualKas, jac=getJacobianKas, control = nls.lm.control(nprint=1))");
                }
                regressionStatements.add("coefs = coef(model)");
            }
            else if (method.equalsIgnoreCase("rlm")) {
                if (includeIntercept) {
                    regressionStatements.add("model = "+method+"("+regressionModel+", method = \"MM\")");
                    //regressionStatements.add("model = "+method+"("+regressionModel+")");
                }
                else {
                    regressionStatements.add("model = "+method+"("+regressionModel+" - 1, method = \"MM\")");
                    //regressionStatements.add("model = "+method+"("+regressionModel+" - 1)");
                }
                regressionStatements.add("coefs = coef(model)");
            }
            else if (method.equalsIgnoreCase("rlm.default")) {
                if (includeIntercept) {
                    regressionStatements.add("model = "+"rlm"+"("+regressionModel+")");
                    //regressionStatements.add("model = "+"rlm"+"("+regressionModel+")");
                }
                else {
                    regressionStatements.add("model = "+"rlm"+"("+regressionModel+" - 1)");
                    //regressionStatements.add("model = "+"rlm"+"("+regressionModel+" - 1)");
                }
                regressionStatements.add("coefs = coef(model)");
            }
            else { // just plain lm (linear model)
                if (includeIntercept) {
                    regressionStatements.add("model = "+method+"("+regressionModel+")");
                }
                else {
                    regressionStatements.add("model = "+method+"("+regressionModel+" - 1)");
                }
                regressionStatements.add("coefs = coef(model)");
            }




            /////////////////////////////////////////////////////////////////////////////////////////////
            // Call RLM multiple times and average (bootstrapping)
            /////////////////////////////////////////////////////////////////////////////////////////////
            ArrayList<RegressionFitData> fitDataArrayList = new ArrayList<RegressionFitData>();

            for (int i=0; i< numRegressions; i++) {


                boolean displayReturnData = (boolean)initFileCache.getBoolean("REGRESSION", "DisplayReturnData");
                out.println("\nEvaluating the following regression statements in order in R:");
                for (String regressionStatement : regressionStatements) {
                    out.println("\t"+regressionStatement);
                    rexp = rengine.eval(regressionStatement);
                    if (rexp == null) {
                        out.println("\nError: rengine.eval(regressionStatement) returned null.\n");
                        return(null);
                    }
                    else if (displayReturnData) {
                        out.println("\t\treturn XT type is "+rexp.getType()+" : "+rexp.toString());
                    }
                }

                out.println("Done.");

                double[] coefs = null;
                if (rexp != null) {
                    coefs = rexp.asDoubleArray();
                }

                if (coefs == null) {
                    out.println("\nError: could not access coefficients after model fit.\n");
                    return(null);
                }

                // add an intercept of 0 if it isn't present in the coefs
                if (!includeIntercept
                    || method.equalsIgnoreCase("lars")
                    || method.equalsIgnoreCase("enet")
                    || method.equalsIgnoreCase("nls.levMar")
                    //|| method.equalsIgnoreCase("glm.poisson")
                    ) {
                    //lars coefs does NOT include the intercept
                    // put 0 as the intercept in the first position
                    double[] newCoefs = new double[coefs.length+1];
                    newCoefs[0] = 0; // intercept set to 0
                    for (int jj=0; jj < coefs.length; jj++) {
                        newCoefs[jj+1] = coefs[jj];
                    }
                    coefs = newCoefs;
                }

                // add a PSAM coef if used OFFSET
                // PSAM coef is after the intercept
                if (psamOffset) {
                    coefs = ArrayUtils.add(coefs, 1, fsamUnivarCoeff);
                }

                // for NNLS must call lm to get the R2, RSS, etc.
                if (method.equalsIgnoreCase("nnls")) {
                    rexp = rengine.eval("predictions.nnls = as.vector(coefs %*% t(x))");
                    //rexp = rengine.eval("lmModel = lm(y ~ offset(predictions.nnls) - 1)");
                    //rexp = rengine.eval("lmModel = lm(y ~ predictions.nnls - 1)");
                    rexp = rengine.eval("lmModel = lm(y ~ predictions.nnls)");
                    rexp = rengine.eval("sum = summary(lmModel)");
                }
                else {
                    rexp = rengine.eval("sum = summary(model)");
                }

                double rSquared = Double.NaN;
                double stddev = Double.NaN;
                double rss = Double.NaN;
                double residuals[] = null;

                // get error statistics
                rexp = rengine.eval("residuals = model$residuals");
                residuals = rexp.asDoubleArray();

                if (method.equalsIgnoreCase("lm") || method.equalsIgnoreCase("glm.poisson")) {
                    rexp = rengine.eval("r_squared = sum$r.squared");
                    rSquared = rexp.asDouble();
                    rexp = rengine.eval("rss = sum(sum$residuals^2)");
                    rss = rexp.asDouble();
                }
                else if (method.equalsIgnoreCase("nls.levMar")) {
                    rexp = rengine.eval("r_squared = sum$r.squared");
                    rSquared = rexp.asDouble();
                    //rexp = rengine.eval("rss = sum(sum$residuals^2)");
                    //rss = rexp.asDouble();
                }
                else if (method.equalsIgnoreCase("nnls")) {
                    rexp = rengine.eval("r_squared = sum$r.squared");
                    rSquared = rexp.asDouble();
                    rexp = rengine.eval("rss = sum(sum$residuals^2)");
                    rss = rexp.asDouble();
                }
                else if (method.equalsIgnoreCase("rlm") || method.equalsIgnoreCase("rlm.default") || method.equalsIgnoreCase("lmrob")) {
                    rexp = rengine.eval("stddev = sum$stddev");
                    stddev = rexp.asDouble();
                }
                else if (method.equalsIgnoreCase("lars") || method.equalsIgnoreCase("enet")) {
                    rexp = rengine.eval("rss = sum$Rss");
                    rss = rexp.asDouble();
                    //rexp = rengine.eval("r2 = sum$R-squared");
                    //rSquared = rexp.asDouble();
                    //              rexp = rengine.eval("r2 = model$R2");
                }

                // Replace any NaN with 0
                coefs = ArrayTools.replaceNaNs(coefs, 0);

                // Replace any infinities
                // get min coef but don't include the intercept
                //double minCoef = MathTools.min(Arrays.copyOfRange(coefs, 1, coefs.length));
                //coefs = ArrayTools.replace(coefs, Double.POSITIVE_INFINITY, minCoef);
                //coefs = ArrayTools.replace(coefs, Double.NEGATIVE_INFINITY, minCoef);

                // square the coefs for certain nonlinear solver results
                if (method.equalsIgnoreCase("nls.levMar")) {
                    coefs = MathTools.pow(coefs, 2);
                }

                LinkedHashMap<Object, Double> motifsToCoefsMap = new LinkedHashMap<Object, Double>(motifsInRegression.size());
                double alpha = Double.NaN;
                double beta = Double.NaN;
                double nonSpecKa = Double.NaN;

                // if there are no startPositions then we need to populate the motifsToCoefsMap
                if (startPositions == null) {

                    // if nls.levMar then we need to set the extra params and remove them from coefs
                    // intercept has been added to first position

                    int k = 1; // skip the intercept
                    for (Object motif : motifsInRegression) {
                        motifsToCoefsMap.put(motif, coefs[k]);
                        k++;

                        //                     boolean isRevCompPalindrome = isSelfReverseComplement(motif, complementTable);
                        //                     if ((strand == WeightMatrixTools.BindingStrand.BOTH) && !isRevCompPalindrome) {
                        //                         // add the coefs for the positive and negative strands
                        //                         //motifsToCoefsMap.put(motif, coefs[k] + coefs[k+1]);
                        //                         motifsToCoefsMap.put(motif, (coefs[k] + coefs[k+1]) / 2);
                        //                         k += 2;
                        //                     }
                        //                     else {
                        //                         // there's just one coef for one strand
                        //                         motifsToCoefsMap.put(motif, coefs[k]);
                        //                         k++;
                        //                     }

                    }
                }
                // there are startPositions
                else {
                    // if nls.levMar then we need to set the extra params and remove them from coefs
                    // intercept has been added to first position
                    if (method.equalsIgnoreCase("nls.levMar")) {
                        alpha = coefs[1];
                        beta = coefs[2];
                        nonSpecKa = coefs[3];
                        coefs = ArrayUtils.remove(coefs, 1);
                        coefs = ArrayUtils.remove(coefs, 2);
                        coefs = ArrayUtils.remove(coefs, 3);
                    }
                }

                RegressionFitData regressionFitData = new RegressionFitData(
                    motifsToCoefsMap,
                    motifsNotInRegression,
                    aMotifsToModsMap,
                    aMotifsToMandatoriesMap,
                    coefs,
                    residuals,
                    alpha,
                    beta,
                    nonSpecKa,
                    startPositions,
                    rSquared,
                    stddev,
                    rss);


                if (regressionFitData != null) {
                    fitDataArrayList.add(regressionFitData);
                }

            }

            RegressionFitData regressionFitData = null;
            if (numRegressions == 1) {
                regressionFitData = fitDataArrayList.get(0);
            }
            else {
                regressionFitData = RegressionFitData.average(fitDataArrayList);
            }


            // free globals
            this.globalY = null;
            this.globalX = null;

            // remove all variables in R
            rexp = rengine.eval("rm(list=ls())");

            // garbage collect
            rexp = rengine.eval("gc()");

            return(regressionFitData);

        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(null);
    }




    double[] fitGradDescentKmerModel(
        int kmerLength,
        double[][] positionalBias,
        char[] wordsArray, // partial array of kmers
        char[][] kmerToProbesMatrix,  // partial array of kmer keys to probes
        double[] residuals,
        KmerMatrix kmerMatrix,
        double fractionDrop,
        boolean includeRevComps,
        String initValFlag,
        int numIters,
        double initTopDrop,
        double initBottomDrop,
        double percentStep,
        double probesPerKmerMultFactor,
        boolean nonNegativeFlag,
        boolean intercept)
    {
        double initValue = 0;
        if (initValFlag.equalsIgnoreCase("avg")) {

            initValue =
                MathTools.trimmedMean(
                    residuals,
                    initTopDrop,
                    initBottomDrop)

                / (kmerMatrix.numKmerWindows);

            if (includeRevComps) {
                initValue /= 2;
            }

            if (positionalBias != null) {
                //initValue *= MathTools.mean(positionalBias);
                //initValue /= MathTools.mean(positionalBias);
            }
        }
        out.println("\n\tIteration 0: initAffinity = "+initValue);

        double[] kmerToAffinityMatrix = null;

        int probesPerKmerStandard = (int) Math.ceil(kmerMatrix.getAvgProbesPerKmer(kmerMatrix.wordToProbesMatrix));

        for (int i=0; i < numIters; i++) {

            //out.println("\n\tIteration "+(i+1)+"");
            //                                 if (i < itersArray[z]/2) {
            //                                     percentStep /= ((itersArray[z]/2) - i);
            //                                     //percentStep = Math.pow(percentStep, .70);
            //                                 }

            kmerToAffinityMatrix = kmerMatrix.makeKmerToAffinityMatrix(
                kmerLength,
                wordsArray,
                // intensities,
                residuals,
                positionalBias,
                kmerMatrix.wordToProbesMatrix,
                kmerToAffinityMatrix,
                fractionDrop,
                nonNegativeFlag,
                percentStep,
                initValue,
                probesPerKmerStandard,
                probesPerKmerMultFactor,
                includeRevComps);

            if (kmerToAffinityMatrix != null) {
                out.print("\n\tIteration "+(i+1)+": step = "+format(percentStep,3,3)+" : Affinities[] = "
                    +"  "+format(kmerToAffinityMatrix[0], 3, 3)
                    +"  "+format(kmerToAffinityMatrix[1], 3, 3)
                    +"  "+format(kmerToAffinityMatrix[2], 3, 3)
                    +"  "+format(kmerToAffinityMatrix[3], 3, 3)
                    +"  "+format(kmerToAffinityMatrix[4], 3, 3)
                    +"  "+format(kmerToAffinityMatrix[5], 3, 3)
                    +"  "+format(kmerToAffinityMatrix[6], 3, 3)
                    +"  "+format(kmerToAffinityMatrix[7], 3, 3)
                    +"...");

            }

        }

        return(kmerToAffinityMatrix);
    }


    double[] fitSparseKmerModel(
        int kmerLength,
        double[][] positionalBias,
        char[] wordsArray, // partial array of kmers
        char[][] kmerToProbesMatrix,  // partial array of kmer keys to probes
        double[] residuals,
        KmerMatrix kmerMatrix,
        int maxColumns, // max columns (kmers) in the regression
        double fractionDrop,
        boolean includeRevComps,
        boolean intercept)
    {
        REXP rexp = null;

        int numWords = wordsArray.length;

        if ((kmerLength > 5) || (numWords > maxColumns)) {
            //if (numWords > maxColumns) {
            // get the top N trimmedMean-intensity K-mers
            out.print("\nCreating the full list of "+numWords+" "+kmerLength+"-mers sorted by abs(trimmed-mean(residuals))...");

            //int[] kmers = MathTools.sequence(0, wordsArray.length-1, 1);
            double[] kmerToTrimmedMeanArray = kmerMatrix.makeKmerToTrimmedMeanMatrix(
                kmerLength,
                residuals,
                positionalBias,
                wordsArray,
                kmerToProbesMatrix,
                fractionDrop,
                includeRevComps);

            // 0 = kmer
            // 1 = trimmedMean(residuals)
            // 2 = abs(trimmedMean(residuals))
            //Table kmersToTrimmedMeanTable(kmers, kmerToTrimmedMeanArray);
            Table kmersToTrimmedMeanTable = new Table();
            for (int index=0; index < wordsArray.length; index++) {

                java.util.List keyValueList = Arrays.asList(
                    new Integer(wordsArray[index]),
                    new Double(kmerToTrimmedMeanArray[index]),
                    new Double(Math.abs(kmerToTrimmedMeanArray[index])),
                    new Integer(kmerToProbesMatrix[index].length),
                    kmerToProbesMatrix[index]);

                kmersToTrimmedMeanTable.add(keyValueList);
            }

            out.println(" Done.");

            int probeRemovalThresh = 10;

            // Remove any of the rows with probe-counts < probeRemovalThresh
            out.print("\nRemoving any of the "+kmerLength+"-mers with probe counts < "+probeRemovalThresh+"...");
            kmersToTrimmedMeanTable.sort(3); //sort by ascending abs(probeCounts)
            int searchIndex = kmersToTrimmedMeanTable.binarySearch(3, new Integer(probeRemovalThresh));
            if (true) {
                //if (searchIndex < 0) {
                // if not found then value = (-(insertion point) - 1)
                kmersToTrimmedMeanTable.removeRows(0, Math.abs(searchIndex+1));
            }
            out.println(" Done.");
            out.println("\nRemoved "+(Math.abs(searchIndex+1))+" "+kmerLength+"-mers with probe counts < "+probeRemovalThresh+".");

            // Remove any of the rows with trimmed-mean residuals close to zero
            out.print("\nRemoving any of the "+kmerLength+"-mers with trimmed-mean residuals close to zero...");
            kmersToTrimmedMeanTable.sort(2); //sort by ascending abs(trimmedMean)
            searchIndex = kmersToTrimmedMeanTable.binarySearch(2, new Double(5.0));
            if (true) {
                //if (searchIndex < 0) {
                // if not found then value = (-(insertion point) - 1)
                kmersToTrimmedMeanTable.removeRows(0, Math.abs(searchIndex+1));
            }
            out.println(" Done.");
            out.println("\nRemoved "+(Math.abs(searchIndex+1))+" "+kmerLength+"-mers with trimmed-mean residuals close to zero.");

            // Now get the top N of what's left
            kmersToTrimmedMeanTable.sort(1); //sort by ascending trimmedMean
            kmersToTrimmedMeanTable.reverse(); // get descending order
            if (kmersToTrimmedMeanTable.rows() > maxColumns) {
                //out.print("\nCreating the sorted list of the top "+maxColumns+" "+kmerLength+"-mers with the highest abs(trimmed-mean(residuals))...");
                out.print("\nCreating the sorted list of the top "+maxColumns+" "+kmerLength+"-mers with the highest trimmed-mean(residuals)...");
                 // remove all but the top N
                kmersToTrimmedMeanTable.removeRows(maxColumns, kmersToTrimmedMeanTable.rows());
                out.println(" Done.");
            }

            if (kmersToTrimmedMeanTable.rows() >= 5) {

                out.println("\nFirst 5 highest TrimmedMeans = "
                    +kmersToTrimmedMeanTable.getElement(0, 1)+", "
                    +kmersToTrimmedMeanTable.getElement(1, 1)+", "
                    +kmersToTrimmedMeanTable.getElement(2, 1)+", "
                    +kmersToTrimmedMeanTable.getElement(3, 1)+", "
                    +kmersToTrimmedMeanTable.getElement(4, 1)
                            );

                out.println("\nLast 5 lowest TrimmedMeans = "
                    +kmersToTrimmedMeanTable.getElement(kmersToTrimmedMeanTable.rows()-5, 1)+", "
                    +kmersToTrimmedMeanTable.getElement(kmersToTrimmedMeanTable.rows()-4, 1)+", "
                    +kmersToTrimmedMeanTable.getElement(kmersToTrimmedMeanTable.rows()-3, 1)+", "
                    +kmersToTrimmedMeanTable.getElement(kmersToTrimmedMeanTable.rows()-2, 1)+", "
                    +kmersToTrimmedMeanTable.getElement(kmersToTrimmedMeanTable.rows()-1, 1)
                            );
            }

            java.util.List kmerIndecesList = kmersToTrimmedMeanTable.getColumn(0);
            java.util.List kmerToProbesMatrixList = kmersToTrimmedMeanTable.getColumn(4);

             // will not change wordsArray in caller!!!
            wordsArray = ArrayTools.toCharArray((ArrayList<Integer>)kmerIndecesList);

            // will not change kmerToProbesMatrix in caller!!!
            kmerToProbesMatrix = ArrayTools.toCharMatrix((ArrayList)kmerToProbesMatrixList);
            //kmerToProbesMatrix = (char[][])kmerToProbesMatrixList.toArray();

        }

        for (int regressionLoop = 1; regressionLoop <= 20; regressionLoop++) {

            out.println("\n"+wordsArray.length+" "+kmerLength+"-kmers will be used to contruct the sparse design matrix.");

            //[0] = row, [1] = column, [2] = value
            double[][] sparseKmerMatrix = kmerMatrix.makeSparseKmerMatrix(
                kmerLength,
                wordsArray,
                kmerToProbesMatrix,
                positionalBias,
                residuals,
                includeRevComps,
                intercept);

            ArrayList<String> regressionStatements = new ArrayList<String>();

            out.print("\nCreating the sparse matrix in R...");
            rengine.assign("rows", sparseKmerMatrix[0]);
            rengine.assign("columns", sparseKmerMatrix[1]);
            rengine.assign("values", sparseKmerMatrix[2]);
            rengine.assign("y", sparseKmerMatrix[3]);
            out.println(" Done.");

            regressionStatements.add("sMatrix = sparseMatrix(rows, columns, x=values)");

            // cholesky requires t(sMatrix) ?
            // qr requires sMatrix ?
            if (kmerLength <= 5) {

                //////////////////////////////////////////////////////////////////////////
                // QR-factorization = more accurate, slower
                //////////////////////////////////////////////////////////////////////////
                //regressionStatements.add("model = Matrix:::lm.fit.sparse(sMatrix, y)", method=\"qr\")");
                //regressionStatements.add("coefs = model$coef");

                //////////////////////////////////////////////////////////////////////////
                // Cholesky -factorization = less accurate, faster
                //////////////////////////////////////////////////////////////////////////
                //regressionStatements.add("model = Matrix:::lm.fit.sparse(t(sMatrix), y, method=\"cholesky\")");
                //regressionStatements.add("coefs = model$coef");

                regressionStatements.add("model = solve(crossprod(sMatrix), crossprod(sMatrix, y))");
                regressionStatements.add("coefs = model[,1]");

            }
            else {
                //////////////////////////////////////////////////////////////////////////
                // QR-factorization = more accurate, slower
                //////////////////////////////////////////////////////////////////////////
                //regressionStatements.add("model = Matrix:::lm.fit.sparse(sMatrix, y, method=\"qr\")");
                //regressionStatements.add("coefs = model$coef");

                //////////////////////////////////////////////////////////////////////////
                // Cholesky-factorization = less accurate, faster
                //////////////////////////////////////////////////////////////////////////
                // regressionStatements.add("model = Matrix:::lm.fit.sparse(t(sMatrix), y, method=\"cholesky\")");
                // regressionStatements.add("coefs = model$coef");

                regressionStatements.add("model = solve(crossprod(sMatrix), crossprod(sMatrix, y))");
                regressionStatements.add("coefs = model[,1]");
            }


            boolean displayReturnData = (boolean)initFileCache.getBoolean("REGRESSION", "DisplayReturnData");
            out.println("\nEvaluating the following regression statements in order in R:");
            for (String regressionStatement : regressionStatements) {
                out.println("\t"+regressionStatement);
                rexp = rengine.eval(regressionStatement);
                if (rexp == null) {
                    out.println("\nError: rengine.eval(regressionStatement) returned null.\n");
                    //return(null);
                    break;
                }
                else if (displayReturnData) {
                    out.println("\t\treturn XT type is "+rexp.getType()+" : "+rexp.toString());
                }
            }

            if (rexp == null) {
                out.println("\nCholesky Sparse Regression Failed: trying again with 5% less "+kmerLength+"-mers.\n");
                int tenPercentLess = (int)Math.round(wordsArray.length *.95);
                wordsArray = Arrays.copyOfRange(wordsArray, 0, tenPercentLess);
            }
            else {
                break;
            }

        }
        out.println("Done.");

        double[] coefs = null;
        if (rexp != null) {
            coefs = rexp.asDoubleArray();
        }

        // remove all variables in R
        rexp = rengine.eval("rm(list=ls())");

        // garbage collect
        rexp = rengine.eval("gc()");

        if (coefs == null) {
            out.println("\nError: could not access coefficients after model fit.\n");
            return(null);
        }

        // if using all the Kmers in the regression then just return all the coefs
        int maxWord = KmerMatrix.multBy4ToPower(1, kmerLength);
        if (wordsArray.length == maxWord) {
            return(coefs);
        }

        // add +1 to array length for the intercept
        double[] fullKmerToAffinityArray = null;
        if (intercept) {
            fullKmerToAffinityArray = new double[maxWord+1];
        }
        else {
            fullKmerToAffinityArray = new double[maxWord];
        }

        // copy over all the kmer affinities
        for (int index=0; index < wordsArray.length; index++) {
            char word = wordsArray[index];
            fullKmerToAffinityArray[word] = coefs[index];

            if (includeRevComps) {
                char revComp = kmerMatrix.revCompMatrix[word];
                fullKmerToAffinityArray[revComp] = coefs[index];
            }

        }

        // copy over the intercept
        if (intercept) {
            fullKmerToAffinityArray[maxWord] = coefs[wordsArray.length];
        }

        return(fullKmerToAffinityArray);
    }


    //////////////////////////////////////////////////////////////////////////
    // parList is in the format (alpha, beta, nonSpecKa, concGamma1pos, concGamma2pos, ...., concGamma1neg, ....)
    // this.globalX is in the format (posStrandKa(1), posStrandKa(2), ...., negStrandKa(1), .....) for each row (probe)
    // lambda is the tuner for the L2 regularization ( - lambda[ (concGamma_I_pos - concGamma_I+1_pos)^2 ] )
    // alpha is an overall signal intensity scaler (models the light source, reflection intensity, light meter sensitivity)
    // beta is an overall signal intensity constant (intercept)
    //
    // return value is really a matrix forced into an array
    // return array is in "banded jacobian" format = (d(e1,v1), d(e2,v1), ....., d(e1,v2), d(e2,v2), ......)
    //
    // if isSelfRevComp then look only at the positive strand
    //////////////////////////////////////////////////////////////////////////
    public double[] getJacobianGammas(double[] parList) {
//         System.out.println("parList = "+StringTools.toString(parList, "\t"));
//         System.out.println("y.length = "+y.length);
//         System.out.println("x.length = "+x.length);
//         System.out.println("x[0].length = "+x[0].length);
//         System.out.println("numWindows = "+numWindows);
//         System.out.println("isSelfRevComp = "+Boolean.toString(isSelfRevComp));

        int numPreVars = 3;
        double[][] kaMatrix = this.globalX[0];
        double[] jacobian = new double[kaMatrix.length*(numPreVars + (2*this.globalNumWindows))];
        //System.out.println("jacobian[0] = "+StringTools.toString(Arrays.copyOfRange(jacobian, 0, numPreVars+this.globalNumWindows), "\t"));
        double alpha = parList[0];
        double beta = parList[1];
        double nonSpecKa = parList[2];

        if (!this.globalIsSelfRevComp) {
            for (int row=0; row < kaMatrix.length; row++) {
//                 if ((row % 5000) == 0) {
//                     System.out.println("kaMatrix.length = "+kaMatrix.length+" : current row = "+row);
//                 }

                for (int i=0; i < this.globalNumWindows; i++) {
                    double concGammaPos = parList[numPreVars + i];
                    double kaPos = kaMatrix[row][i];
                    double concGammaNeg = parList[numPreVars + this.globalNumWindows + i];
                    double kaNeg = kaMatrix[row][this.globalNumWindows + i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns);
                    double bigProd = (concGammaPos*kaPos) + (concGammaNeg*kaNeg) + ((concGammaPos+concGammaNeg)*nonSpecKa);
                    double bigProdPlus1Sq = Math.pow((bigProd +1),2);

                    // d(alpha)
                    jacobian[row] += -1 * (bigProd / (bigProd +1));

                    // d(beta)
                    jacobian[kaMatrix.length + row] += -1;

                    // d(nonSpecKa)
                    jacobian[(2*kaMatrix.length) + row] += (-1 * alpha * (concGammaPos + concGammaNeg)) / bigProdPlus1Sq;

                    // d(posGamma)
                    jacobian[ ((numPreVars+i) * kaMatrix.length) + row] = (-1 * alpha * (kaPos+nonSpecKa)) / bigProdPlus1Sq;

                    // d(negGamma )
                    jacobian[ (kaMatrix.length*(numPreVars + this.globalNumWindows)) + ((numPreVars+i) * kaMatrix.length) + row] = (-1 * alpha * (kaNeg+nonSpecKa)) / bigProdPlus1Sq;
                }
            }
        }

        else { // isSelfRevComp == 1
            for (int row=0; row < kaMatrix.length; row++) {
//                 if ((row % 5000) == 0) {
//                     System.out.println("kaMatrix.length = "+kaMatrix.length+" : current row = "+row);
//                 }

                for (int i=0; i < this.globalNumWindows; i++) {
                    double concGammaPos = parList[numPreVars + i];
                    double kaPos = kaMatrix[row][i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd = (concGammaPos*kaPos) + (concGammaPos*kaPos) + ((concGammaPos+concGammaPos)*nonSpecKa);
                    double bigProdPlus1Sq = Math.pow((bigProd +1),2);

                    // d(alpha)
                    jacobian[row] += -1 * (bigProd / (bigProd +1));

                    // d(beta)
                    jacobian[kaMatrix.length + row] += -1;

                    // d(nonSpecKa)
                    jacobian[(2*kaMatrix.length) + row] += (-1 * alpha * (concGammaPos + concGammaPos)) / bigProdPlus1Sq ;

                    // d(posGamma)
                    jacobian[ ((numPreVars+i) * kaMatrix.length) + row] = (-1 * alpha * ((2*kaPos) + (2*nonSpecKa))) / bigProdPlus1Sq;

                    // d(negGamma )
                    // jacobian[ (kaMatrix.length*(numPreVars + this.globalNumWindows)) + ((numPreVars+i) * kaMatrix.length) + row] = 0;
                }
            }
        }
        //System.out.println("jacobian[0] = "+StringTools.toString(Arrays.copyOfRange(jacobian, 0, numPreVars+this.globalNumWindows), "\t"));
        return(jacobian);
    }

    public double[] getResidualGammas(double[] parList) {
        return(MathTools.subtract(this.globalY, getPredictedGammas(parList)));
    }

    public double getMinimizationGammas(double[] parList) {
        return(MathTools.L2Distance(getResidualGammas(parList)));
    }

    //////////////////////////////////////////////////////////////////////////
    // if isSelfRevComp then look only at the positive strand
    //////////////////////////////////////////////////////////////////////////
    public double[] getPredictedGammas(double[] parList) {
        //System.out.println("parList = "+StringTools.toString(parList, "\t"));

        int numPreVars = 3;
        double[][] kaMatrix = this.globalX[0];
        double[] predicted = new double[this.globalY.length];
        double alpha = parList[0];
        double beta = parList[1];
        double nonSpecKa = parList[2];

        if (!this.globalIsSelfRevComp) {
            for (int row=0; row < this.globalY.length; row++ ) {
                for (int i=0; i < this.globalNumWindows; i++ ) {
                    double concGammaPos = parList[numPreVars + i];
                    double kaPos = kaMatrix[row][i];

                    double concGammaNeg = parList[numPreVars + this.globalNumWindows + i];
                    double kaNeg = kaMatrix[row][this.globalNumWindows + i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd = (concGammaPos*kaPos) + (concGammaNeg*kaNeg) + ((concGammaPos+concGammaNeg)*nonSpecKa);

                    // predicted = (alpha * (bigProd / (bigProd + 1)) + beta
                    predicted[row] += (alpha * (bigProd / (bigProd + 1))) + beta;
                }
            }
        }
        else { // isSelfRevComp == 1
            for (int row=0; row < this.globalY.length; row++ ) {
                for (int i=0; i < this.globalNumWindows; i++ ) {
                        double concGammaPos = parList[numPreVars + i];
                        double kaPos = kaMatrix[row][i];

                        // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                        double bigProd = (2*(concGammaPos*kaPos)) + ((2*concGammaPos)*nonSpecKa);

                        // predicted = (alpha * (bigProd / (bigProd + 1)) + beta
                        predicted[row] += (alpha * (bigProd / (bigProd + 1))) + beta;
                    }
                }

        }
        return(predicted);
    }

    //////////////////////////////////////////////////////////////////////////
    // parList is in the format (alpha, beta, nonSpecKa, concGamma1pos, concGamma2pos, ...., concGamma1neg, ....)
    // kaMatrix is in the format (posStrandKa(1), posStrandKa(2), ...., negStrandKa(1), .....) for each row (probe)
    // lambda is the tuner for the L2 regularization ( - lambda[ (concGamma_I_pos - concGamma_I+1_pos)^2 ] )
    // alpha is an overall signal intensity scaler (models the light source, reflection intensity, light meter sensitivity)
    // beta is an overall signal intensity constant (intercept)
    //
    // return value is really a matrix forced into an array
    // return array is in "banded jacobian" format = (d(e1,v1), d(e2,v1), ....., d(e1,v2), d(e2,v2), ......)
    //
    // if isSelfRevComp then look only at the positive strand
    //////////////////////////////////////////////////////////////////////////
    public double[] getJacobianGammaSqs(double[] parList) {
//         System.out.println("parList = "+StringTools.toString(parList, "\t"));
//         System.out.println("y.length = "+y.length);
//         System.out.println("x.length = "+x.length);
//         System.out.println("x[0].length = "+x[0].length);
//         System.out.println("numWindows = "+numWindows);
//         System.out.println("isSelfRevComp = "+Boolean.toString(isSelfRevComp));

        int numPreVars = 3;
        double[][] kaMatrix = this.globalX[0];
        int numParameters = numPreVars + (this.globalNumStrands*this.globalNumWindows);
        int numResiduals = kaMatrix.length;

        double[] jacobian = new double[numResiduals * numParameters];

        double alpha = parList[0];
        double beta = parList[1];
        double nonSpecKa = parList[2];

        double alphaSq = Math.pow(alpha,2);
        double betaSq = Math.pow(beta,2);
        double nonSpecKaSq = Math.pow(nonSpecKa,2);

        // HACK!!!!!!!
        // If the parList contains any negative values then return all zeros!
        if (MathTools.hasLessThan(parList, 0)) {
            return(jacobian);
        }
        else {
            globalLastParList = parList;
        }

        // If numStrands=2 and !selfRevComp
        if ((this.globalNumStrands==2) && (!this.globalIsSelfRevComp)) {
            for (int row=0; row < kaMatrix.length; row++) {
//                 if ((row % 5000) == 0) {
//                     System.out.println("kaMatrix.length = "+kaMatrix.length+" : current row = "+row);
//                 }

                for (int i=0; i < this.globalNumWindows; i++) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];
                    double concGammaNeg = parList[numPreVars + this.globalNumWindows + i];
                    double concGammaNegSq = Math.pow(concGammaNeg,2);
                    double kaNeg = kaMatrix[row][this.globalNumWindows + i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns);
                    double bigProd = (concGammaPosSq*kaPos) + (concGammaNegSq*kaNeg) + ((concGammaPosSq+concGammaNegSq)*nonSpecKaSq);
                    double bigProdPlus1Sq = Math.pow((bigProd +1),2);

                    // d(alpha)
                    jacobian[row] += -1 * (2*alpha) * (bigProd / (bigProd +1));

                    // d(beta)
                    jacobian[numResiduals + row] += -1 * (2*beta);

                    // d(nonSpecKa)
                    jacobian[(2*numResiduals) + row] += (-1 * alphaSq * (concGammaPosSq + concGammaNegSq) * (2*nonSpecKa)) / bigProdPlus1Sq;

                    // d(posGamma)
                    jacobian[ ((numPreVars+i) * numResiduals) + row] = (-1 * alphaSq * (kaPos+nonSpecKaSq) * (2*concGammaPos)) / bigProdPlus1Sq;

                    // d(negGamma )
                    jacobian[ (numResiduals*(numPreVars + this.globalNumWindows)) + ((numPreVars+i) * numResiduals) + row] = (-1 * alphaSq * (kaNeg+nonSpecKaSq) * (2*concGammaNeg)) / bigProdPlus1Sq;
                }
            }
        }

        // (numStrands=2 and selfRevComp) or numStrands=1
        else {
            for (int row=0; row < kaMatrix.length; row++) {
//                 if ((row % 5000) == 0) {
//                     System.out.println("kaMatrix.length = "+kaMatrix.length+" : current row = "+row);
//                 }

                for (int i=0; i < this.globalNumWindows; i++) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd = (2 *(concGammaPosSq*kaPos)) + ((2 * concGammaPosSq)*nonSpecKaSq);
                    double bigProdPlus1Sq = Math.pow((bigProd +1),2);

                    // d(alpha)
                    jacobian[row] += -1 * (2*alpha) * (bigProd / (bigProd +1));

                    // d(beta)
                    jacobian[numResiduals + row] += -1 * (2*beta);

                    // d(nonSpecKa)
                    jacobian[(2*numResiduals) + row] += (-1 * alphaSq * (2*concGammaPosSq) * (2*nonSpecKa)) / bigProdPlus1Sq;

                    // d(posGamma)
                    jacobian[ ((numPreVars+i) * numResiduals) + row] = (-1 * alphaSq * ((2*kaPos) + (2*nonSpecKaSq)) * (2*concGammaPos)) / bigProdPlus1Sq;

                    // d(negGamma )
                    // jacobian[ (numResiduals*(numPreVars + this.globalNumWindows)) + ((numPreVars+i) * numResiduals) + row] = 0;
                }
            }
        }
        return(jacobian);
    }

    public double[] getResidualGammaSqs(double[] parList) {
        return(MathTools.subtract(this.globalY, getPredictedGammaSqs(parList)));
    }

    public double getMinimizationGammaSqs(double[] parList) {
        return(MathTools.L2Distance(getResidualGammaSqs(parList)));
    }

    //////////////////////////////////////////////////////////////////////////
    // if isSelfRevComp then look only at the positive strand
    //////////////////////////////////////////////////////////////////////////
    public double[] getPredictedGammaSqs(double[] parList) {
//         System.out.println("parList = "+StringTools.toString(parList, "\t"));
//         System.out.println("y.length = "+y.length);
//         System.out.println("x.length = "+x.length);
//         System.out.println("x[0].length = "+x[0].length);
//         System.out.println("numWindows = "+numWindows);
//         System.out.println("isSelfRevComp = "+Boolean.toString(isSelfRevComp));

        int numPreVars = 3;
        double[][] kaMatrix = this.globalX[0];
        double[] predicted = new double[this.globalY.length];
        double alpha = parList[0];
        double beta = parList[1];
        double nonSpecKa = parList[2];

        double alphaSq = Math.pow(alpha,2);
        double betaSq = Math.pow(beta,2);
        double nonSpecKaSq = Math.pow(nonSpecKa,2);

        if (!this.globalIsSelfRevComp) {
            for (int row=0; row < this.globalY.length; row++ ) {
                for (int i=0; i < this.globalNumWindows; i++ ) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];

                    double concGammaNeg = parList[numPreVars + this.globalNumWindows + i];
                    double concGammaNegSq = Math.pow(concGammaNeg,2);
                    double kaNeg = kaMatrix[row][this.globalNumWindows + i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd = (concGammaPosSq*kaPos) + (concGammaNegSq*kaNeg) + ((concGammaPosSq+concGammaNegSq)*nonSpecKaSq);

                    // predicted = (alpha * (bigProd / (bigProd + 1)) + beta
                    predicted[row] += (alphaSq * (bigProd / (bigProd + 1))) + betaSq;
                }
            }
        }
        else { // isSelfRevComp == 1
            for (int row=0; row < this.globalY.length; row++ ) {
                for (int i=0; i < this.globalNumWindows; i++ ) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd = (2*(concGammaPosSq*kaPos)) + ((2*concGammaPosSq)*nonSpecKaSq);

                    // predicted = (alpha * (bigProd / (bigProd + 1)) + beta
                    predicted[row] += (alphaSq * (bigProd / (bigProd + 1))) + betaSq;
                }
            }

        }

//         System.out.println("L1(observed) = " + MathTools.getL1(this.globalY));
//         System.out.println("L1(predicted) = " + MathTools.getL1(predicted));
//         System.out.println("L2(observed - predicted) = " + MathTools.L2Distance(this.globalY, predicted));

        return(predicted);
    }

    //////////////////////////////////////////////////////////////////////////
    // parList is in the format (alpha, beta, nonSpecKa, concGamma1pos, concGamma2pos, ...., concGamma1neg, ....)
    // kaMatrix is in the format (posStrandKa(1), posStrandKa(2), ...., negStrandKa(1), .....) for each row (probe)
    // lambda is the tuner for the L2 regularization ( - lambda[ (concGamma_I_pos - concGamma_I+1_pos)^2 ] )
    // alpha is an overall signal intensity scaler (models the light source, reflection intensity, light meter sensitivity)
    // beta is an overall signal intensity constant (intercept)
    //
    // return value is really a matrix forced into an array
    // return array is in "banded jacobian" format = (d(e1,v1), d(e2,v1), ....., d(e1,v2), d(e2,v2), ......)
    //
    // if isSelfRevComp then look only at the positive strand
    //////////////////////////////////////////////////////////////////////////
    public double[] getJacobianGammaSqsRg(double[] parList) {
//         System.out.println("parList = "+StringTools.toString(parList, "\t"));
//         System.out.println("y.length = "+y.length);
//         System.out.println("x.length = "+x.length);
//         System.out.println("x[0].length = "+x[0].length);
//         System.out.println("numWindows = "+numWindows);
//         System.out.println("isSelfRevComp = "+Boolean.toString(isSelfRevComp));

        int numPreVars = 3;
        double[][] kaMatrix = this.globalX[0];
        int numResiduals = kaMatrix.length + 2*(this.globalNumWindows-1);
        int numParameters = numPreVars + (2*this.globalNumWindows);
        double lambdaSqrt = Math.pow(this.globalLambda, .5);

        double[] jacobian = new double[numResiduals * numParameters];
        double alpha = parList[0];
        double beta = parList[1];
        double nonSpecKa = parList[2];

        double alphaSq = Math.pow(alpha,2);
        double betaSq = Math.pow(beta,2);
        double nonSpecKaSq = Math.pow(nonSpecKa,2);

        if (!this.globalIsSelfRevComp) {
            for (int row=0; row < kaMatrix.length; row++) {
//                 if ((row % 5000) == 0) {
//                     System.out.println("kaMatrix.length = "+kaMatrix.length+" : current row = "+row);
//                 }

                for (int i=0; i < this.globalNumWindows; i++) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];
                    double concGammaNeg = parList[numPreVars + this.globalNumWindows + i];
                    double concGammaNegSq = Math.pow(concGammaNeg,2);
                    double kaNeg = kaMatrix[row][this.globalNumWindows + i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns);
                    double bigProd = (concGammaPosSq*kaPos) + (concGammaNegSq*kaNeg) + ((concGammaPosSq+concGammaNegSq)*nonSpecKaSq);
                    double bigProdPlus1Sq = Math.pow((bigProd +1),2);

                    // d(alpha)
                    jacobian[row] += -1 * (2*alpha) * (bigProd / (bigProd +1));

                    // d(beta)
                    jacobian[numResiduals + row] += -1 * (2*beta);

                    // d(nonSpecKa)
                    jacobian[(2*numResiduals) + row] += (-1 * alphaSq * (concGammaPosSq + concGammaNegSq)*(2*nonSpecKa)) / bigProdPlus1Sq;

                    // d(posGamma)
                    jacobian[ ((numPreVars+i) * numResiduals) + row] = (-1 * alphaSq * (kaPos+nonSpecKaSq) * (2*concGammaPos)) / bigProdPlus1Sq;

                    // d(negGamma )
                    jacobian[ ((numPreVars+this.globalNumWindows) * numResiduals) + ((numPreVars+i) * numResiduals) + row] = (-1 * alphaSq * (kaNeg+nonSpecKaSq) * (2*concGammaNeg)) / bigProdPlus1Sq;
                }
            }

            for (int j=0; j < this.globalNumWindows; j++) {
                double concGammaPos = parList[numPreVars + j];
                double concGammaNeg = parList[numPreVars + this.globalNumWindows + j];

                if (j == 0) {
                    jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + j] = lambdaSqrt*2*concGammaPos;

                    jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + (this.globalNumWindows-1) + j] = lambdaSqrt*2*concGammaNeg;
                }
                else if (j < this.globalNumWindows-1) {
                    jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + (j-1)] = -lambdaSqrt*2*concGammaPos;
                    jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + j] = lambdaSqrt*2*concGammaPos;

                    jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + (this.globalNumWindows-1) + (j-1)] = -lambdaSqrt*2*concGammaNeg;
                    jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + (this.globalNumWindows-1) + j] = lambdaSqrt*2*concGammaNeg;
                }
                else { // j == this.globalNumWindows-1
                    jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + (j-1)] = -lambdaSqrt*2*concGammaPos;

                    jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + (this.globalNumWindows-1) + (j-1)] = -lambdaSqrt*2*concGammaNeg;
                }
            }
        }

        else { // isSelfRevComp == 1
            for (int row=0; row < kaMatrix.length; row++) {
//                 if ((row % 5000) == 0) {
//                     System.out.println("kaMatrix.length = "+kaMatrix.length+" : current row = "+row);
//                 }

                for (int i=0; i < this.globalNumWindows; i++) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd = (2 *(concGammaPosSq*kaPos)) + ((2 * concGammaPosSq)*nonSpecKaSq);
                    double bigProdPlus1Sq = Math.pow((bigProd +1),2);

                    // d(alpha)
                    jacobian[row] += -1 * (2*alpha) * (bigProd / (bigProd +1));

                    // d(beta)
                    jacobian[numResiduals + row] += -1 * (2*beta);

                    // d(nonSpecKa)
                    jacobian[(2*numResiduals) + row] += (-1 * alphaSq * (2*concGammaPosSq) * (2*nonSpecKa)) / bigProdPlus1Sq;

                    // d(posGamma)
                    jacobian[ ((numPreVars+i) * numResiduals) + row] = (-1 * alphaSq * ((2*kaPos) + (2*nonSpecKaSq)) * (2*concGammaPos)) / bigProdPlus1Sq;

                    // d(negGamma )
                    // jacobian[ (numResiduals*(numPreVars + this.globalNumWindows)) + ((numPreVars+i) * numResiduals) + row] = 0;
                }
            }

            for (int j=0; j < this.globalNumWindows; j++) {
                double concGammaPos = parList[numPreVars + j];
                double concGammaNeg = parList[numPreVars + this.globalNumWindows + j];

                if (j == 0) {
                    jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + j] = lambdaSqrt*2*concGammaPos;

                    //jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + (this.globalNumWindows-1) + j] = lambdaSqrt*2*concGammaNeg;
                }
                else if (j < this.globalNumWindows-1) {
                    jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + (j-1)] = -lambdaSqrt*2*concGammaPos;
                    jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + j] = lambdaSqrt*2*concGammaPos;

                    //jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + (this.globalNumWindows-1) + (j-1)] = -lambdaSqrt*2*concGammaNeg;
                    //jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + (this.globalNumWindows-1) + j] = lambdaSqrt*2*concGammaNeg;
                }
                else { // j == this.globalNumWindows-1
                    jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + (j-1)] = -lambdaSqrt*2*concGammaPos;

                    //jacobian[((numPreVars+j) * numResiduals) + kaMatrix.length + (this.globalNumWindows-1) + (j-1)] = -lambdaSqrt*2*concGammaNeg;
                }
            }
        }
        return(jacobian);
    }

    public double[] getResidualGammaSqsRg(double[] parList) {
        int numPreVars = 3;
        double[] parSqList = MathTools.pow(parList, 2.0);

        double[] residuals = MathTools.subtract(this.globalY, getPredictedGammaSqsRg(parList));
        double[] posStrandRegs = getGammaRegularization(Arrays.copyOfRange(parSqList, numPreVars, numPreVars+this.globalNumWindows));
        double[] addedPosStrandRegs = ArrayUtils.addAll(residuals, posStrandRegs);
        double[] negStrandRegs;

        if (this.globalIsSelfRevComp) {
            negStrandRegs = new double[this.globalNumWindows-1];
        }
        else {
            negStrandRegs = getGammaRegularization(Arrays.copyOfRange(parSqList, numPreVars+this.globalNumWindows, parList.length));
        }

//         System.out.println("residuals.length = "+residuals.length);
//         System.out.println("posStrandRegs.length = "+posStrandRegs.length);
//         System.out.println("negStrandRegs.length = "+negStrandRegs.length);

        return(ArrayUtils.addAll(addedPosStrandRegs, negStrandRegs));
    }

    public double getMinimizationGammaSqsRg(double[] parList) {
        return(MathTools.L2Distance(getResidualGammaSqsRg(parList)));
    }

    //////////////////////////////////////////////////////////////////////////
    // if isSelfRevComp then look only at the positive strand
    //////////////////////////////////////////////////////////////////////////
    public double[] getPredictedGammaSqsRg(double[] parList) {
//         System.out.println("parList = "+StringTools.toString(parList, "\t"));
//         System.out.println("y.length = "+y.length);
//         System.out.println("x.length = "+x.length);
//         System.out.println("x[0].length = "+x[0].length);
//         System.out.println("numWindows = "+numWindows);
//         System.out.println("isSelfRevComp = "+Boolean.toString(isSelfRevComp));

        int numPreVars = 3;
        double[][] kaMatrix = this.globalX[0];
        double[] predicted = new double[this.globalY.length];
        double alpha = parList[0];
        double beta = parList[1];
        double nonSpecKa = parList[2];

        double alphaSq = Math.pow(alpha,2);
        double betaSq = Math.pow(beta,2);
        double nonSpecKaSq = Math.pow(nonSpecKa,2);

        if (!this.globalIsSelfRevComp) {
            for (int row=0; row < this.globalY.length; row++ ) {
                for (int i=0; i < this.globalNumWindows; i++ ) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];

                    double concGammaNeg = parList[numPreVars + this.globalNumWindows + i];
                    double concGammaNegSq = Math.pow(concGammaNeg,2);
                    double kaNeg = kaMatrix[row][this.globalNumWindows + i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd = (concGammaPosSq*kaPos) + (concGammaNegSq*kaNeg) + ((concGammaPosSq+concGammaNegSq)*nonSpecKaSq);

                    // predicted = (alpha * (bigProd / (bigProd + 1)) + beta
                    predicted[row] += (alphaSq * (bigProd / (bigProd + 1))) + betaSq;
                }
            }
        }
        else { // isSelfRevComp == 1
            for (int row=0; row < this.globalY.length; row++ ) {
                for (int i=0; i < this.globalNumWindows; i++ ) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd = (2*(concGammaPosSq*kaPos)) + ((2*concGammaPosSq)*nonSpecKaSq);

                    // predicted = (alpha * (bigProd / (bigProd + 1)) + beta
                    predicted[row] += (alphaSq * (bigProd / (bigProd + 1))) + betaSq;
                }
            }

        }

//         System.out.println("L1(observed) = " + MathTools.getL1(this.globalY));
//         System.out.println("L1(predicted) = " + MathTools.getL1(predicted));
//         System.out.println("L2(observed - predicted) = " + MathTools.L2Distance(this.globalY, predicted));

        return(predicted);
    }

    public double[] getGammaRegularization(double[] gammas) {
        double[] residuals = new double[gammas.length - 1];
        double lambdaSqrt = Math.pow(this.globalLambda, .5);
        for (int i=0; i < gammas.length - 1; i++) {
            residuals[i] = lambdaSqrt * (gammas[i] - gammas[i+1]);
        }
        return(residuals);
    }

    //////////////////////////////////////////////////////////////////////////
    // parList is in the format (alpha, beta, nonSpecKa, concGamma1pos, concGamma2pos, ...., concGamma1neg, ....)
    // kaMatrix is in the format (posStrandKa(1), posStrandKa(2), ...., negStrandKa(1), .....) for each row (probe)
    // lambda is the tuner for the L2 regularization ( - lambda[ (concGamma_I_pos - concGamma_I+1_pos)^2 ] )
    // alpha is an overall signal intensity scaler (models the light source, reflection intensity, light meter sensitivity)
    // beta is an overall signal intensity constant (intercept)
    //
    // return value is really a matrix forced into an array
    // return array is in "banded jacobian" format = (d(e1,v1), d(e2,v1), ....., d(e1,v2), d(e2,v2), ......)
    //
    // if isSelfRevComp then look only at the positive strand
    //////////////////////////////////////////////////////////////////////////
    public double[] getJacobianGammaSqsL1(double[] parList) {
//         System.out.println("parList = "+StringTools.toString(parList, "\t"));
//         System.out.println("y.length = "+y.length);
//         System.out.println("x.length = "+x.length);
//         System.out.println("x[0].length = "+x[0].length);
//         System.out.println("numWindows = "+numWindows);
//         System.out.println("isSelfRevComp = "+Boolean.toString(isSelfRevComp));

        int numPreVars = 3;
        double[][] kaMatrix = this.globalX[0];
        int numParameters = numPreVars + (this.globalNumStrands*this.globalNumWindows);
        int numResiduals = kaMatrix.length + numParameters;
        double lambdaSqrt = Math.pow(this.globalLambda, .5);

        double[] jacobian = new double[numResiduals * numParameters];

        double alpha = parList[0];
        double beta = parList[1];
        double nonSpecKa = parList[2];

        double alphaSq = Math.pow(alpha,2);
        double betaSq = Math.pow(beta,2);
        double nonSpecKaSq = Math.pow(nonSpecKa,2);

        // HACK!!!!!!!
        // If the parList contains any negative values then return all zeros!
        if (MathTools.hasLessThan(parList, 0)) {
            return(jacobian);
        }
        else {
            globalLastParList = parList;
        }

        // If numStrands=2 and !selfRevComp
        if ((this.globalNumStrands==2) && (!this.globalIsSelfRevComp)) {
            for (int row=0; row < kaMatrix.length; row++) {
//                 if ((row % 5000) == 0) {
//                     System.out.println("kaMatrix.length = "+kaMatrix.length+" : current row = "+row);
//                 }

                for (int i=0; i < this.globalNumWindows; i++) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];
                    double concGammaNeg = parList[numPreVars + this.globalNumWindows + i];
                    double concGammaNegSq = Math.pow(concGammaNeg,2);
                    double kaNeg = kaMatrix[row][this.globalNumWindows + i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns);
                    double bigProd = (concGammaPosSq*kaPos) + (concGammaNegSq*kaNeg) + ((concGammaPosSq+concGammaNegSq)*nonSpecKaSq);
                    double bigProdPlus1Sq = Math.pow((bigProd +1),2);

                    // d(alpha)
                    jacobian[row] += -1 * (2*alpha) * (bigProd / (bigProd +1));

                    // d(beta)
                    jacobian[numResiduals + row] += -1 * (2*beta);

                    // d(nonSpecKa)
                    jacobian[(2*numResiduals) + row] += (-1 * alphaSq * (concGammaPosSq + concGammaNegSq) * (2*nonSpecKa)) / bigProdPlus1Sq;

                    // d(posGamma)
                    jacobian[ ((numPreVars+i) * numResiduals) + row] = (-1 * alphaSq * (kaPos+nonSpecKaSq) * (2*concGammaPos)) / bigProdPlus1Sq;

                    // d(negGamma )
                    jacobian[ ((numPreVars+this.globalNumWindows) * numResiduals) + ((numPreVars+i) * numResiduals) + row] = (-1 * alphaSq * (kaNeg+nonSpecKaSq) * (2*concGammaNeg)) / bigProdPlus1Sq;
                }
            }

            for (int j=0; j < parList.length; j++) {
                jacobian[(j * numResiduals) + kaMatrix.length + j] = lambdaSqrt*2*parList[j];
            }
        }

        // (numStrands=2 and selfRevComp) or numStrands=1
        else {
            for (int row=0; row < kaMatrix.length; row++) {
//                 if ((row % 5000) == 0) {
//                     System.out.println("kaMatrix.length = "+kaMatrix.length+" : current row = "+row);
//                 }

                for (int i=0; i < this.globalNumWindows; i++) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd = (2 *(concGammaPosSq*kaPos)) + ((2 * concGammaPosSq)*nonSpecKaSq);
                    double bigProdPlus1Sq = Math.pow((bigProd +1),2);

                    // d(alpha)
                    jacobian[row] += -1 * (2*alpha) * (bigProd / (bigProd +1));

                    // d(beta)
                    jacobian[numResiduals + row] += -1 * (2*beta);

                    // d(nonSpecKa)
                    jacobian[(2*numResiduals) + row] += (-1 * alphaSq * (2*concGammaPosSq) * (2*nonSpecKa)) / bigProdPlus1Sq;

                    // d(posGamma)
                    jacobian[ ((numPreVars+i) * numResiduals) + row] = (-1 * alphaSq * ((2*kaPos) + (2*nonSpecKaSq)) * (2*concGammaPos)) / bigProdPlus1Sq;

                    // d(negGamma )
                    // jacobian[ (numResiduals*(numPreVars + this.globalNumWindows)) + ((numPreVars+i) * numResiduals) + row] = 0;
                }
            }

            // do not loop over gammaNegatives
            for (int j=0; j < numPreVars + (1*this.globalNumWindows); j++) {
                jacobian[(j * numResiduals) + kaMatrix.length + j] = lambdaSqrt*2*parList[j];
            }
        }

        return(jacobian);
    }

    public double[] getResidualGammaSqsL1(double[] parList) {

        // HACK!!!!!!!
        // If the parList contains any negative values then return all zeros!
        if (MathTools.hasLessThan(parList, 0)) {
            return(new double[parList.length]);
        }
        else {
            globalLastParList = parList;
        }

        int numPreVars = 3;
        double[] parSqList = MathTools.pow(parList, 2.0);

        double[] residuals = MathTools.subtract(this.globalY, getPredictedGammaSqsL1(parList));
        double[] posStrandRegs = getGammaL1Regularization(Arrays.copyOfRange(parSqList, 0, numPreVars+this.globalNumWindows));
        residuals = ArrayUtils.addAll(residuals, posStrandRegs);

        if ((this.globalNumStrands==2) && (!this.globalIsSelfRevComp)) {
            double[] negStrandRegs = getGammaL1Regularization(Arrays.copyOfRange(parSqList, numPreVars+this.globalNumWindows, parList.length));
            residuals = ArrayUtils.addAll(residuals, negStrandRegs);
        }

//         System.out.println("residuals.length = "+residuals.length);
//         System.out.println("posStrandRegs.length = "+posStrandRegs.length);
//         System.out.println("negStrandRegs.length = "+negStrandRegs.length);

        return(residuals);
    }

    public double getMinimizationGammaSqsL1(double[] parList) {
        return(MathTools.L2Distance(getResidualGammaSqsL1(parList)));
    }

    //////////////////////////////////////////////////////////////////////////
    // if isSelfRevComp then look only at the positive strand
    //////////////////////////////////////////////////////////////////////////
    public double[] getPredictedGammaSqsL1(double[] parList) {
//         System.out.println("parList = "+StringTools.toString(parList, "\t"));
//         System.out.println("y.length = "+y.length);
//         System.out.println("x.length = "+x.length);
//         System.out.println("x[0].length = "+x[0].length);
//         System.out.println("numWindows = "+numWindows);
//         System.out.println("isSelfRevComp = "+Boolean.toString(isSelfRevComp));

        int numPreVars = 3;
        double[][] kaMatrix = this.globalX[0];
        double[] predicted = new double[this.globalY.length];
        double alpha = parList[0];
        double beta = parList[1];
        double nonSpecKa = parList[2];

        double alphaSq = Math.pow(alpha,2);
        double betaSq = Math.pow(beta,2);
        double nonSpecKaSq = Math.pow(nonSpecKa,2);

        if ((this.globalNumStrands==2) && (!this.globalIsSelfRevComp)) {
            for (int row=0; row < this.globalY.length; row++ ) {
                for (int i=0; i < this.globalNumWindows; i++ ) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];

                    double concGammaNeg = parList[numPreVars + this.globalNumWindows + i];
                    double concGammaNegSq = Math.pow(concGammaNeg,2);
                    double kaNeg = kaMatrix[row][this.globalNumWindows + i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd = (concGammaPosSq*kaPos) + (concGammaNegSq*kaNeg) + ((concGammaPosSq+concGammaNegSq)*nonSpecKaSq);

                    // predicted = (alpha * (bigProd / (bigProd + 1)) + beta
                    predicted[row] += (alphaSq * (bigProd / (bigProd + 1))) + betaSq;
                }
            }
        }
        else { // isSelfRevComp == 1
            for (int row=0; row < this.globalY.length; row++ ) {
                for (int i=0; i < this.globalNumWindows; i++ ) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd = (2*(concGammaPosSq*kaPos)) + ((2*concGammaPosSq)*nonSpecKaSq);

                    // predicted = (alpha * (bigProd / (bigProd + 1)) + beta
                    predicted[row] += (alphaSq * (bigProd / (bigProd + 1))) + betaSq;
                }
            }

        }

//         System.out.println("L1(observed) = " + MathTools.getL1(this.globalY));
//         System.out.println("L1(predicted) = " + MathTools.getL1(predicted));
//         System.out.println("L2(observed - predicted) = " + MathTools.L2Distance(this.globalY, predicted));

        return(predicted);
    }

    public double[] getGammaL1Regularization(double[] gammas) {

        double lambdaSqrt = Math.pow(this.globalLambda, .5);
        double[] gammasSqrt = MathTools.pow(gammas, .5);
        double[] l1TermsArray = MathTools.multiply(lambdaSqrt, gammasSqrt);
        return(l1TermsArray);
        //return(new double[gammas.length]);
    }

    //////////////////////////////////////////////////////////////////////////
    // parList is in the format (alpha, beta, nonSpecKa, concGamma1pos, concGamma2pos, ...., concGamma1neg, ....)
    // kaMatrix is in the format (posStrandKa(1), posStrandKa(2), ...., negStrandKa(1), .....) for each row (probe)
    // lambda is the tuner for the L2 regularization ( - lambda[ (concGamma_I_pos - concGamma_I+1_pos)^2 ] )
    // alpha is an overall signal intensity scaler (models the light source, reflection intensity, light meter sensitivity)
    // beta is an overall signal intensity constant (intercept)
    //
    // return value is really a matrix forced into an array
    // return array is in "banded jacobian" format = (d(e1,v1), d(e2,v1), ....., d(e1,v2), d(e2,v2), ......)
    //
    // if isSelfRevComp then look only at the positive strand
    //////////////////////////////////////////////////////////////////////////
    public double[] getJacobianGammaSqsL2(double[] parList) {
//         System.out.println("parList = "+StringTools.toString(parList, "\t"));
//         System.out.println("y.length = "+y.length);
//         System.out.println("x.length = "+x.length);
//         System.out.println("x[0].length = "+x[0].length);
//         System.out.println("numWindows = "+numWindows);
//         System.out.println("isSelfRevComp = "+Boolean.toString(isSelfRevComp));

        int numPreVars = 3;
        double[][] kaMatrix = this.globalX[0];
        int numParameters = numPreVars + (2*this.globalNumWindows);
        int numResiduals = kaMatrix.length + numParameters;
        double lambdaSqrt = Math.pow(this.globalLambda, .5);

        double[] jacobian = new double[numResiduals * numParameters];
        double alpha = parList[0];
        double beta = parList[1];
        double nonSpecKa = parList[2];

        double alphaSq = Math.pow(alpha,2);
        double betaSq = Math.pow(beta,2);
        double nonSpecKaSq = Math.pow(nonSpecKa,2);

        if (!this.globalIsSelfRevComp) {
            for (int row=0; row < kaMatrix.length; row++) {
//                 if ((row % 5000) == 0) {
//                     System.out.println("kaMatrix.length = "+kaMatrix.length+" : current row = "+row);
//                 }

                for (int i=0; i < this.globalNumWindows; i++) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];
                    double concGammaNeg = parList[numPreVars + this.globalNumWindows + i];
                    double concGammaNegSq = Math.pow(concGammaNeg,2);
                    double kaNeg = kaMatrix[row][this.globalNumWindows + i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns);
                    double bigProd = (concGammaPosSq*kaPos) + (concGammaNegSq*kaNeg) + ((concGammaPosSq+concGammaNegSq)*nonSpecKaSq);
                    double bigProdPlus1Sq = Math.pow((bigProd +1),2);

                    // d(alpha)
                    jacobian[row] += -1 * (2*alpha) * (bigProd / (bigProd +1));

                    // d(beta)
                    jacobian[numResiduals + row] += -1 * (2*beta);

                    // d(nonSpecKa)
                    jacobian[(2*numResiduals) + row] += (-1 * alphaSq * (concGammaPosSq + concGammaNegSq) * (2*nonSpecKa)) / bigProdPlus1Sq;

                    // d(posGamma)
                    jacobian[ ((numPreVars+i) * numResiduals) + row] = (-1 * alphaSq * (kaPos+nonSpecKaSq) * (2*concGammaPos)) / bigProdPlus1Sq;

                    // d(negGamma )
                    jacobian[ ((numPreVars+this.globalNumWindows) * numResiduals) + ((numPreVars+i) * numResiduals) + row] = (-1 * alphaSq * (kaNeg+nonSpecKaSq) * (2*concGammaNeg)) / bigProdPlus1Sq;
                }
            }

            for (int j=0; j < parList.length; j++) {
                jacobian[(j * numResiduals) + kaMatrix.length + j] = lambdaSqrt*2*parList[j];
            }
        }

        else { // isSelfRevComp == 1
            for (int row=0; row < kaMatrix.length; row++) {
//                 if ((row % 5000) == 0) {
//                     System.out.println("kaMatrix.length = "+kaMatrix.length+" : current row = "+row);
//                 }

                for (int i=0; i < this.globalNumWindows; i++) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd = (2 *(concGammaPosSq*kaPos)) + ((2 * concGammaPosSq)*nonSpecKaSq);
                    double bigProdPlus1Sq = Math.pow((bigProd +1),2);

                    // d(alpha)
                    jacobian[row] += -1 * (2*alpha) * (bigProd / (bigProd +1));

                    // d(beta)
                    jacobian[numResiduals + row] += -1 * (2*beta);

                    // d(nonSpecKa)
                    jacobian[(2*numResiduals) + row] += (-1 * alphaSq * (2*concGammaPosSq) * (2*nonSpecKa)) / bigProdPlus1Sq;

                    // d(posGamma)
                    jacobian[ ((numPreVars+i) * numResiduals) + row] = (-1 * alphaSq * ((2*kaPos) + (2*nonSpecKaSq)) * (2*concGammaPos)) / bigProdPlus1Sq;

                    // d(negGamma )
                    // jacobian[ (numResiduals*(numPreVars + this.globalNumWindows)) + ((numPreVars+i) * numResiduals) + row] = 0;
                }
            }

            // do not loop over gammaNegatives
            for (int j=0; j < numPreVars + (1*this.globalNumWindows); j++) {
                jacobian[(j * numResiduals) + kaMatrix.length + j] = lambdaSqrt*2*parList[j];
            }
        }
        return(jacobian);
    }

    public double[] getResidualGammaSqsL2(double[] parList) {
        int numPreVars = 3;
        double[] parSqList = MathTools.pow(parList, 2.0);

        double[] residuals = MathTools.subtract(this.globalY, getPredictedGammaSqsL2(parList));
        double[] posStrandRegs = getGammaL2Regularization(Arrays.copyOfRange(parSqList, 0, numPreVars+this.globalNumWindows));
        double[] addedPosStrandRegs = ArrayUtils.addAll(residuals, posStrandRegs);
        double[] negStrandRegs;

        if (this.globalIsSelfRevComp) {
            negStrandRegs = new double[this.globalNumWindows];
        }
        else {
            negStrandRegs = getGammaL2Regularization(Arrays.copyOfRange(parSqList, numPreVars+this.globalNumWindows, parList.length));
        }

//         System.out.println("residuals.length = "+residuals.length);
//         System.out.println("posStrandRegs.length = "+posStrandRegs.length);
//         System.out.println("negStrandRegs.length = "+negStrandRegs.length);

        return(ArrayUtils.addAll(addedPosStrandRegs, negStrandRegs));
    }

    public double getMinimizationGammaSqsL2(double[] parList) {
        return(MathTools.L2Distance(getResidualGammaSqsL2(parList)));
    }

    //////////////////////////////////////////////////////////////////////////
    // if isSelfRevComp then look only at the positive strand
    //////////////////////////////////////////////////////////////////////////
    public double[] getPredictedGammaSqsL2(double[] parList) {
//         System.out.println("parList = "+StringTools.toString(parList, "\t"));
//         System.out.println("y.length = "+y.length);
//         System.out.println("x.length = "+x.length);
//         System.out.println("x[0].length = "+x[0].length);
//         System.out.println("numWindows = "+numWindows);
//         System.out.println("isSelfRevComp = "+Boolean.toString(isSelfRevComp));

        int numPreVars = 3;
        double[][] kaMatrix = this.globalX[0];
        double[] predicted = new double[this.globalY.length];
        double alpha = parList[0];
        double beta = parList[1];
        double nonSpecKa = parList[2];

        double alphaSq = Math.pow(alpha,2);
        double betaSq = Math.pow(beta,2);
        double nonSpecKaSq = Math.pow(nonSpecKa,2);

        if (!this.globalIsSelfRevComp) {
            for (int row=0; row < this.globalY.length; row++ ) {
                for (int i=0; i < this.globalNumWindows; i++ ) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];

                    double concGammaNeg = parList[numPreVars + this.globalNumWindows + i];
                    double concGammaNegSq = Math.pow(concGammaNeg,2);
                    double kaNeg = kaMatrix[row][this.globalNumWindows + i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd = (concGammaPosSq*kaPos) + (concGammaNegSq*kaNeg) + ((concGammaPosSq+concGammaNegSq)*nonSpecKaSq);

                    // predicted = (alpha * (bigProd / (bigProd + 1)) + beta
                    predicted[row] += (alphaSq * (bigProd / (bigProd + 1))) + betaSq;
                }
            }
        }
        else { // isSelfRevComp == 1
            for (int row=0; row < this.globalY.length; row++ ) {
                for (int i=0; i < this.globalNumWindows; i++ ) {
                    double concGammaPos = parList[numPreVars + i];
                    double concGammaPosSq = Math.pow(concGammaPos,2);
                    double kaPos = kaMatrix[row][i];

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd = (2*(concGammaPosSq*kaPos)) + ((2*concGammaPosSq)*nonSpecKaSq);

                    // predicted = (alpha * (bigProd / (bigProd + 1)) + beta
                    predicted[row] += (alphaSq * (bigProd / (bigProd + 1))) + betaSq;
                }
            }

        }

//         System.out.println("L1(observed) = " + MathTools.getL1(this.globalY));
//         System.out.println("L1(predicted) = " + MathTools.getL1(predicted));
//         System.out.println("L2(observed - predicted) = " + MathTools.L2Distance(this.globalY, predicted));

        return(predicted);
    }

    public double[] getGammaL2Regularization(double[] gammas) {
        double [] l1TermsArray = new double[gammas.length];

        double lambdaSqrt = Math.pow(this.globalLambda, .5);
        //double[] gammasSqrt = MathTools.pow(gammas, .5);

        for (int i=0; i < gammas.length; i++) {
            l1TermsArray[i] = lambdaSqrt * gammas[i];
        }
        return(l1TermsArray);
    }

    public double[] getJacobianKas(double[] parList) {
        return(new double[0]);
    }

    public double[] getResidualKas(double[] parList) {
        return(new double[0]);
    }

    public double[] getPredictedKas(double[] parList) {
        return(new double[0]);
    }

    //////////////////////////////////////////////////////////////////////////
    // parList is in the format (Ka(A), Ka(C), Ka(G), Ks(T))
    // kaMatrix is in the format (posStrandKa(1), posStrandKa(2), ...., negStrandKa(1), .....) for each row (probe)
    // alpha is an overall signal intensity scaler (models the light source, reflection intensity, light meter sensitivity)
    // beta is an overall signal intensity constant (intercept)
    //
    // return value is really a matrix forced into an array
    // return array is in "banded jacobian" format = (d(e1,v1), d(e2,v1), ....., d(e1,v2), d(e2,v2), ......)
    //
    // if isSelfRevComp then look only at the positive strand
    //////////////////////////////////////////////////////////////////////////
    public double[] getJacobianKaSqs(double[] parList) {
//         System.out.println("parList = "+StringTools.toString(parList, "\t"));
//         System.out.println("y.length = "+y.length);
//         System.out.println("x.length = "+x.length);
//         System.out.println("x[0].length = "+x[0].length);
//         System.out.println("numWindows = "+numWindows);
//         System.out.println("isSelfRevComp = "+Boolean.toString(isSelfRevComp));

        int numPreVars = 0;
        int numParameters = 4;
        int numResiduals = this.globalY.length;

        double[] jacobian = new double[numResiduals * numParameters];

        double alpha = this.globalAlpha;
        double beta = this.globalBeta;
        double nonSpecKa = this.globalNonSpecKa;

        // HACK!!!!!!!
        // If the parList contains any negative values then return all zeros!
        if (MathTools.hasLessThan(parList, 0)) {
            return(jacobian);
        }
        else {
            globalLastParList = parList;
        }

        // If numStrands=2 and !selfRevComp
        if ((this.globalNumStrands==2) && (!this.globalIsSelfRevComp)) {
            for (int row=0; row < this.globalX[0].length; row++) {
//                 if ((row % 5000) == 0) {
//                     System.out.println("kaMatrix.length = "+kaMatrix.length+" : current row = "+row);
//                 }

                for (int i=0; i < this.globalNumWindows; i++) {
                    double concGammaPos = this.globalConcGammas[i];
                    double concGammaNeg = this.globalConcGammas[this.globalNumWindows + i];
                    double bigSum = 0;

                    for (int nucl=0; nucl < this.globalX.length; nucl++ ) {
                        double kaPos = this.globalX[nucl][row][i];
                        double kaNeg = this.globalX[nucl][row][this.globalNumWindows + i];
                        double kaNucSq = Math.pow(parList[nucl],2);
                        bigSum += kaNucSq * ((concGammaPos*kaPos) + (concGammaNeg*kaNeg));
                    }

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd =  bigSum + ((concGammaPos+concGammaNeg)*nonSpecKa);
                    double bigProdPlus1Sq = Math.pow((bigProd +1),2);

                    for (int nucl=0; nucl < this.globalX.length; nucl++ ) {
                        double kaPos = this.globalX[nucl][row][i];
                        double kaNeg = this.globalX[nucl][row][this.globalNumWindows + i];

                        // d(Ka(nucl))
                        jacobian[(nucl*numResiduals) + row] += (-1 * alpha * (2*parList[nucl]) * ((concGammaPos*kaPos) + (concGammaNeg*kaNeg))) / bigProdPlus1Sq;
                    }
                }
            }
        }

        // (numStrands=2 and selfRevComp) or numStrands=1
        else {
            for (int row=0; row < this.globalX[0].length; row++) {
//                 if ((row % 5000) == 0) {
//                     System.out.println("kaMatrix.length = "+kaMatrix.length+" : current row = "+row);
//                 }

                for (int i=0; i < this.globalNumWindows; i++) {
                    double concGammaPos = this.globalConcGammas[i];
                    double bigSum = 0;

                    for (int nucl=0; nucl < this.globalX.length; nucl++ ) {
                        double kaPos = this.globalX[nucl][row][i];
                        double kaNucSq = Math.pow(parList[nucl],2);
                        bigSum += kaNucSq * (2*(concGammaPos*kaPos));
                    }

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd =  bigSum + ((2*concGammaPos)*nonSpecKa);
                    double bigProdPlus1Sq = Math.pow((bigProd +1),2);

                    for (int nucl=0; nucl < this.globalX.length; nucl++ ) {
                        double kaPos = this.globalX[nucl][row][i];

                        // d(Ka(nucl))
                        jacobian[(nucl*numResiduals) + row] += (-1 * alpha * (2*parList[nucl]) * ((2*concGammaPos)*kaPos)) / bigProdPlus1Sq;
                    }
                }
            }
        }

        return(jacobian);
    }

    public double[] getResidualKaSqs(double[] parList) {

        // HACK!!!!!!!
        // If the parList contains any negative values then return all zeros!
        if (MathTools.hasLessThan(parList, 0)) {
            return(new double[parList.length]);
        }
        else {
            globalLastParList = parList;
        }

        double[] residuals = MathTools.subtract(this.globalY, getPredictedKaSqs(parList));
        return(residuals);

//         int numPreVars = 3;
//         double[] parSqList = MathTools.pow(parList, 2.0);

//         double[] posStrandRegs = getKaRegularization(Arrays.copyOfRange(parSqList, 0, numPreVars+this.globalNumWindows));
//         residuals = ArrayUtils.addAll(residuals, posStrandRegs);

//         if ((this.globalNumStrands==2) && (!this.globalIsSelfRevComp)) {
//             double[] negStrandRegs = getKaRegularization(Arrays.copyOfRange(parSqList, numPreVars+this.globalNumWindows, parList.length));
//             residuals = ArrayUtils.addAll(residuals, negStrandRegs);
//         }

// //         System.out.println("residuals.length = "+residuals.length);
// //         System.out.println("posStrandRegs.length = "+posStrandRegs.length);
// //         System.out.println("negStrandRegs.length = "+negStrandRegs.length);

//         return(residuals);
    }

    public double getMinimizationKaSqs(double[] parList) {
        return(MathTools.L2Distance(getResidualKaSqs(parList)));
    }

    //////////////////////////////////////////////////////////////////////////
    // if isSelfRevComp then look only at the positive strand
    //////////////////////////////////////////////////////////////////////////
    public double[] getPredictedKaSqs(double[] parList) {
//         System.out.println("parList = "+StringTools.toString(parList, "\t"));
//         System.out.println("y.length = "+y.length);
//         System.out.println("x.length = "+x.length);
//         System.out.println("x[0].length = "+x[0].length);
//         System.out.println("numWindows = "+numWindows);
//         System.out.println("isSelfRevComp = "+Boolean.toString(isSelfRevComp));

        int numPreVars = 0;
        double[] predicted = new double[this.globalY.length];
//         double alpha = parList[0];
//         double beta = parList[1];
//         double nonSpecKa = parList[2];
        double alpha = this.globalAlpha;
        double beta = this.globalBeta;
        double nonSpecKa = this.globalNonSpecKa;

        if ((this.globalNumStrands==2) && (!this.globalIsSelfRevComp)) {
            for (int row=0; row < this.globalY.length; row++ ) {
                for (int i=0; i < this.globalNumWindows; i++ ) {
                    double concGammaPos = this.globalConcGammas[i];
                    double concGammaNeg = this.globalConcGammas[this.globalNumWindows + i];
                    double bigSum = 0;

                    for (int nucl=0; nucl < this.globalX.length; nucl++ ) {
                        double kaPos = this.globalX[nucl][row][i];
                        double kaNeg = this.globalX[nucl][row][this.globalNumWindows + i];
                        double kaNucSq = Math.pow(parList[nucl],2);
                        bigSum += kaNucSq * ((concGammaPos*kaPos) + (concGammaNeg*kaNeg));
                    }

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd =  bigSum + ((concGammaPos+concGammaNeg)*nonSpecKa);

                    // predicted = (alpha * (bigProd / (bigProd + 1)) + beta
                    predicted[row] += (alpha * (bigProd / (bigProd + 1))) + beta;
                }
            }
        }
        else { // isSelfRevComp == 1
            for (int row=0; row < this.globalY.length; row++ ) {
                for (int i=0; i < this.globalNumWindows; i++ ) {
                    double concGammaPos = this.globalConcGammas[i];
                    double bigSum = 0;

                    for (int nucl=0; nucl < this.globalX.length; nucl++ ) {
                        double kaPos = this.globalX[nucl][row][i];
                        double kaNucSq = Math.pow(parList[nucl],2);
                        bigSum += kaNucSq * (2*(concGammaPos*kaPos));
                    }

                    // bigProd = (X_i,pos * Ka_i,pos) + (X_i,neg * Ka_i,neg) + ((X_i,pos + X_i,neg) * Ka_ns)
                    double bigProd =  bigSum + ((2*concGammaPos)*nonSpecKa);

                    // predicted = (alpha * (bigProd / (bigProd + 1)) + beta
                    predicted[row] += (alpha * (bigProd / (bigProd + 1))) + beta;
                }
            }

        }

//         System.out.println("L1(observed) = " + MathTools.getL1(this.globalY));
//         System.out.println("L1(predicted) = " + MathTools.getL1(predicted));
//         System.out.println("L2(observed - predicted) = " + MathTools.L2Distance(this.globalY, predicted));

        return(predicted);
    }

    public double[] getKaL1Regularization(double[] kas) {

        double lambdaSqrt = Math.pow(this.globalLambda, .5);
        double[] kasSqrt = MathTools.pow(kas, .5);
        double[] l1TermsArray = MathTools.multiply(lambdaSqrt, kasSqrt);
        return(l1TermsArray);
        //return(new double[kas.length]);
    }

    public static java.util.List<String> getRStatements(String file) {
        String evalsString = FileTools.readString(file, 0, null);
        //String evalsString = FileTools.readString(file, "#", true);

//         evalsString = evalsString.replace("}\n\n\n", "} evalsSeparator ");
//         evalsString = evalsString.replace("\n", " ");
//         evalsString = evalsString.replace("} evalsSeparator ", "}\n");

        String[] evalsStrings = evalsString.split("# New Eval");
        java.util.List<String> evalList = Arrays.asList(evalsStrings);
        return(evalList);
    }

    public static String symListCollectionDump(Collection<SymbolList> aSymListCollection) {
        StringBuffer output = new StringBuffer();
        try {
            for (SymbolList symList : aSymListCollection) {
                output.append(symList.seqString()+"\n");
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(output.toString());
    }


    public static void main(String[] args) throws IOException {
        FeatureReduce monster = new FeatureReduce();
        monster.parseCommandLine(args);
    }


}


class TextConsole implements RMainLoopCallbacks
{
    public void rWriteConsole(Rengine re, String text, int oType) {
        System.out.print(text);
    }

    public void rBusy(Rengine re, int which) {
        System.out.println("rBusy("+which+")");
    }

    public String rReadConsole(Rengine re, String prompt, int addToHistory) {
        System.out.print(prompt);
        try {
            BufferedReader br=new BufferedReader(new InputStreamReader(System.in));
            String s=br.readLine();
            return (s==null||s.length()==0)?s:s+"\n";
        } catch (Exception e) {
            System.out.println("jriReadConsole exception: "+e.getMessage());
        }
        return null;
    }

    public void rShowMessage(Rengine re, String message) {
        System.out.println("rShowMessage \""+message+"\"");
    }

    public String rChooseFile(Rengine re, int newFile) {
	FileDialog fd = new FileDialog(new java.awt.Frame(), (newFile==0)?"Select a file":"Select a new file", (newFile==0)?FileDialog.LOAD:FileDialog.SAVE);
	fd.show();
	String res=null;
	if (fd.getDirectory()!=null) res=fd.getDirectory();
	if (fd.getFile()!=null) res=(res==null)?fd.getFile():(res+fd.getFile());
	return res;
    }

    public void   rFlushConsole (Rengine re) {
    }

    public void   rLoadHistory  (Rengine re, String filename) {
    }

    public void   rSaveHistory  (Rengine re, String filename) {
    }
}

