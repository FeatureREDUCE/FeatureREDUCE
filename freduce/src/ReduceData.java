/*
 * ReduceData.java - Todd Riley
 *
 */

//package org.biojava.bio.mydp;

import java.util.*;
import java.util.regex.Pattern;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;
import java.io.File;

import org.biojava.bio.*;
import org.biojava.bio.dist.*;
import org.biojava.bio.dp.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.db.*;
import org.biojava.bio.seq.io.*;
import org.biojava.bio.symbol.*;
import org.biojava.utils.*;
import org.biojava.bio.seq.impl.SimpleSequence;

import org.biojava.bio.mydp.*;
import org.biojava.bio.mydp.FileTools;

/*
The ReduceData drives all access to the probe sequence, intensitities, and motif counts data.

There are 3 (LinkedHashMap) caches:

Cache on ProbeID: (retrieves ProbeData)
   Normalized probe intensity
   word-count suffix tree for pos strand
   word-count suffix tree for neg strand

Cache on a Word (SymboList): (retrives MotifData)
   word-count hash map for pos strand (key=ProbeID)
   word-count hash map for neg strand (key=ProbeID)

Cache on a WeightMatrix Label (String): (retrieves MotifData)
   word-count hash map for pos strand (key=ProbeID)
   word-count hash map for neg strand (key=ProbeID)

Hash Map objects are saved as (human readable) ascii key-value files
and as serialized binary files
Sequence Data for the suffix trees can be retrieved from a fasta file of a Probe Map File.

Current Probe Map Files supported:
 BPMAP text file - Affymettrix

 *.bpmap.txt file = PMX [tab] PMY [tab] MMX [tab] MMY [tab] Seq [tab] Pos [tab] ProbeSequence
                    0         1         2         3         4         5         6

You can create the *.bpmap.txt file from a *.bpmap file with the ProbExporter Tool (from the Affy website).

The ReduceData class can also separate the data into *.fasta and *.bpmap.intensities files to conform to the
MatrixREDUCE input file formats.

 */

public class ReduceData {

    ////////////////////////////////////////////////////////////////////////////////
    // Globals
    ////////////////////////////////////////////////////////////////////////////////

    // init file cache
    KeyValueMap initFileCache = null;

    // intensities Linked HashMaps
    LinkedHashMap<Integer, LinkedHashMap<String, Double>> expNumToOrderedUsedIntensities = null;

    // SeqIDs Set
    List<Set<String>> allProbeIDs = null;

    // symbol matrix hash map
    LinkedHashMap<Integer, SymbolMatrix> expNumToSymbolMatrix = null;

    // kmer counts hash map
    LinkedHashMap<Integer, KmerMatrix> expNumKmerToKmerMatrix= null;

    // SymbolList to count array hash map
    LinkedHashMap<SymbolList, double[]> wordCountsHM;

    // WeightMatrix to count array hash map
    LinkedHashMap<WeightMatrix, double[]> weightMatrixCountsHM;

    // flag to cache or not cache counts arrays
    boolean cacheCountsArrays = false;

    // used by the PackedSymbolLists in the SymbolMatrix objects
    private Packing packing;

    int numExperiments = 0;

    // There can be only one probeMapFile per experiment (for all replicates)
    // OR there is exactly one probeMapFile per replicate
    String  probeMapFilePathNames[] = null;

    // There are potentially many test and control intensities files per experiment
    String  probeIntensityFiles[] = null;

    private static String assigner = "=";
    private static String delimiter = ", ";
    private static String comma = ",";
    private static String semicolon = ";";
    private static String commaWhiteSpace = "[,\\s]";
    private static String whiteSpace = "\\s+";
    private static String tabs = "\\t+";
    private static String tab = "\\t";
    private static Pattern assignerPattern = Pattern.compile(assigner);
    private static Pattern delimiterPattern = Pattern.compile(delimiter);
    private static Pattern commaPattern = Pattern.compile(comma);
    private static Pattern semicolonPattern = Pattern.compile(semicolon);
    private static Pattern commaWhiteSpacePattern = Pattern.compile(commaWhiteSpace);
    private static Pattern whiteSpacePattern = Pattern.compile(whiteSpace);
    private static Pattern tabsPattern = Pattern.compile(tabs);
    private static Pattern tabPattern = Pattern.compile(tab);

    Output out = null;
    String outputLevelString = null;

    String alphabetName = "DNA";
    Alphabet alphabet = null;
    SymbolTokenization symbolTokenization = null;
    AlphabetIndex indexer;
    ReversibleTranslationTable complementTable = null;

    String genomeDirString = null;
    UCSCGenome ucscGenome = null;

    // The columns for the Probe Map File (e.g. - bmpap file)
    int mapSeqLabelColumn = -1;
    int mapSeqIndexColumn = -1;
    int mapXColumn = -1;
    int mapYColumn = -1;

    // The columns for the Intensity File (e.g. - cel file)
    int intensIntensityColumn = -1;
    int intensSeqLabelColumn = -1;
    int intensSeqIndexColumn = -1;
    int intensXColumn = -1;
    int intensYColumn = -1;

    ////////////////////////////////////////////////////////////////////////////////
    // Constructors
    ////////////////////////////////////////////////////////////////////////////////

    public ReduceData(KeyValueMap anInitFileCache) {

        try {
            initialize(anInitFileCache);

            LinkedHashMap<String, Double> allProbeMapIntensities = null;

            String intensInputFileType = (String)initFileCache.get("INTENSITIES", "InputFileType");
            String intensOutputFilePrefix = (String)initFileCache.get("INTENSITIES", "OutputFilePrefix");
            String intensKeySeparator = (String)initFileCache.get("INTENSITIES", "KeySeparator");
            String intensGrouping = (String)initFileCache.get("INTENSITIES", "Grouping");
            String probeIntensitiesDir = (String)initFileCache.get("DATA", "ProbeIntensitiesDir");
            String numberOfProbesString = (String)initFileCache.get("REGRESSION", "NumberOfProbes");

            String probeLocsDir = (String)initFileCache.get("DATA", "ProbeLocsDir");
            String seqKeySeparator = (String)initFileCache.get("SEQUENCE", "KeySeparator");

            int experimentNumber = 1;
            String probeMapFilePrefix;
            while ((probeMapFilePrefix = (String)initFileCache.get("SEQUENCE", Integer.toString(experimentNumber))) != null) {

                numExperiments++;

                if (initFileCache.get("DREAM", "Label") != null) {
                    String idString = (String)initFileCache.get("DREAM", "ID");
                    out.print("\nRetrieving all probe IDs for protein "+idString+" from "+probeMapFilePrefix+" file....");

                    String[] allLines = FileTools.readStrings(
                        probeMapFilePrefix,
                        1,
                        (String)null,
                        -1,
                        "\\t",
                        initFileCache.getInt("DREAM", "LabelColumn"), // match column
                        idString);

                    //out.println("\nallLines length = " + allLines.length);

                    ArrayList<String> probeIDs = StringTools.getColumn(
                        allLines,
                        tabPattern,
                        null,
                        initFileCache.getInt("DREAM", "SequenceColumn"), // get column
                        initFileCache.getInt("DREAM", "LabelColumn"), // match column
                        idString); // match String

                    //out.println("\nprobeIDs size = " + probeIDs.size());
                    //out.println("\nprobeID[0] = " + probeIDs.get(0));

                    setAllProbeIDs(experimentNumber, new LinkedHashSet(probeIDs));
                    out.println("Done.");
                    break;
                }

                //out.println("\nHere 1!.\n");
                if (intensInputFileType.equalsIgnoreCase(".unified") || intensInputFileType.equalsIgnoreCase(".unifiedLabels")) {
                    out.print("\nRetrieving all probe IDs from "+probeMapFilePrefix+" file....");
                    String[] allLines = FileTools.readStrings(probeMapFilePrefix, 1, null);
                    out.println("\nallLines[0] = "+allLines[0]+"\n");
                    ArrayList<String> probeIDs = StringTools.getColumn(allLines, "\t", null, 0);
                    setAllProbeIDs(experimentNumber, new LinkedHashSet(probeIDs));
                    out.println("\nprobeIDs[0] = "+probeIDs.get(0)+"\n");
                    out.println("\nprobeIDs[1] = "+probeIDs.get(1)+"\n");
                    out.println("Done.");
                    break;
                }

                //String probeMapFilePrefix = (String)initFileCache.get("SEQUENCE", Integer.toString(experimentNumber));
                String intensFilePrefixes = FileTools.stripPath((String)initFileCache.get("INTENSITIES", Integer.toString(experimentNumber)));
                String probeIntensitiesFilePrefix = probeIntensitiesDir+File.separator+intensOutputFilePrefix+"."+Integer.toString(experimentNumber);

                String allProbeMapIntensitiesFile = probeIntensitiesFilePrefix+".bpmap.norm.all.intensities";
                String useProbeMapIntensitiesFile = probeIntensitiesFilePrefix+".bpmap.norm.use.intensities";

                // Find out what Probe Map File we have
                String bpmapFileString = probeLocsDir+File.separator+probeMapFilePrefix+".bpmap.txt";
                if (FileTools.exists(bpmapFileString)) {
                    this.probeMapFilePathNames = new String[1];
                    this.probeMapFilePathNames[experimentNumber-1] = bpmapFileString;
                    this.mapSeqLabelColumn = 4;
                    this.mapSeqIndexColumn = 5;
                    this.mapXColumn = 0;
                    this.mapYColumn = 1;
                }

                ////////////////////////////////////////////////////////////////////////////////
                // create <experimentNumber>.bpmap.all.intensities file if it doesn't exist
                ////////////////////////////////////////////////////////////////////////////////
                if (!FileTools.exists(allProbeMapIntensitiesFile)) {
                    out.println("\nDid not find existing "+allProbeMapIntensitiesFile+" file. Creating one.");

                    if (intensInputFileType.equalsIgnoreCase(".bar.txt")) {
                        // write out <experimentNumber>.bpmap.intensities file from the one *.bar.txt for the experiment
                        allProbeMapIntensities = getProbeMapIntensities(probeIntensitiesDir, intensFilePrefixes, seqKeySeparator, probeIntensitiesFilePrefix);
                    }
                    else {
                        String[] intensFilePrefixesArray = whiteSpacePattern.split(intensFilePrefixes);
                        LinkedHashMap<String, Double> XYtoIntensity = getAverageXYIntensities(probeIntensitiesDir, intensFilePrefixesArray, intensGrouping, intensKeySeparator, intensInputFileType);
                        // write out <experimentNumber>.xy.intensities file
                        FileTools.write(XYtoIntensity, probeIntensitiesFilePrefix+".xy.all.intensities", "\t", "");

                        if (this.probeMapFilePathNames != null) {
                            allProbeMapIntensities = getProbeMapIntensities(XYtoIntensity, probeIntensitiesFilePrefix, seqKeySeparator, intensKeySeparator, experimentNumber);
                        }
                        else {
                            allProbeMapIntensities = XYtoIntensity;
                        }

                    }
                    //allProbeIDs = allProbeMapIntensities.keySet();
                    setAllProbeIDs(experimentNumber, allProbeMapIntensities.keySet());

                    // write out <experimentNumber>.bpmap.intensities file
                    FileTools.write(allProbeMapIntensities, allProbeMapIntensitiesFile, "\t", "");
                }
                else {
                    out.println("\nFound existing "+allProbeMapIntensitiesFile+" file. Not creating a new one.");
                    String[] allLines = FileTools.readStrings(allProbeMapIntensitiesFile, 1, null);
                    ArrayList<String> probeIDs = StringTools.getColumn(allLines, "\t", null, 0);
                    setAllProbeIDs(experimentNumber, new LinkedHashSet(probeIDs));
                }

//                 ////////////////////////////////////////////////////////////////////////////////
//                 // create <prefix>.<experimentNumber>.bpmap.norm.use.intensities file if it doesn't exist
//                 ////////////////////////////////////////////////////////////////////////////////
//                 if (!FileTools.exists(useProbeMapIntensitiesFile)) {
//                     out.print("\nDid not find existing "+useProbeMapIntensitiesFile+" file. Creating one...");

//                     // Sort the allProbeMapIntensities LinkedHashMap by populating a Table and then sorting
//                     if (allProbeMapIntensities == null) {
//                         allProbeMapIntensities = FileTools.readStringDoubleLHM(allProbeMapIntensitiesFile, "\t", 0);
//                         allProbeIDs = allProbeMapIntensities.keySet();
//                     }

//                     Table sortedProbeMapIntensities = new Table(allProbeMapIntensities);
//                     sortedProbeMapIntensities.sort(1); //sort by ascending values
//                     sortedProbeMapIntensities.reverse(); // get descending order

//                     int numberOfUsedProbes = 0;
//                     if (numberOfProbesString.equalsIgnoreCase("ALL")) {
//                         numberOfUsedProbes = allProbeMapIntensities.size();
//                     }
//                     else {
//                         numberOfUsedProbes = Integer.parseInt(numberOfProbesString);
//                     }

//                     FileTools.write(sortedProbeMapIntensities.toString(0, numberOfUsedProbes, "\t"), useProbeMapIntensitiesFile, false);
//                     out.println(" Done.");
//                     out.println("\nThe top "+numberOfUsedProbes+" sorted probe intensities were written to "+useProbeMapIntensitiesFile+".");
//                 }
//                 else {
//                     out.println("\nFound existing "+useProbeMapIntensitiesFile+" file. Not creating a new one.");
//                 }

                // Free allProbeMapIntensities
                allProbeMapIntensities = null;
                //allProbeIDs = null;

                experimentNumber++;

            } // end-while

            out.println("\nFound "+numExperiments+" experiments to fit. The coefficients from each fit(experiment) will be averaged.");

        }

        catch (Exception ex) {
            ex.printStackTrace();
        }

    }

    public void putOrderedUsedIntensities(int experimentNumber, int unifiedColumn, LinkedHashMap<String, Double> orderedUsedIntensities) {
        if (expNumToOrderedUsedIntensities == null) {
             expNumToOrderedUsedIntensities = new LinkedHashMap<Integer, LinkedHashMap<String, Double>>();
        }
        expNumToOrderedUsedIntensities.put(experimentNumber*10000 + unifiedColumn, orderedUsedIntensities);
    }

    // Intensity files must have a column labels header line
    public LinkedHashMap<String, Double> getOrderedUsedIntensities(int experimentNumber, int unifiedColumn) {
        try {

            LinkedHashMap<String, Double> orderedUsedIntensities;

            if (expNumToOrderedUsedIntensities == null) {
                expNumToOrderedUsedIntensities = new LinkedHashMap<Integer, LinkedHashMap<String, Double>>();
            }

            if ((orderedUsedIntensities = expNumToOrderedUsedIntensities.get(experimentNumber*10000 + unifiedColumn)) == null) {

                ///////////////////////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////
                // Clear the HashMap if we don't have enough RAM to keep them all in the cache!!!!
                ///////////////////////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////
                expNumToOrderedUsedIntensities.clear();

                ////////////////////////////////////////////////////////////////////////////////
                // Load the <prefix>.<experimentNumber>.bpmap.norm.use.intensities file into orderedUsedIntensities
                ////////////////////////////////////////////////////////////////////////////////
                String inputFileType = (String)initFileCache.get("INTENSITIES", "InputFileType");
                String numberOfProbesString = (String)initFileCache.get("REGRESSION", "NumberOfProbes");
                String allProbeMapIntensitiesFile;
                int keyColumn = 0;
                int valueColumn = 1;
                int matchColumn = -1;
                String matchString = null;

                out.println("\nunifiedColumn="+unifiedColumn+".");

                if (initFileCache.get("DREAM", "Label") != null) {
                    keyColumn = initFileCache.getInt("DREAM", "SequenceColumn");
                    valueColumn = initFileCache.getInt("DREAM", "IntensityColumn");
                    matchColumn = initFileCache.getInt("DREAM", "LabelColumn");
                    matchString = (String)initFileCache.get("DREAM", "ID");
                    allProbeMapIntensitiesFile = (String)initFileCache.get("INTENSITIES", Integer.toString(experimentNumber));
                }
                else if (inputFileType.equalsIgnoreCase(".unified") || inputFileType.equalsIgnoreCase(".unifiedLabels")) {
                    if (unifiedColumn >= 0) {
                        valueColumn = unifiedColumn;
                    }
                    allProbeMapIntensitiesFile = (String)initFileCache.get("INTENSITIES", Integer.toString(experimentNumber));
                }
                else {
                    String probeIntensitiesDir = (String)initFileCache.get("DATA", "ProbeIntensitiesDir");
                    String intensOutputFilePrefix = FileTools.stripPath((String)initFileCache.get("INTENSITIES", "OutputFilePrefix"));
                    String probeIntensitiesFilePrefix = probeIntensitiesDir+File.separator+intensOutputFilePrefix+"."+Integer.toString(experimentNumber);
                    allProbeMapIntensitiesFile = probeIntensitiesFilePrefix+".bpmap.norm.all.intensities";
                }

                //String useProbeMapIntensitiesFile = probeIntensitiesFilePrefix+".bpmap.norm.use.intensities";

//                 // initialize as same size as keyIntensityPairs
//                 orderedUsedIntensities = FileTools.readStringDoubleLHM(useProbeMapIntensitiesFile, "\t", 0);
                // Sort the allProbeMapIntensities LinkedHashMap by populating a Table and then sorting

                ////////////////////////////////////////////////////////////////////////////////
                // create sorted list of intensities and keep only top N
                ////////////////////////////////////////////////////////////////////////////////
                out.print("\nCaching all probe intensities from column "+valueColumn+" in file "+allProbeMapIntensitiesFile+"...");
                //orderedUsedIntensities = FileTools.readStringDoubleLHM(allProbeMapIntensitiesFile, "\t", 1);

                orderedUsedIntensities = FileTools.readStringDoubleLHM(
                    allProbeMapIntensitiesFile,
                    "\t",
                    1,
                    -1,
                    -1,
                    keyColumn,
                    valueColumn,
                    matchColumn,
                    matchString);

                out.println(" Done.");
                out.println("\n"+orderedUsedIntensities.size()+" probe intensities were read from "+allProbeMapIntensitiesFile+".");

                out.print("\nSorting the intensities...");
                Table sortedProbeMapIntensities = new Table(orderedUsedIntensities);
                sortedProbeMapIntensities.sort(1); //sort by ascending values
                sortedProbeMapIntensities.reverse(); // get descending order
                out.println(" Done.");

                List<Double> intensitiesColumn = (List<Double>)(List)sortedProbeMapIntensities.getColumn(1);
                // double sum = MathTools.sum(intensitiesColumn);
                // if (MathTools.isInteger(sum)) {
                //     initFileCache.put("INTENSITIES", "FinalRoundTotalCounts", (int)sum);
                // }

                int numberOfUsedProbes = 0;
                if (numberOfProbesString.equalsIgnoreCase("ALL")) {
                    numberOfUsedProbes = orderedUsedIntensities.size();
                }
                else {
                    numberOfUsedProbes = Integer.parseInt(numberOfProbesString);
                }

                if (numberOfUsedProbes < orderedUsedIntensities.size()) {
                    out.print("\nKeeping only the top "+numberOfUsedProbes+" probe intensities...");
                    sortedProbeMapIntensities.removeRows(numberOfUsedProbes, sortedProbeMapIntensities.rows());

//                     List aRow = sortedProbeMapIntensities.getRow(numberOfUsedProbes);
//                     double minIntensity = ((Double)aRow.get(1)).doubleValue();

//                     for (Iterator iter = orderedUsedIntensities.values().iterator(); iter.hasNext(); ) {
//                         double intensity = ((Double)iter.next()).doubleValue();
//                         if (intensity < minIntensity) {
//                             iter.remove();// avoids ConcurrentModificationException
//                         }
//                     }
                    out.println(" Done.");
                }
                orderedUsedIntensities = sortedProbeMapIntensities.toMap(0, 1);
                putOrderedUsedIntensities(experimentNumber, unifiedColumn, orderedUsedIntensities);
            }
            return(orderedUsedIntensities);
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(null);
    }

    private void initialize(KeyValueMap anInitFileCache) {
        try {
            // Get parameters from the init file cache
            this.initFileCache = anInitFileCache;

            this.outputLevelString = (String)initFileCache.get("OUTPUT", "Level");
            if (outputLevelString.equalsIgnoreCase("Verbose")) {
                out = new Output(Output.Destination.STDOUT, (byte)25, false);
            }
            else {
                out = new Output(Output.Destination.STDOUT, (byte)100, false);
            }

            this.alphabetName = (String)initFileCache.get("SEQUENCE", "Alphabet");
            // Set the alphabet and complementTable
            if (alphabetName == null || alphabetName.equalsIgnoreCase("DNA")) {
                this.alphabet = DNATools.getDNA();
                this.complementTable = DNATools.complementTable();
            }
            else if (alphabetName.equalsIgnoreCase("RNA")) {
                this.alphabet = RNATools.getRNA();
                this.complementTable = RNATools.complementTable();
            }
            else if (alphabetName.equalsIgnoreCase("Protein")) {
                this.alphabet = ProteinTools.getAlphabet();
            }
            this.symbolTokenization = alphabet.getTokenization("token");
            this.indexer = AlphabetManager.getAlphabetIndex((FiniteAlphabet)alphabet);

            this.packing = PackingFactory.getPacking((FiniteAlphabet)this.alphabet, false);

        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Methods
    ////////////////////////////////////////////////////////////////////////////////

    public void setCacheCountsArraysFlag(boolean aFlag) {
        //this.cacheCountsArrays = aFlag;
    }

    public ArrayList<Sequence> getSequences(String seqKey, SequenceDB aSequenceDB) {
        ArrayList<Sequence> sequences = new ArrayList<Sequence>(1);
        Sequence probeSeq = null;

        try {
            probeSeq = aSequenceDB.getSequence(seqKey);
            sequences.add(probeSeq);
        }
        catch (Exception ex) {
            ex.printStackTrace();
            System.exit(-1);
        }

        return(sequences);
    }

    public ArrayList<Sequence> getSequences(String seqKey, LinkedHashMap<String, String> aSequenceDB) {
        boolean reverseSequence = (boolean)initFileCache.getBoolean("SEQUENCE", "ReverseSequence");
        ArrayList<Sequence> sequences = new ArrayList<Sequence>(1);
        Sequence probeSeq = null;

        try {
            String seqValue = aSequenceDB.get(seqKey);
            if (reverseSequence) {
                seqValue = StringTools.reverse(seqValue);
            }

            //if (seqValue == null) {
            //    out.println("\nError: Could not get sequence for Key:"+seqKey+"\n");
           // }
           // else if (seqValue.equals("")) {
           //     out.println("\nError: Retrieved sequence for Key="+seqKey+" is an empty string.\n");
           // }

            // out.println("\nseqKey="+seqKey+"\n");
            // out.println("\nseqValue="+seqValue+"\n");

            // String randomDNA = "ACGT'";
            // String randomDNA = "";
            // out.println("\nrandomDNA="+randomDNA+"\n");

            // DNATools.createDNASequence(sequence, name);
            if( !seqValue.equals("") && !(seqValue == null)){
            	probeSeq = DNATools.createDNASequence(seqValue, seqKey);
            //probeSeq = DNATools.createDNASequence(seqKey, seqKey);
            //probeSeq = DNATools.createDNASequence(seqValue, seqValue);
            //probeSeq = DNATools.createDNASequence(randomDNA, randomDNA);

            	if (probeSeq == null) {
                	out.println("\nseqKey="+seqKey+"\n");
                	out.println("\nseqValue="+seqValue+"\n");
           	 }

            	sequences.add(probeSeq);
		}
        }
        catch (Exception ex) {
            ex.printStackTrace();
            System.exit(-1);
        }

        return(sequences);
    }

    public enum ProbeLoc {CENTER, LEFT, RIGHT}

    public ArrayList<Sequence> getSequences(String[] probeMapFileLine, ProbeLoc aProbeLoc, boolean anIncludeProbeFlag, int aBandWidth, int aProbeWidth) {
        //return(getSequences(probeMapFileLine[mapSeqLabelColumn], probeMapFileLine[mapSeqIndexColumn], aProbeLoc, anIncludeProbeFlag, aBandWidth, aProbeWidth));

        return(getSequences(probeMapFileLine[0], probeMapFileLine[1], aProbeLoc, anIncludeProbeFlag, aBandWidth, aProbeWidth));
    }

    // aBandWidth is always conserved (the returned sequences are always within the bandwidth of the probe
    // seqIndex is assumed to point to the left-hand starting nucleotide of the probe-sequence
    public ArrayList<Sequence> getSequences(String seqLabel, String seqIndex, ProbeLoc aProbeLoc, boolean anIncludeProbeFlag, int aBandWidth, int aProbeWidth) {

        ArrayList<Sequence> sequences = new ArrayList<Sequence>(2);
        Sequence probeSeq = null;

        try {
            // If this.ucscGenome == null then we must create one
            this.genomeDirString = (String)initFileCache.get("DATA", "GenomeDir");
            if (this.ucscGenome == null) {
                out.print("\nCreating UCSCGenome object to access sequence from "+this.genomeDirString+" ...");
                this.ucscGenome = new UCSCGenome(this.genomeDirString);
                out.println(" Done.");
            }

            int start, end;

            // seqLabel can be either chr[1-9+,X,Y] or <species>:<build>;chr[1-9+,X,Y]
            String[] seqLabelParts = semicolonPattern.split(seqLabel);
            String chrNumString = null;
            if (seqLabelParts.length == 1) {
                // seqLabel is chr[1-9,X,Y]
                chrNumString = seqLabel;
            }
            else {
            // seqLabel is <species>:<build>;chr[1-9,X,Y]
                chrNumString = seqLabelParts[1];
            }

            // if the chrNumString does not start with "chr" then it is an affy control probe!!
            if (!chrNumString.startsWith("chr")) {
                //out.println("chrNumString="+chrNumString);
                return(sequences);
            }

            // remove the "chr"
            String chrRemovedString = chrNumString.substring(3); // remove "chr"

            int probeIndex = Integer.parseInt(seqIndex);

            switch (aProbeLoc) {
            case CENTER:
                if (anIncludeProbeFlag == true) {
                    start = (Integer.parseInt(seqIndex) +1) + (int)Math.round(Math.floor((float) aProbeWidth /2)) - (int)Math.round(Math.floor((float) aBandWidth /2));
                    end = (Integer.parseInt(seqIndex) +1) + (int)Math.round(Math.floor((float) aProbeWidth /2)) + (int)Math.round(Math.floor((float) aBandWidth /2));
                    probeSeq = this.ucscGenome.getSequence("Chromosome", chrRemovedString, start, end, this.alphabet, false);

                    if (probeSeq == null) {
                        out.println("\nError: Could not get sequence from UCSCGenome object: "+seqLabel+", "+start+"-"+end+".");
                    }
                    else {
                        sequences.add(probeSeq);
                    }
                }
                else { // Don't include the probe sequence
                    // Get the left-hand-side
                    start = (Integer.parseInt(seqIndex) +1) + (int)Math.round(Math.floor((float) aProbeWidth /2)) - (int)Math.round(Math.floor((float) aBandWidth /2));
                    end = (Integer.parseInt(seqIndex) +1);
                    probeSeq = this.ucscGenome.getSequence("Chromosome", chrRemovedString, start, end, this.alphabet, false);

                    if (probeSeq == null) {
                        out.println("\nError: Could not get sequence from UCSCGenome object: "+seqLabel+", "+start+"-"+end+".");
                    }
                    else {
                        sequences.add(probeSeq);
                    }

                    // Get the right-hand-side
                    start = (Integer.parseInt(seqIndex) +1) + aProbeWidth;
                    end = (Integer.parseInt(seqIndex) +1) + (int)Math.round(Math.floor((float) aProbeWidth /2)) + (int)Math.round(Math.floor((float) aBandWidth /2));
                    probeSeq = this.ucscGenome.getSequence("Chromosome", chrRemovedString, start, end, this.alphabet, false);

                    if (probeSeq == null) {
                        out.println("\nError: Could not get sequence from UCSCGenome object: "+seqLabel+", "+start+"-"+end+".");
                    }
                    else {
                        sequences.add(probeSeq);
                    }
                }
                break;
            case LEFT:
                if (anIncludeProbeFlag == true) {
                    start = (Integer.parseInt(seqIndex) +1);
                    end = (Integer.parseInt(seqIndex) +1) + aBandWidth;
                }
                else { // Don't include the probe sequence
                    start = (Integer.parseInt(seqIndex) +1) + aProbeWidth;
                    end = (Integer.parseInt(seqIndex) +1) + aBandWidth;
                }

                probeSeq = this.ucscGenome.getSequence("Chromosome", chrRemovedString, start, end, this.alphabet, false);

                if (probeSeq == null) {
                    out.println("\nError: Could not get sequence from UCSCGenome object: "+seqLabel+", "+start+"-"+end+".");
                }
                else {
                    sequences.add(probeSeq);
                }

                break;
            case RIGHT:
                if (anIncludeProbeFlag == true) {
                    start = (Integer.parseInt(seqIndex) +1) + aProbeWidth - aBandWidth;
                    end = (Integer.parseInt(seqIndex) +1) + aProbeWidth;
                }
                else { // Don't include the probe sequence
                    start = (Integer.parseInt(seqIndex) +1) + aProbeWidth - aBandWidth;
                    end = (Integer.parseInt(seqIndex) +1);
                }

                probeSeq = this.ucscGenome.getSequence("Chromosome", chrRemovedString, start, end, this.alphabet, false);

                if (probeSeq == null) {
                    out.println("\nError: Could not get sequence from UCSCGenome object: "+seqLabel+", "+start+"-"+end+".");
                }
                else {
                    sequences.add(probeSeq);
                }

                break;
            default: //Error
                out.println("\nError : could not get sequence because of unsupported probe location: "+aProbeLoc.toString()+".");
                return(null);
            }

        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

        return(sequences);
    }

    // Input: *.bmap.txt file
    // Output: *.fasta file

    // IF YOU REMOVE THE PROBE SEQUENCE FROM A CENTERED BANDWIDTH THEN THE FASTA FILE WILL HAVE TWO SEQUENCES WITH THE SAME SEQID
    // WHICH IS NOT GOOD
    //
    // DON'T REMOVE THE PROBE SEQUENCE FROM A CENTERED BANDWIDTH
    public void createFastaFile(String aProbeMapFilePathName, String outDirString, String outFilePrefixString, int aBandWidth, int aProbeWidth, String probeLocation, String keySeparator, boolean includeProbeSequence) {

        String outputFullPathName = new String(outDirString+File.separator+outFilePrefixString+".fasta");

        try {
            // get the probeLoc
            ProbeLoc probeLoc;
            if (probeLocation.equalsIgnoreCase("Center")) {
                probeLoc = ProbeLoc.CENTER;
            }
            else if (probeLocation.equalsIgnoreCase("Left")) {
                probeLoc = ProbeLoc.LEFT;
            }
            else if (probeLocation.equalsIgnoreCase("Right")) {
                probeLoc = ProbeLoc.RIGHT;
            }
            else { //Error
                out.println("\nError : could not get sequence because of unsupported [SEQUENCE]->ProbeLocation "+probeLocation+".");
                return;
            }

            String lineString;
            String[] lineFields;
            ArrayList<Sequence> probeSeqs;

            BufferedReader inputBuffer = new BufferedReader(new FileReader(aProbeMapFilePathName));
            BufferedWriter outputBuffer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFullPathName)));

            out.print("\nCreating "+outputFullPathName+" from "+aProbeMapFilePathName+": BandWidth="+aBandWidth+" probeLocation="+probeLocation+" includeProbeSequence="+includeProbeSequence+" ...");

            int numSeqsRetrieved = 0;
            while((lineString = inputBuffer.readLine()) != null) {
                // include this line if line doesn't start with commentToken "#"
                if (!lineString.trim().startsWith("#") && !lineString.trim().startsWith("PMX")) {

                    lineFields = tabsPattern.split(lineString);
                    probeSeqs = this.getSequences(lineFields, probeLoc, includeProbeSequence, aBandWidth, aProbeWidth);

                    for (Sequence probeSeq : probeSeqs) {
                        if (probeSeq != null) {
                            // Write sequence in FASTA format
                            //outputBuffer.write("\n"+">"+lineFields[mapSeqLabelColumn].trim()+keySeparator+lineFields[mapSeqIndexColumn].trim()+" BandWidth="+aBandWidth+" probeLocation="+probeLocation+" includeProbeSequence="+includeProbeSequence+"\n"+probeSeq.seqString()+"");
                            outputBuffer.write("\n"+">"+lineFields[0].trim()+keySeparator+lineFields[1].trim()+" BandWidth="+aBandWidth+" probeLocation="+probeLocation+" includeProbeSequence="+includeProbeSequence+"\n"+probeSeq.seqString()+"");

                            numSeqsRetrieved++;
                            if((numSeqsRetrieved % 1000) == 0 ) {
                                out.println("\nRetrieved "+numSeqsRetrieved+" probe-associated sequences so far...");
                            }

                        }
                    }
                }
            }

            out.println(" Done.");

            outputBuffer.close();
            inputBuffer.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }

    public double[] getUsedIntensities(int experimentNumber, int unifiedColumn) {
        LinkedHashMap<String, Double> orderedSubsetIntensities = getOrderedUsedIntensities(experimentNumber, unifiedColumn);

        double[] intensities = new double[orderedSubsetIntensities.size()];
        int i = 0;
        for (Double intensity : orderedSubsetIntensities.values()) {
            intensities[i] = intensity;
            i++;
        }
        return(intensities);
    }

    public int numberOfUsedProbes(int experimentNumber, int unifiedColumn) {
        LinkedHashMap<String, Double> orderedSubsetIntensities = getOrderedUsedIntensities(experimentNumber, unifiedColumn);

        return(orderedSubsetIntensities.size());
    }

    public Set<String> getUsedProbeIDs(int experimentNumber, int unifiedColumn) {
        LinkedHashMap<String, Double> orderedSubsetIntensities = getOrderedUsedIntensities(experimentNumber, unifiedColumn);

        return(orderedSubsetIntensities.keySet());
    }

    public Set<String> getAllProbeIDs(int experimentNumber) {
        return(allProbeIDs.get(experimentNumber-1));
    }

    public void setAllProbeIDs(int experimentNumber, Set<String> probeIDs) {
        if (allProbeIDs == null) {
            allProbeIDs = new ArrayList<Set<String>>();
        }
        if (experimentNumber == allProbeIDs.size()+1) {
            allProbeIDs.add(probeIDs);
        }
        else {
            allProbeIDs.set(experimentNumber-1, probeIDs);
        }
    }

    public int getNumExperiments() {
        return(numExperiments);
    }

    public int getSeqLength() {
        return(getSeqLength(1, -1));
    }

    public int getSeqLength(int experimentNumber, int unifiedColumn) {
        int seqLength = initFileCache.getInt("SEQUENCE", "SeqLength");

        if (seqLength <= 0) {
            seqLength = getSymbolMatrix(experimentNumber, unifiedColumn).getSeqLength();
        }

        return(seqLength);

        // return(35);

        // SymbolMatrix aSymbolMatrix = expNumToSymbolMatrix.values().iterator().next();
        // return(aSymbolMatrix.getSeqLength());

        // return(getSymbolMatrix(1,-1).getSeqLength());
        // return(getSymbolMatrix(1,2).getSeqLength());
    }

    public WeightMatrixTools.BindingStrand getStrand() {
            String strandString = (String)initFileCache.get("REGRESSION", "Strand");
            return(getStrand(strandString));
    }

    public static WeightMatrixTools.BindingStrand getStrand(String strandString) {
            WeightMatrixTools.BindingStrand strand;

            if (strandString.equals("+") || strandString.equalsIgnoreCase("pos")) {
                strand = WeightMatrixTools.BindingStrand.POS;
            }
            else if (strandString.equals("-")  || strandString.equalsIgnoreCase("neg")) {
                strand = WeightMatrixTools.BindingStrand.NEG;
            }
            else {
                strand = WeightMatrixTools.BindingStrand.BOTH;
            }

            return(strand);
    }

    public WeightMatrixTools.BothStrandsCalc getBothStrandsCalc() {
            String calcString = (String)initFileCache.get("REGRESSION", "BothStrandsCalc");
            return(getBothStrandsCalc(calcString));
    }

    public static WeightMatrixTools.BothStrandsCalc getBothStrandsCalc(String calcString) {
            WeightMatrixTools.BothStrandsCalc calc = null;

            if (calcString.equalsIgnoreCase("MAX")) {
                calc = WeightMatrixTools.BothStrandsCalc.MAX;
            }
            else if (calcString.equalsIgnoreCase("SUM")) {
                calc = WeightMatrixTools.BothStrandsCalc.SUM;
            }
            else if (calcString.equalsIgnoreCase("UNION")) {
                calc = WeightMatrixTools.BothStrandsCalc.UNION;
            }
            else if (calcString.equalsIgnoreCase("NORMED_SUM")) {
                calc = WeightMatrixTools.BothStrandsCalc.NORMED_SUM;
            }

            return(calc);
    }

    //
    // a motif can be either a SymbolList or a WeightMatrix
    //
    // startPos - non-null means look only at the motif-window starting at startPos
    //
    // positionalWeights - binding bias based on position in the probe sequence
    // positionalWeights[0] = positive strand positionalWeights
    // positionalWeights[1] = negative strand positionalWeights
    //
    public double[] getCountsArray(
        Object aMotif,
        LinkedHashSet<SymbolList> hammingDistHashSet,
        Integer startPos,
        boolean[][] mandatoryColumns,
        double[][] positionalWeights,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        int experimentNumber,
        int unifiedColumn) {

        if (aMotif instanceof SymbolList) {
            return(getCountsArray((SymbolList)aMotif, startPos, positionalWeights, strand, experimentNumber, unifiedColumn));
        }
        else if (aMotif instanceof WeightMatrix) {
            return(getCountsArray((WeightMatrix)aMotif, hammingDistHashSet, startPos, mandatoryColumns, positionalWeights, strand, calc, experimentNumber, unifiedColumn));
        }
        else { // FeaturedWeightMatrix
            double threshold = initFileCache.getDouble("AFFINITY MODELS", "RelativeFeatureAffinityThreshold");

            return(getCountsArray((FeaturedWeightMatrix)aMotif, hammingDistHashSet, startPos, mandatoryColumns, positionalWeights, strand, calc, experimentNumber, unifiedColumn, threshold));
        }
    }

    // mandatoryColumns[0] = positive strand mandatoryColumns
    // mandatoryColumns[1] = negative strand mandatoryColumns
    // positionalWeights[0] = positive strand positionalWeights
    // positionalWeights[1] = negative strand positionalWeights
    public double[] getCountsArray(
        Object aMotif,
        Integer startPos,
        boolean[][] mandatoryColumns,
        double[][] positionalWeights,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc aCalc,
        int experimentNumber,
        int unifiedColumn) {

        if (aMotif instanceof SymbolList) {
            return(getCountsArray((SymbolList)aMotif, startPos, positionalWeights, strand, experimentNumber, unifiedColumn));
        }
        else { // weightMatrix
            return(getCountsArray((WeightMatrix)aMotif, null, startPos, mandatoryColumns, positionalWeights, strand, aCalc, experimentNumber, unifiedColumn));
        }
    }

    public double[] getCountsArray(
        WeightMatrix aPWM,
        LinkedHashSet<SymbolList> hammingDistHashSet,
        Integer startPos,
        boolean[][] mandatoryColumns,
        double[][] positionalWeights,
        WeightMatrixTools.BindingStrand aStrand,
        WeightMatrixTools.BothStrandsCalc aCalc,
        int experimentNumber,
        int unifiedColumn) {

        double[] affinitiesForMotif = null;
        boolean found = true;

        if (cacheCountsArrays) { // use the cache
            if (weightMatrixCountsHM == null) {
                weightMatrixCountsHM = new LinkedHashMap<WeightMatrix, double[]>();
            }
            affinitiesForMotif = weightMatrixCountsHM.get(aPWM);
        }

        if (affinitiesForMotif == null) { // if not caching or not found
            found = false;

            if (hammingDistHashSet == null) {
                SymbolMatrix symbolMatrix = getSymbolMatrix(experimentNumber, unifiedColumn);
                if ((positionalWeights == null) && (startPos == null)) {
                    // no positionalWeights nor a specific start position
                    affinitiesForMotif = symbolMatrix.getProbSums(aPWM, aStrand, aCalc, 0, 0);
                }
                else if (positionalWeights != null) {
                    // has positionalWeights and no specific start position
                    if (mandatoryColumns == null) {
                        mandatoryColumns = new boolean[2][];
                    }
                    affinitiesForMotif = symbolMatrix.getProbSums(aPWM, positionalWeights, mandatoryColumns, aStrand, aCalc, 0, 0);
                }
                else {
                    // has a specific start position
                    // and no positionalWeights and no mandatory columns !!!
                    affinitiesForMotif = symbolMatrix.getProbs(aPWM, aStrand, aCalc, 0, 0, startPos.intValue());
                }
            }
            else { // use hammingDistHashSet for the quick approximation for counts(WeightMatrix)
                //out.print("\nHere5!");
                affinitiesForMotif = new double[numberOfUsedProbes(experimentNumber, 0)];

                WeightMatrix negStrandPWM = null;
                double wmDistance = 1.0;

                if ((aStrand == WeightMatrixTools.BindingStrand.NEG) || (aStrand == WeightMatrixTools.BindingStrand.BOTH)) {
                    negStrandPWM = WeightMatrixTools.reverseComplement(aPWM, complementTable);
                    //wmDistance = WeightMatrixTools.overlap(aPWM, negStrandPWM);
                    wmDistance = WeightMatrixTools.similarity(aPWM, negStrandPWM);
                }

                for (SymbolList aSymList : hammingDistHashSet) {
                    double score = WeightMatrixTools.score(aPWM, negStrandPWM, aSymList, aStrand, aCalc, 0, 0, wmDistance);
                    //out.print("\nscore="+score);

                    double[] affinities = getCountsArray(aSymList, startPos, positionalWeights, aStrand, experimentNumber, unifiedColumn);
                    double[] scaledAffinities = MathTools.multiply(score, affinities);
                    affinitiesForMotif = MathTools.add(affinitiesForMotif, scaledAffinities);
                }

            }
        }

        if (cacheCountsArrays && !found) {
            // cache the motif Counts for this Weight Matrix
            weightMatrixCountsHM.put(aPWM, affinitiesForMotif);
        }

        return(affinitiesForMotif);

    }

    public double[] getCountsArray(
        FeaturedWeightMatrix anFSAM,
        LinkedHashSet<SymbolList> hammingDistHashSet,
        Integer startPos,
        boolean[][] mandatoryColumns,
        double[][] positionalWeights,
        WeightMatrixTools.BindingStrand aStrand,
        WeightMatrixTools.BothStrandsCalc aCalc,
        int experimentNumber,
        int unifiedColumn,
        double threshold) {

        anFSAM.setMandatoryColumns(mandatoryColumns);

        double[] affinitiesForMotif = null;
        SymbolMatrix symbolMatrix = getSymbolMatrix(experimentNumber, unifiedColumn);

        if ((positionalWeights == null) && (startPos == null)) {
            // no positionalWeights nor a specific start position
            affinitiesForMotif = symbolMatrix.getProbSums(anFSAM, aStrand, aCalc, 0, 0, threshold);
        }
        else if (positionalWeights != null) {
            // has positionalWeights and no specific start position
            affinitiesForMotif = symbolMatrix.getProbSums(anFSAM, positionalWeights, aStrand, aCalc, 0, 0, threshold);
        }
        else {
            // has a specific start position
            // and no positionalWeights
            // use the mandatory columns for each feature!!!
            affinitiesForMotif = symbolMatrix.getProbs(anFSAM, aStrand, aCalc, 0, 0, threshold, startPos.intValue());
        }

        return(affinitiesForMotif);
    }


    public double[] getCountsArray(SymbolList aWord,
        Integer startPos,
        double[][] positionalWeights,
        WeightMatrixTools.BindingStrand aStrand,
        int experimentNumber,
        int unifiedColumn) {

        double[] motifCounts = null;
        boolean found = true;

        if (cacheCountsArrays) { // use the cache
            if (wordCountsHM == null) {
                wordCountsHM = new LinkedHashMap<SymbolList, double[]>();
            }
            motifCounts = wordCountsHM.get(aWord);
        }

        if (motifCounts == null) { // if not caching or not found
            found = false;
            KmerMatrix kmerMatrix = getKmerMatrix(aWord.length(), experimentNumber, unifiedColumn);

            if ((positionalWeights == null) && (startPos == null)) {
                motifCounts = kmerMatrix.getCountsArray(aWord, aStrand);
            }
            else if (positionalWeights != null) {
                motifCounts = kmerMatrix.getCountsArray(aWord, positionalWeights, aStrand);
            }
            else {
                motifCounts = kmerMatrix.getCountsArray(aWord, startPos.intValue(), aStrand);
            }
        }

        if (cacheCountsArrays && !found) {
            wordCountsHM.put(aWord, motifCounts);
        }

        return(motifCounts);
    }

    public void putCountsArray(WeightMatrix aPWM, double[] countsArray) {
        weightMatrixCountsHM.put(aPWM, countsArray);
    }

    public void putSymbolMatrix (int experimentNumber, int unifiedColumn, SymbolMatrix aSymbolMatrix) {

        if (expNumToSymbolMatrix == null) {
            expNumToSymbolMatrix = new LinkedHashMap<Integer, SymbolMatrix>();
        }

        int key = experimentNumber*10000 + unifiedColumn;
        expNumToSymbolMatrix.put(key, aSymbolMatrix);
    }

    public SymbolMatrix getSymbolMatrix (int experimentNumber, int unifiedColumn) {
        SymbolMatrix aSymbolMatrix;
        int key = experimentNumber*10000 + unifiedColumn;

        if (expNumToSymbolMatrix == null) {
            expNumToSymbolMatrix = new LinkedHashMap<Integer, SymbolMatrix>();
        }
        if ((aSymbolMatrix = expNumToSymbolMatrix.get(key)) == null) {
            //return(null);

            ///////////////////////////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////
            // Clear the HashMap if we don't have enough RAM to keep them all in the cache!!!!
            ///////////////////////////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////
            expNumToSymbolMatrix.clear();

            aSymbolMatrix = cacheSymbolMatrixFile(experimentNumber, unifiedColumn);
        }
        return(aSymbolMatrix);
    }

    public void putKmerMatrix (int kmerLength, int experimentNumber, int unifiedColumn, KmerMatrix aKmerMatrix) {

        if (expNumKmerToKmerMatrix == null) {
            expNumKmerToKmerMatrix = new LinkedHashMap<Integer, KmerMatrix>();
        }

        int key = (100000*experimentNumber) + (1000*unifiedColumn) + kmerLength;
        expNumKmerToKmerMatrix.put(key, aKmerMatrix);
    }

    public KmerMatrix getKmerMatrix (int kmerLength, int experimentNumber, int unifiedColumn) {
        KmerMatrix aKmerMatrix;
        int key = (100000*experimentNumber) + (1000*unifiedColumn) + kmerLength;

        if (expNumKmerToKmerMatrix == null) {
            expNumKmerToKmerMatrix = new LinkedHashMap<Integer, KmerMatrix>();
        }
        if ((aKmerMatrix = expNumKmerToKmerMatrix.get(key)) == null) {
            //return(null);

            ///////////////////////////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////
            // Clear the HashMap if we don't have enough RAM to keep them all in the cache!!!!
            ///////////////////////////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////
            expNumKmerToKmerMatrix.clear();

            aKmerMatrix = cacheKmerMatrixFile(kmerLength, experimentNumber, unifiedColumn);
        }
        //out.println("Here4! - kmerLength="+kmerLength);
        return(aKmerMatrix);
    }

    // Input: seqIndexToProbeMapData
    // Output: *.fasta file
        enum SeqInput {UCSCGENOME, FASTA, UNIFIED, TABLE}

    //
    // creates the kmersByProbe KmerMatrix File for each kmer-length
    //
    // The ALL probeIDs file exists in the ProbeLocsDir/wordCountsDir
    //
    public KmerMatrix createKmerMatrixFile(int kmerLength, int experimentNumber, int unifiedColumn) {
        //        return(createKmerMatrixFile(kmerLength, getSymbolMatrix(experimentNumber, unifiedColumn), experimentNumber));
        return(createKmerMatrixFile(kmerLength, getSymbolMatrix(experimentNumber, unifiedColumn), experimentNumber));
    }

    public KmerMatrix createKmerMatrixFile(int kmerLength, SymbolMatrix aSymbolMatrix, int experimentNumber) {
        KmerMatrix aKmerMatrix = null;
        try {
            int numSeqIDs = aSymbolMatrix.getNumSeqIDs();

            String probeLocsDir = (String)initFileCache.get("DATA", "ProbeLocsDir");
            String outputDir = (String)initFileCache.get("OUTPUT", "Directory");
            String wordCountsDir = (String)initFileCache.get("OUTPUT", "WordCountsDir");

            // The ALL probeIDs file exists in the ProbeLocsDir/wordCountsDir
            //String allWordCountsDir = probeLocsDir + File.separator + wordCountsDir;

            String proteinID = (String)initFileCache.get("DREAM", "ID");
            String label = (String)initFileCache.get("DREAM", "Label");

            String allWordCountsDir = null;
            if (proteinID == null) {
                allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir;
            }
            else {
                //allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir + File.separator + proteinID;
                allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir + File.separator + label + "." + proteinID;
            }



            String probeMapFilePrefix = FileTools.stripPath((String)initFileCache.get("SEQUENCE", Integer.toString(experimentNumber)));
            String allKmerMatrixObjectFile = allWordCountsDir + File.separator + probeMapFilePrefix + ".all." + kmerLength + "mers.ser";
            String allKmerMatrixTableFile = allWordCountsDir + File.separator + probeMapFilePrefix + ".all." + kmerLength + "mers.tsv";

            if (initFileCache.get("DREAM", "Label") != null) {
                // do nothing here!
            }
            else if (!FileTools.exists(allKmerMatrixObjectFile)) {
                out.println("\nDid not find Kmer Matrix binary object file "+allKmerMatrixObjectFile+" for all ("+numSeqIDs+") "+kmerLength+"mers.");
            }
            else {
                out.println("\nFound existing Kmer Matrix binary object file "+allKmerMatrixObjectFile+" for all ("+numSeqIDs+") "+kmerLength+"mers. Not creating a new one.");
                return(null);
            }

            // init a new KmerMatrix object
            aKmerMatrix = new KmerMatrix(
                aSymbolMatrix.getNumSeqIDs(),
                kmerLength, aSymbolMatrix.getSeqLength(),
                alphabetName,
                alphabet,
                symbolTokenization);

            FileTools.makeDir(allWordCountsDir);

            out.println("\nCreating the Kmer Matrix for all "+kmerLength+"mers in all "+numSeqIDs+" probe-associated sequences.");

            for (int seqIdIndex=0; seqIdIndex < numSeqIDs; seqIdIndex++) {
                SymbolList probeSeq = aSymbolMatrix.getSymbolList(seqIdIndex);
                aKmerMatrix.addSymbolLists(seqIdIndex, probeSeq, kmerLength, false);

                if((seqIdIndex % 10000) == 0 ) {
                    out.println("\n\tParsed "+seqIdIndex+" out of "+numSeqIDs+" probe-associated sequences so far...");
                }
            }
            out.println(" Done.");

            // Save the KmerMatrix object file and table file
            out.print("\nSaving the Kmer Matrix binary object file for all ("+numSeqIDs+") "+kmerLength+"mers to "+allKmerMatrixObjectFile+"...");
            aKmerMatrix.writeSerializedFile(allKmerMatrixObjectFile);
            out.println(" Done.");


//             out.print("\nSaving the Kmer Matrix TSV-table file for all ("+numSeqIDs+") "+kmerLength+"mers to "+allKmerMatrixTableFile+"...");
//             aKmerMatrix.writeTable(allKmerMatrixTableFile);
//             out.println(" Done.");
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return(aKmerMatrix);

    }


    public SymbolMatrix createSymbolMatrixFile(int experimentNumber) {
        return(createSymbolMatrixFile(getAllProbeIDs(experimentNumber), experimentNumber));
    }

    public SymbolMatrix createSymbolMatrixFile(Set<String> probeSet, int experimentNumber) {

        SymbolMatrix aSymbolMatrix = null;
        try {
            String probeLocsDir = (String)initFileCache.get("DATA", "ProbeLocsDir");
            String outputDir = (String)initFileCache.get("OUTPUT", "Directory");
            String wordCountsDir = (String)initFileCache.get("OUTPUT", "WordCountsDir");

            // The ALL probeIDs file exists in the ProbeLocsDir/wordCountsDir
            //String allWordCountsDir = probeLocsDir + File.separator + wordCountsDir;
            //String allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir;

            String proteinID = (String)initFileCache.get("DREAM", "ID");
            String label = (String)initFileCache.get("DREAM", "Label");

            String allWordCountsDir = null;
            if (proteinID == null) {
                allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir;
            }
            else {
                //allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir + File.separator + proteinID;
                allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir + File.separator + label + "." + proteinID;
            }

            String probeMapFilePrefix = FileTools.stripPath((String)initFileCache.get("SEQUENCE", Integer.toString(experimentNumber)));
            String allSymbolMatrixObjectFile = allWordCountsDir + File.separator + probeMapFilePrefix + ".all.ser";
            String allSymbolMatrixTableFile = allWordCountsDir + File.separator + probeMapFilePrefix + ".all.tsv";

            String bandWidthString = (String)initFileCache.get("SEQUENCE", "BandWidth");
            int bandWidth = Integer.parseInt(bandWidthString);

            if (initFileCache.get("DREAM", "Label") != null) {
                // do nothing here!
            }
            else if (!FileTools.exists(allSymbolMatrixObjectFile)) {
                out.println("\nDid not find Symbol Matrix binary object file "+allSymbolMatrixObjectFile+" for all ("+probeSet.size()+") sequences.");
            }
            else {
                out.println("\nFound existing Symbol Matrix binary object file "+allSymbolMatrixObjectFile+" for all ("+probeSet.size()+") sequences. Not creating a new one.");
                return(null);
            }

            // init a new SymbolMatrix object
            aSymbolMatrix = new SymbolMatrix(probeSet.toArray(new String[0]), alphabetName, alphabet, symbolTokenization, packing);

            FileTools.makeDir(allWordCountsDir);

            String seqInputType = (String)initFileCache.get("SEQUENCE", "InputFileType");
            String seqKeySeparator = (String)initFileCache.get("SEQUENCE", "KeySeparator");
            String probeLocation = (String)initFileCache.get("SEQUENCE", "ProbeLocation");
            String probeWidthString = (String)initFileCache.get("SEQUENCE", "ProbeWidth");
            int probeWidth = Integer.parseInt(probeWidthString);
            boolean includeProbeSeq = (boolean)initFileCache.getBoolean("SEQUENCE","IncludeProbeSequence");
            SeqInput seqInput;
            ProbeLoc probeLoc;

            int varRegionStart = 0;
            int varRegionLength = initFileCache.getInt("SEQUENCE", "SeqLength");
            if (varRegionLength <= 0) {
                varRegionStart = -1;
            }

            out.println("\nCreating the Symbol Matrix for all "+probeSet.size()+" probe-associated sequences.");

            SequenceDB sequenceDB = null;
            LinkedHashMap<String, String> seqStringLHM = null;

            // Seqs are in a fasta file
            if (initFileCache.get("DREAM", "Label") != null) {
                seqInput = SeqInput.UNIFIED;
                probeLoc = null;
                String idString = (String)initFileCache.get("DREAM", "ID");
                // Unified files must have a column labels header line
                String unifiedFile = (String)initFileCache.get("SEQUENCE", Integer.toString(experimentNumber));
                //String unifiedPathName = probeLocsDir + File.separator + unifiedFile;
                String unifiedPathName = unifiedFile;

                out.print("\nCaching all Probe Sequences for protein "+idString+" from "+unifiedPathName+" ...");

                seqStringLHM = FileTools.readStringStringLHM(
                    unifiedPathName,
                    "\\t", // delimiter
                    1, // numLinesSkip
                    -1, // keySubStringStart
                    -1, // keySubStringLength
                    varRegionStart, // valueSubStringStart
                    varRegionLength, // valueSubStringLength
                    // -1, // valueSubStringStart
                    // -1, // valueSubStringLength
                    initFileCache.getInt("DREAM", "SequenceColumn"), // keyColumn
                    initFileCache.getInt("DREAM", "SequenceColumn"), // valueColumn
                    initFileCache.getInt("DREAM", "LabelColumn"), // idColumn
                    idString); //idString

                //seqStringLHM = FileTools.readStringStringLHM(unifiedPathName, "\t", 1);
                out.println(" Done.");
            }
            else if (seqInputType.equalsIgnoreCase(".fasta")) {
                seqInput = SeqInput.FASTA;
                probeLoc = null;

                String seqFastaFile = (String)initFileCache.get("SEQUENCE", Integer.toString(experimentNumber)) + ".fasta";
                // String seqFastaPathName = probeLocsDir + File.separator + seqFastaFile;
                String seqFastaPathName = seqFastaFile;
                //sequenceDB = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(probeLocsDir + File.separator + seqFastaFile)));
                out.print("\nCaching all ProbeMap Sequences from "+seqFastaPathName+" ...");
                sequenceDB = SeqIOTools.readFasta(new FileInputStream(seqFastaPathName), this.alphabet);
                out.println(" Done.");
            }
            // seqs are in a <label> <seq> file
            else if (seqInputType.equalsIgnoreCase(".table")) {
                seqInput = SeqInput.TABLE;
                probeLoc = null;

                // Table files must have a column labels header line
                String tableFile = (String)initFileCache.get("SEQUENCE", Integer.toString(experimentNumber));
                // String tablePathName = probeLocsDir + File.separator + tableFile;
                String tablePathName = tableFile;
                out.print("\nCaching all Probe Sequences from "+tablePathName+" ...");
                //seqStringLHM = FileTools.readStringStringLHM(tablePathName, "\t", 1, -1, -1, 0, varRegionLength, 0, 1);
                seqStringLHM = FileTools.readStringStringLHM(tablePathName, "\t", 1, -1, -1, 0, 1);
                out.println(" Done.");
            }
            // seqs are in a <seq> <intensity> file
            else if (seqInputType.equalsIgnoreCase(".unified")) {
                seqInput = SeqInput.UNIFIED;
                probeLoc = null;

                // Unified files must have a column labels header line
                String unifiedFile = (String)initFileCache.get("SEQUENCE", Integer.toString(experimentNumber));
                // String unifiedPathName = probeLocsDir + File.separator + unifiedFile;
                String unifiedPathName = unifiedFile;
                out.print("\nCaching all Probe Sequences from "+unifiedPathName+" ...");
                //seqStringLHM = FileTools.readStringStringLHM(unifiedPathName, "\t", 1, -1, -1, 0, varRegionLength, 0, 1);
                // seqStringLHM = FileTools.readStringStringLHM(unifiedPathName, "\t", 1, -1, -1, 0, 0);
                seqStringLHM = FileTools.readStringStringLHM(
                    unifiedPathName, // aFilePathName
                    "\t",            // delimiter
                    1,               // numLinesSkip
                    -1,              // keySubStringStart
                    -1,              // keySubStringLength
                    varRegionStart,  // valueSubStringStart
                    varRegionLength, // valueSubStringLength
                    // -1,              // valueSubStringStart
                    // -1,              // valueSubStringLength
                    0,               // keyColumn
                    0);              // valueColumn
                out.println(" Done.");
            }
            // seqs are in a <label> <seq> <intensity> file
            else if (seqInputType.equalsIgnoreCase(".unifiedLabels")) {
                seqInput = SeqInput.UNIFIED;
                probeLoc = null;

                // Unified files must have a column labels header line
                String unifiedFile = (String)initFileCache.get("SEQUENCE", Integer.toString(experimentNumber));
                // String unifiedPathName = probeLocsDir + File.separator + unifiedFile;
                String unifiedPathName = unifiedFile;
                out.print("\nCaching all Probe Sequences from "+unifiedPathName+" ...");
                seqStringLHM = FileTools.readStringStringLHM(
                    unifiedPathName, // aFilePathName
                    "\t",            // delimiter
                    1,               // numLinesSkip
                    -1,              // keySubStringStart
                    -1,              // keySubStringLength
                    varRegionStart,  // valueSubStringStart
                    varRegionLength, // valueSubStringLength
                    // -1,              // valueSubStringStart
                    // -1,              // valueSubStringLength
                    0,               // keyColumn
                    1);              // valueColumn

                //seqStringLHM = FileTools.readStringStringLHM(unifiedPathName, "\t", 1);
                out.println(" Done.");
            }
            else {
                out.println("\nParameters: probeLocation="+probeLocation+" bandWidth="+bandWidthString+" probeWidth="+probeWidthString+" includeProbeSequence="+Boolean.toString(includeProbeSeq)+".");

                seqInput = SeqInput.UCSCGENOME;
                if (probeLocation.equalsIgnoreCase("Center")) {
                    probeLoc = ProbeLoc.CENTER;
                }
                else if (probeLocation.equalsIgnoreCase("Left")) {
                    probeLoc = ProbeLoc.LEFT;
                }
                else if (probeLocation.equalsIgnoreCase("Right")) {
                    probeLoc = ProbeLoc.RIGHT;
                }
                else { //Error
                    probeLoc = null;
                    out.println("\nError : could not get sequence because of unsupported [SEQUENCE]->ProbeLocation "+probeLocation+".");
                    return(null);
                }
            }

            int seqIdIndex = 0;
            for (String seqLabelString : probeSet) {

                // Get the sequence from *.fasta or UCSCGenome
                ArrayList<Sequence> probeSeqs = null;
                if (seqInput == SeqInput.FASTA) {
                    probeSeqs = this.getSequences(seqLabelString, sequenceDB);
                }
                else if (seqInput == SeqInput.TABLE) {
                    probeSeqs = this.getSequences(seqLabelString, seqStringLHM);
                }
                else if (seqInput == SeqInput.UNIFIED) {
                    // String[] seqLabelStringArray = seqLabelString.split(seqKeySeparator);
                    // out.println("\nseqLabelStringArray[0]="+seqLabelStringArray[0]+"\n");
                    // probeSeqs = this.getSequences(seqLabelStringArray[0], seqLabelStringArray[1], probeLoc, includeProbeSeq, bandWidth, probeWidth);
                    // out.println("\nseqLabelString="+seqLabelString+"\n");
                    probeSeqs = this.getSequences(seqLabelString, seqStringLHM);
                }
                else {
                    String[] seqLabelStringArray = seqLabelString.split(seqKeySeparator);
                    probeSeqs = this.getSequences(seqLabelStringArray[0], seqLabelStringArray[1], probeLoc, includeProbeSeq, bandWidth, probeWidth);
                    //probeSeq.setName(seqLabelString);
                }

                // add sequence to cache
                //out.print("\nCreating Symbol Matrix file for probe-associated sequence "+seqLabelString+"...");
                for (Sequence probeSeq : probeSeqs) {
                    // out.println("\nProbe sequence for "+seqLabelString+" is "+probeSeq.seqString());
                    aSymbolMatrix.addSymbolList(seqIdIndex, probeSeq);
                }
                //out.println(" Done.");

                seqIdIndex++;
                if((seqIdIndex % 10000) == 0 ) {
                    out.println("\n\tParsed "+seqIdIndex+" out of "+probeSet.size()+" probe-associated sequences so far...");
                }

                //return(null);
            }
            out.println(" Done.");

            // Save the SymbolMatrix object file and table file
            out.print("\nSaving the Symbol Matrix binary object file for all ("+probeSet.size()+") sequences to "+allSymbolMatrixObjectFile+"...");
            aSymbolMatrix.writeSerializedFile(allSymbolMatrixObjectFile);
            out.println(" Done.");


            out.print("\nSaving the Symbol Matrix TSV-table file for all ("+probeSet.size()+") sequences to "+allSymbolMatrixTableFile+"...");
            aSymbolMatrix.writeTable(allSymbolMatrixTableFile);
            out.println(" Done.");
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return(aSymbolMatrix);

    }

    // caches the kmersByProbe KmerMatrix File for each kmer-length
    //
    // The ALL probeIDs file exists in the ProbeLocsDir/wordCountsDir
    // The USE probeIDs file exists in the ProbeIntensitiesDir/outputDir/wordCountsDir
    //
    public KmerMatrix cacheKmerMatrixFile(int kmerLength, int experimentNumber, int unifiedColumn) {

        KmerMatrix aKmerMatrix = null;
        try {

//             if ((aKmerMatrix = getKmerMatrix(kmerLength, experimentNumber)) != null) {
//                 return(aKmerMatrix);
//             }

            //out.println("Here1!");

            String intensitiesDir = (String)initFileCache.get("DATA", "ProbeIntensitiesDir");
            String probeLocsDir = (String)initFileCache.get("DATA", "ProbeLocsDir");
            String outputDir = (String)initFileCache.get("OUTPUT", "Directory");
            String wordCountsDir = (String)initFileCache.get("OUTPUT", "WordCountsDir");

            // The ALL probeIDs file exists in the ProbeLocsDir/wordCountsDir
            //String allWordCountsDir = probeLocsDir + File.separator + wordCountsDir;
            //String allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir;

            String proteinID = (String)initFileCache.get("DREAM", "ID");
            String label = (String)initFileCache.get("DREAM", "Label");

            String allWordCountsDir = null;
            if (proteinID == null) {
                allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir;
            }
            else {
                //allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir + File.separator + proteinID;
                allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir + File.separator + label + "." + proteinID;
            }

            // The USE probeIDs file exists in the ProbeIntensitiesDir/outputDir/wordCountsDir
            //String useWordCountsDir = intensitiesDir + File.separator + outputDir + File.separator + wordCountsDir;

            String probeMapFilePrefix = FileTools.stripPath((String)initFileCache.get("SEQUENCE", Integer.toString(experimentNumber)));
            String allKmerCountsObjectFile = allWordCountsDir + File.separator + probeMapFilePrefix + ".all." + kmerLength + "mers.ser";
            String allKmerCountsTableFile = allWordCountsDir + File.separator + probeMapFilePrefix + ".all." + kmerLength + "mers.tsv";
            //String useKmerCountsObjectFile = useWordCountsDir + File.separator + probeMapFilePrefix + ".use." + kmerLength + "mers.ser";
            //String useKmerCountsTableFile = useWordCountsDir + File.separator + probeMapFilePrefix + ".use." + kmerLength + "mers.tsv";

            SymbolMatrix symbolMatrix = getSymbolMatrix(experimentNumber, unifiedColumn);

            // This returns null and does nothing if counts-file already exists
            if (initFileCache.get("DREAM", "Label") == null) {
                aKmerMatrix = createKmerMatrixFile(kmerLength, symbolMatrix, experimentNumber);
                //out.println("HI!!!!!");
            }

            if (aKmerMatrix == null) {
                out.print("\nCaching Kmer Matrix binary object file for all "+kmerLength+"mers from "+allKmerCountsObjectFile+" ...");
                aKmerMatrix = (KmerMatrix)FileTools.readSerializedFile(allKmerCountsObjectFile);
                out.println(" Done.");
            }

            // get the Ordered Subset of word counts
            out.print("\nSorting the Kmer Matrix for all "+kmerLength+"mers by intensity-value ...");
            aKmerMatrix.setOrderedSubset(
                symbolMatrix.getOriginalSeqIDs(),
                symbolMatrix.getSeqIDs());
            //getUsedProbeIDs(experimentNumber).toArray(new String[0]));
            out.println(" Done.");


//             if (!FileTools.exists(useKmerCountsObjectFile)) {
//                 out.println("\nDid not find the intensity-ordered Kmer Matrix binary object file "+useKmerCountsObjectFile+" for all "+kmerLength+"mers in all "+numberOfUsedProbes(experimentNumber)+" intensity-ordered, probe-associated sequences.");

//                 if (aKmerMatrix == null) {
//                     out.print("\nCaching Kmer Matrix binary object file for all "+kmerLength+"mers from "+allKmerCountsObjectFile+" ...");
//                     aKmerMatrix = new KmerMatrix(alphabetName, alphabet, symbolTokenization, allKmerCountsObjectFile);
//                     out.println(" Done.");
//                 }

//                 // get the Ordered Subset of word counts
//                 out.print("\nSorting the Kmer Matrix for all "+kmerLength+"mers by intensity-value ...");
//                 aKmerMatrix.setOrderedSubset(
//                     symbolMatrix.getOriginalSeqIDs(),
//                     symbolMatrix.getSeqIDs());
//                 //getUsedProbeIDs(experimentNumber).toArray(new String[0]));
//                 out.println(" Done.");

//                 // save used ordered subset
//                 FileTools.makeDir(useWordCountsDir);

//                 out.print("\nSaving the intensity-ordered Kmer Matrix binary object file for all "+kmerLength+"mers in all "+numberOfUsedProbes(experimentNumber)+" intensity-ordered, probe-associated sequences to "+useKmerCountsObjectFile+" ...");
//                 aKmerMatrix.writeSerializedFile(useKmerCountsObjectFile);
//                 out.println(" Done.");

// //                 out.print("\nSaving the intensity-ordered Kmer Matrix TSV-table file for all "+kmerLength+"mers in all "+numberOfUsedProbes(experimentNumber)+" intensity-ordered, probe-associated sequences to "+useKmerCountsTableFile+" ...");
// //                 aKmerMatrix.writeTable(useKmerCountsTableFile);
// //                 out.println(" Done.");

//             }
//             else {
//                 out.print("\nCaching the intensity-ordered Kmer Matrix binary object file for all "+kmerLength+"mers in all "+numberOfUsedProbes(experimentNumber)+" intensity-ordered, probe-associated sequences from "+useKmerCountsObjectFile+" ...");
//                 aKmerMatrix = new KmerMatrix(alphabetName, alphabet, symbolTokenization, useKmerCountsObjectFile);
//                 out.println(" Done.");
//             }


//             // If .ser object file exists then load it
//             if (FileTools.exists(kmerCountsObjectFile)) {

//                 out.print("\nCaching Kmer Matrix binary object file for all "+kmerLength+"mers in all "+numberOfUsedProbes(experimentNumber)+" probe-associated sequences from "+kmerCountsObjectFile+" ...");
//                 aKmerMatrix = new KmerMatrix(alphabetName, alphabet, symbolTokenization, kmerCountsObjectFile);
//                 out.println(" Done.");
//                 //aKmerMatrix.writeTable(kmerCountsTableFile);

//             }
//             else if (FileTools.exists(kmerCountsTableFile)) {

//                 out.print("\nCaching Kmer Matrix TSV-table file for all "+kmerLength+"mers in all "+numberOfUsedProbes(experimentNumber)+" probe-associated sequences from "+kmerCountsTableFile+" ...");
//                 aKmerMatrix = new KmerMatrix(numberOfUsedProbes(experimentNumber), kmerLength, bandWidth, alphabetName, alphabet, symbolTokenization, kmerCountsTableFile);
//                 out.println(" Done.");

//             }

//             // else create the WordCount objects from *.fasta or UCSCGenome and save them
//             else {

//                 out.println("\nDid not find Kmer Matrix binary object file "+kmerCountsObjectFile+" nor TSV-table file "+kmerCountsTableFile+" for all "+kmerLength+"mers.");
//                 aKmerMatrix = createKmerMatrixFile(kmerLength, getUsedProbeIDs(experimentNumber), experimentNumber);
//             }

        }
        catch (Exception e) {
            e.printStackTrace();
        }

        putKmerMatrix(kmerLength, experimentNumber, unifiedColumn, aKmerMatrix);
        return(aKmerMatrix);
    }


    public SymbolMatrix loadAllSymbolMatrixFile(int experNum) {
        String probeLocsDir = (String)initFileCache.get("DATA", "ProbeLocsDir");
        String wordCountsDir = (String)initFileCache.get("OUTPUT", "WordCountsDir");
        String outputDir = (String)initFileCache.get("OUTPUT", "Directory");

        //String allWordCountsDir = probeLocsDir + File.separator + wordCountsDir;
        //String allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir;

        String proteinID = (String)initFileCache.get("DREAM", "ID");
        String label = (String)initFileCache.get("DREAM", "Label");

        String allWordCountsDir = null;
        if (proteinID == null) {
            allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir;
        }
        else {
            //allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir + File.separator + proteinID;
            allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir + File.separator + label + "." + proteinID;
        }

        String probeMapFilePrefix = FileTools.stripPath((String)initFileCache.get("SEQUENCE", Integer.toString(experNum)));
        String allSymbolMatrixObjectFile = allWordCountsDir + File.separator + probeMapFilePrefix + ".all.ser";

        out.print("\nLoading Symbol Matrix binary object file from "+allSymbolMatrixObjectFile+" ...");
        SymbolMatrix symbolMatrix = new SymbolMatrix(alphabetName, alphabet, symbolTokenization, allSymbolMatrixObjectFile, packing);
        out.println(" Done.");

        return(symbolMatrix);
    }

    public double[] getOrderedInitPoolFreqs(int experimentNumber, int unifiedColumn) {
      String initPoolExpectedFreqsFile = (String)initFileCache.get("INTENSITIES", "InitPoolExpectedFreqsFile");
      LinkedHashMap<String, Double> initPoolFreqs = FileTools.readStringDoubleLHM(initPoolExpectedFreqsFile, "\t", 1);

      // order the initPoolFreqs by the Intensities
      Set<String> usedProbeIDs = getUsedProbeIDs(experimentNumber, unifiedColumn);

      double[] orderedInitPoolFreqs = new double[usedProbeIDs.size()];
      int i = 0;
      //for (int i=0; i < usedProbeIDs.size(); i++) {
      for (String probeID : usedProbeIDs) {
          //orderedInitPoolFreqs[i] = initPoolFreqs.get(usedProbeIDs.get(i));
          orderedInitPoolFreqs[i] = initPoolFreqs.get(probeID);
          i++;
      }

      return(orderedInitPoolFreqs);
    }


    // caches the SymbolListBySeqID SymbolMatrix File
    //
    // The ALL probeIDs file exists in the ProbeLocsDir/wordCountsDir
    // The USE probeIDs file exists in the ProbeIntensitiesDir/outputDir/wordCountsDir
    //
    public SymbolMatrix cacheSymbolMatrixFile(int experimentNumber, int unifiedColumn) {

        SymbolMatrix aSymbolMatrix = null;
        try {

//             if ((aSymbolMatrix = getSymbolMatrix(experimentNumber)) != null) {
//                 return(aSymbolMatrix);
//             }

            //out.println("Here1!");

            String intensitiesDir = (String)initFileCache.get("DATA", "ProbeIntensitiesDir");
            String probeLocsDir = (String)initFileCache.get("DATA", "ProbeLocsDir");
            String outputDir = (String)initFileCache.get("OUTPUT", "Directory");
            String wordCountsDir = (String)initFileCache.get("OUTPUT", "WordCountsDir");

            // The ALL probeIDs file exists in the ProbeLocsDir/wordCountsDir
            //String allWordCountsDir = probeLocsDir + File.separator + wordCountsDir;
            //String allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir;

            String proteinID = (String)initFileCache.get("DREAM", "ID");
            String label = (String)initFileCache.get("DREAM", "Label");

            String allWordCountsDir = null;
            if (proteinID == null) {
                allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir;
            }
            else {
                //allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir + File.separator + proteinID;
                allWordCountsDir = probeLocsDir + File.separator + outputDir + File.separator + wordCountsDir + File.separator + label + "." + proteinID;
            }


            // The USE probeIDs file exists in the ProbeIntensitiesDir/outputDir/wordCountsDir
            String useWordCountsDir = intensitiesDir + File.separator + outputDir + File.separator + wordCountsDir;

            String probeMapFilePrefix = FileTools.stripPath((String)initFileCache.get("SEQUENCE", Integer.toString(experimentNumber)));
            String allSymbolMatrixObjectFile = allWordCountsDir + File.separator + probeMapFilePrefix + ".all.ser";
            String allSymbolMatrixTableFile = allWordCountsDir + File.separator + probeMapFilePrefix + ".all.tsv";

            //String useSymbolMatrixObjectFile = useWordCountsDir + File.separator + probeMapFilePrefix + ".use.ser";
            //String useSymbolMatrixTableFile = useWordCountsDir + File.separator + probeMapFilePrefix + ".use.tsv";

            // This returns null and does nothing if counts-file already exists
            //aSymbolMatrix = createSymbolMatrixFile(getUsedProbeIDs(experimentNumber), experimentNumber);
            if (initFileCache.get("DREAM", "Label") == null) {
                aSymbolMatrix = createSymbolMatrixFile(getAllProbeIDs(experimentNumber), experimentNumber);
            }

            if (aSymbolMatrix == null) {
                out.print("\nCaching Symbol Matrix binary object file from "+allSymbolMatrixObjectFile+" ...");
                aSymbolMatrix = new SymbolMatrix(alphabetName, alphabet, symbolTokenization, allSymbolMatrixObjectFile, packing);
                out.println(" Done.");
            }

            // get the Ordered Subset of word counts
            //            if (unifiedColumn >= 0) {
                out.print("\nSorting the Symbol Matrix by intensity-value ...");
                aSymbolMatrix.setOrderedSubset(getUsedProbeIDs(experimentNumber, unifiedColumn).toArray(new String[0]));
                out.println(" Done.");
                System.out.flush();
                //            }

//             if (!FileTools.exists(useSymbolMatrixObjectFile)) {
//                 out.println("\nDid not find the intensity-ordered Symbol Matrix binary object file "+useSymbolMatrixObjectFile+" for all "+numberOfUsedProbes(experimentNumber)+" intensity-ordered, probe-associated sequences.");

//                 if (aSymbolMatrix == null) {
//                     out.print("\nCaching Symbol Matrix binary object file from "+allSymbolMatrixObjectFile+" ...");
//                     aSymbolMatrix = new SymbolMatrix(alphabetName, alphabet, symbolTokenization, allSymbolMatrixObjectFile, packing);
//                     out.println(" Done.");
//                 }

//                 // get the Ordered Subset of word counts
//                 out.print("\nSorting the Symbol Matrix by intensity-value ...");
//                 aSymbolMatrix.setOrderedSubset(getUsedProbeIDs(experimentNumber).toArray(new String[0]));
//                 out.println(" Done.");

//                 // save used ordered subset
//                 FileTools.makeDir(useWordCountsDir);

//                 out.print("\nSaving the intensity-ordered Symbol Matrix binary object file for all "+numberOfUsedProbes(experimentNumber)+" intensity-ordered, probe-associated sequences to "+useSymbolMatrixObjectFile+" ...");
//                 aSymbolMatrix.writeSerializedFile(useSymbolMatrixObjectFile);
//                 out.println(" Done.");

//                 out.print("\nSaving the intensity-ordered Symbol Matrix TSV-table file for all "+numberOfUsedProbes(experimentNumber)+" intensity-ordered, probe-associated sequences to "+useSymbolMatrixTableFile+" ...");
//                 aSymbolMatrix.writeTable(useSymbolMatrixTableFile);
//                 out.println(" Done.");

//             }
//             else {
//                 out.print("\nCaching the intensity-ordered Symbol Matrix binary object file for all "+numberOfUsedProbes(experimentNumber)+" intensity-ordered, probe-associated sequences from "+useSymbolMatrixObjectFile+" ...");
//                 aSymbolMatrix = new SymbolMatrix(alphabetName, alphabet, symbolTokenization, useSymbolMatrixObjectFile, packing);
//                 out.println(" Done.");
//             }


//             // If .ser object file exists then load it
//             if (FileTools.exists(kmerCountsObjectFile)) {

//                 out.print("\nCaching Symbol Matrix binary object file for all "+numberOfUsedProbes(experimentNumber)+" probe-associated sequences from "+kmerCountsObjectFile+" ...");
//                 aSymbolMatrix = new SymbolMatrix(alphabetName, alphabet, symbolTokenization, kmerCountsObjectFile);
//                 out.println(" Done.");
//                 //aSymbolMatrix.writeTable(kmerCountsTableFile);

//             }
//             else if (FileTools.exists(kmerCountsTableFile)) {

//                 out.print("\nCaching Symbol Matrix TSV-table file for all "+numberOfUsedProbes(experimentNumber)+" probe-associated sequences from "+kmerCountsTableFile+" ...");
//                 aSymbolMatrix = new SymbolMatrix(numberOfUsedProbes(experimentNumber), kmerLength, bandWidth, alphabetName, alphabet, symbolTokenization, kmerCountsTableFile);
//                 out.println(" Done.");

//             }

//             // else create the WordCount objects from *.fasta or UCSCGenome and save them
//             else {

//                 out.println("\nDid not find Symbol Matrix binary object file "+kmerCountsObjectFile+" nor TSV-table file "+kmerCountsTableFile+".");
//                 aSymbolMatrix = createSymbolMatrixFile(getUsedProbeIDs(experimentNumber), experimentNumber);
//             }

        }
        catch (Exception e) {
            e.printStackTrace();
        }

        putSymbolMatrix(experimentNumber, unifiedColumn, aSymbolMatrix);
        return(aSymbolMatrix);
    }

    // Loads a *.cel file into a linked list of intensities
    //
    // Input: *.cel files
    // Output: LinkedHashMap<"X,Y", intensity>
    public LinkedHashMap<String, Double> getCelIntensitiesCache(String dirString, String inputFilePrefixString, String keySeparator) {

        LinkedHashMap<String, Double> intensitiesCache = null;

        try {
            String inputFullPathName;
            String lineString;
            String[] lineFields;
            BufferedReader inputBuffer;

            String celFileString1 = dirString+File.separator+inputFilePrefixString+".CEL";
            String celFileString2 = dirString+File.separator+inputFilePrefixString+".cel";

            int maxColumn = -1;

            // *.CEL file = X	Y	MEAN	STDV	NPIXELS
            intensXColumn = 0;
            intensYColumn = 1;
            intensIntensityColumn = 2;
            maxColumn = 2;

            if (FileTools.exists(celFileString1)) {
                inputFullPathName = celFileString1;
            }
            else if (FileTools.exists(celFileString2)) {
                inputFullPathName = celFileString2;
            }
            else { // Error!
                out.println("\nError : could not cache intensity values because couldn't find either "+celFileString1+" or "+celFileString2+".");
                return(null);
            }

            int celFileLineCount = FileTools.getLineCount(inputFullPathName);
            //int celFileLineCount = FileTools.getLineCount(inputFullPathName,"#", true);

            inputBuffer = new BufferedReader(new InputStreamReader(new FileInputStream(inputFullPathName)));

            out.print("\nCreating probe-intensity values cache from "+inputFullPathName+" ...");

            // Go thru the the cel file until you hit "CellHeader=" at the begining of the line
            int headerLineCount = 0;
            while((lineString = inputBuffer.readLine()) != null) {
                headerLineCount++;
                if (lineString.trim().startsWith("CellHeader=")) {
                    break;
                }
            }

            // initialize the intensities cache with approximate size
            intensitiesCache = new LinkedHashMap<String, Double>(celFileLineCount - headerLineCount);

            // read all the probe intensity values in the [INTENSITY] section
            while((lineString = inputBuffer.readLine()) != null) {

                // if line starts with "[" then we have entered a new section so we are done
                if (lineString.trim().startsWith("[")) {
                    break;
                }

                // include this line if line doesn't start with commentToken "#"
                if (!lineString.trim().equals("") && !lineString.trim().startsWith("#")) {
                    lineFields = tabsPattern.split(lineString);

                    //Only cache the intensity value if it exists
                    if (maxColumn < lineFields.length) {
                        Double probeIntensity = new Double(lineFields[intensIntensityColumn]);
                        //out.println(lineFields[intensXColumn].trim() + keySeparator + lineFields[intensYColumn].trim() + "\t" + probeIntensity.toString());
                        intensitiesCache.put(lineFields[intensXColumn].trim() + keySeparator + lineFields[intensYColumn].trim(), probeIntensity);
                    }
                }
            }

            out.println(" Done.");

            inputBuffer.close();
        }
        catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }

        return(intensitiesCache);
    }

    enum NormalizeAlgorithm {RATIO, LOGRATIO, DIFF}
    //
    // performs normalization on replicates of test and control intensities
    // the result is one set of averaged intensities
    //
    // Intensity files must have a column labels header line
    public LinkedHashMap<String, Double> getAverageXYIntensities(String dirString, String[] inputFilePrefixStrings, String grouping, String keySeparator, String intensInputFileType) {
        int testFileNum = -1;
        int controlFileNum = -1;

        String algorithmString = (String)initFileCache.get("NORMALIZATION", "Algorithm");

        NormalizeAlgorithm algorithm;
        if (algorithmString.equalsIgnoreCase("ratio")) {
            out.println("\nNormalization of probe intensities is the ratio avg(test)/avg(ctrl).");
            algorithm = NormalizeAlgorithm.RATIO;
        }
        else if (algorithmString.equalsIgnoreCase("logratio")) {
            out.println("\nNormalization of probe intensities is the log-ratio log(avg(test)/avg(ctrl)).");
            algorithm = NormalizeAlgorithm.LOGRATIO;
        }
        else if (algorithmString.equalsIgnoreCase("diff")) {
            out.println("\nNormalization of probe intensities is the difference avg(test) - avg(ctrl).");
            algorithm = NormalizeAlgorithm.DIFF;
        }
        else { //Error
            out.println("\nError : could not normalize probe intensities because of unsupported [NORMALIZATION]->Algorithm "+algorithmString+".");
            return(null);
        }

        // grouping is a string of 1's (test) and 0's (control). (e.g. - 111000)
        LinkedHashSet<LinkedHashMap> testIntensitiesCaches = new LinkedHashSet<LinkedHashMap>(grouping.length());
        LinkedHashSet<LinkedHashMap> controlIntensitiesCaches = new LinkedHashSet<LinkedHashMap>(grouping.length());

        LinkedHashMap<String, Double> normalizedCache = null;
        LinkedHashMap<String, Double> anIntensitiesCache = null;

        //Cache all the test and control Intensity Files
        for (int i=0; i < grouping.length(); i++) {

            char groupID = grouping.charAt(i);

            if (intensInputFileType.equalsIgnoreCase(".cel")) {
                anIntensitiesCache = getCelIntensitiesCache(dirString, inputFilePrefixStrings[i], keySeparator);
            }
            else {
                anIntensitiesCache = FileTools.readStringDoubleLHM(dirString+File.separator+inputFilePrefixStrings[i], "\t", 1);
            }

            if (groupID == '1') {
                testIntensitiesCaches.add(anIntensitiesCache);
            }
            else {
                controlIntensitiesCaches.add(anIntensitiesCache);
            }
        }

        // Go through the test Caches and create average
        Double value = null;
        List<Double> replicates = null;

        // Grab the first test intensity cache (there has to be at least one)
        anIntensitiesCache = (LinkedHashMap)testIntensitiesCaches.iterator().next();

        // initialize the normalized intensities cache with approximate size
        normalizedCache = new LinkedHashMap<String, Double>(anIntensitiesCache.size());

        // average the testCache values
        for (String key : anIntensitiesCache.keySet()) {
            replicates = new ArrayList<Double>();
            Iterator testCachesIter = testIntensitiesCaches.iterator();
            while (testCachesIter.hasNext()) {
                LinkedHashMap<String, Double> aTestCache = (LinkedHashMap<String, Double>)testCachesIter.next();
                replicates.add(aTestCache.get(key));
            }
            double testAvg = MathTools.mean(replicates);
            normalizedCache.put(key, new Double(testAvg));
        }

        // Go through the control Caches and create average
        out.print("\nCreating normalized probe-intensity averages cache ...");
        for (String key : anIntensitiesCache.keySet()) {
            double testAvg = ((Double)normalizedCache.get(key)).doubleValue();
            normalizedCache.remove(key);

            replicates = new ArrayList<Double>();
            Iterator controlCachesIter = controlIntensitiesCaches.iterator();
            while (controlCachesIter.hasNext()) {
                LinkedHashMap<String, Double> aControlCache = (LinkedHashMap<String, Double>)controlCachesIter.next();
                Double replicate = aControlCache.get(key);
                if (replicate != null) {
                    replicates.add(replicate);
                }
            }

            // if there are no control intensities then skip this key
            if (replicates.size() == 0) {
                continue;
            }

            double controlAvg = MathTools.mean(replicates);
            double normalizedIntensity = -1;

            switch (algorithm) {
            case RATIO:
                normalizedIntensity = testAvg / controlAvg;
                break;
            case LOGRATIO:
                normalizedIntensity = Math.log(testAvg / controlAvg);
                break;
            case DIFF:
                normalizedIntensity = testAvg - controlAvg;
                break;
            default: //Error
                out.println("\nError : could not normalize probe intensities because of unsupported [NORMALIZATION]->Algorithm "+algorithm.toString()+".");
                return(null);
            }

            normalizedCache.put(key, new Double(normalizedIntensity));
        }
        out.println(" Done.");

        return(normalizedCache);
    }

    // takes a *.bar.txt files and returns a map of label to intensity
    //
    // Input: *.bar.txt files
    // Output: LinkedHashMap<Label-Index, Intensities> of cached intensities
    public LinkedHashMap<String, Double> getProbeMapIntensities(String dirString, String inputFilePrefixString, String keySeparator, String outputFilePrefixString) {

        String inputFullPathName;
        LinkedHashMap<String, Double> keyToIntensityMap = null;

        try {
            String lineString;
            String[] lineFields;
            BufferedReader inputBuffer;

            String barFileString = dirString+File.separator+inputFilePrefixString+".bar.txt";

            // *.bar.txt file = Chr [tab] coordinate [tab] signal_value
            if (FileTools.exists(barFileString)) {
                inputFullPathName = barFileString;
                intensSeqLabelColumn = 0;
                intensSeqIndexColumn = 1;
                intensIntensityColumn = 2;
            }
            else { // Error!
                out.println("\nError : could not generate the ProbeMap intensities cache because couldn't find "+barFileString+".");
                return(null);
            }

            //int lineCount = FileTools.getLineCount(barFileString);
            //keyToIntensityMap = new LinkedHashMap<String, Double>(lineCount-1);
            keyToIntensityMap = new LinkedHashMap<String, Double>();

            inputBuffer = new BufferedReader(new InputStreamReader(new FileInputStream(inputFullPathName)));

            out.print("\nCreating the ProbeMap intensities cache from "+inputFullPathName+" ...");

            while((lineString = inputBuffer.readLine()) != null) {
                // include this line if line doesn't start with commentToken "#"
                if (!lineString.trim().startsWith("#")) {
                    lineFields = tabsPattern.split(lineString);

                    // put this intensity into the LinkedHashMap
                    keyToIntensityMap.put(
                        lineFields[intensSeqLabelColumn].trim() + keySeparator + lineFields[intensSeqIndexColumn].trim(),
                        Double.parseDouble(lineFields[intensIntensityColumn]));
                }
            }

            out.println(" Done.");

            inputBuffer.close();
        }
        catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }
        return(keyToIntensityMap);
    }

    // takes a *.bpmap.txt files and returns a map of label to intensity
    //
    // Input: *.bmap.txt file
    //        LinkedHashMap<XY, Intensities> of cached intensities
    // Output: LinkedHashMap<Label-Index, Intensities> of cached intensities
    public LinkedHashMap<String, Double> getProbeMapIntensities(LinkedHashMap<String, Double> anXYtoIntensityHM, String outputFilePrefixString, String mapKeySeparator, String intensKeySeparator, int experimentNumber) {

        LinkedHashMap<String, Double> keyToIntensityMap = null;

        try {
            String lineString;
            String[] lineFields;
            double probeIntensity;
            String probeMapFilePathName = this.probeMapFilePathNames[experimentNumber-1];

            //int lineCount = FileTools.getLineCount(this.probeMapFilePathNames[experimentNumber]);
            //keyToIntensityMap = new LinkedHashMap<String, Double>(lineCount-1);
            keyToIntensityMap = new LinkedHashMap<String, Double>();

            BufferedReader inputBuffer = new BufferedReader(new FileReader(probeMapFilePathName));

            out.print("\nCreating ProbeMap intensities cache from "+probeMapFilePathName+" ...");

            while((lineString = inputBuffer.readLine()) != null) {
                // include this line if line doesn't start with commentToken "#" or "PMX"
                if (!lineString.trim().startsWith("#") && !lineString.trim().startsWith("PMX")) {
                    lineFields = tabsPattern.split(lineString);
                    probeIntensity =
                        ((Double) anXYtoIntensityHM.get(lineFields[mapXColumn]+intensKeySeparator+lineFields[mapYColumn]))
                        .doubleValue();

                    if (!Double.isNaN(probeIntensity)) {
                        // put this probeIntensity into the LinkedHashMap<Label-Index, Intensities>
                        keyToIntensityMap.put(
                            lineFields[mapSeqLabelColumn]+mapKeySeparator+lineFields[mapSeqIndexColumn],
                            probeIntensity);
                    }
                }
            }

            out.println(" Done.");

            inputBuffer.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
        return(keyToIntensityMap);
    }

    // takes a *.bpmap.txt files and returns a set of labels
    //
    // Input: *.bmap.txt file
    // Output: LinkedHashSet<Label-Index> of cached ProbeIDs
    public LinkedHashSet<String> getProbeMapProbeIDs(String mapKeySeparator, int experimentNumber) {

        LinkedHashSet<String> probeIDsSet = null;

        try {
            String lineString;
            String[] lineFields;
            String probeMapFilePathName = this.probeMapFilePathNames[experimentNumber-1];

            //int lineCount = FileTools.getLineCount(this.probeMapFilePathName);
            //probeIDsSet = new LinkedHashSet<String>(lineCount-1);
            probeIDsSet = new LinkedHashSet<String>();

            BufferedReader inputBuffer = new BufferedReader(new FileReader(probeMapFilePathName));

            out.print("\nCaching all ProbeMap ProbeIDs from "+probeMapFilePathName+" ...");

            while((lineString = inputBuffer.readLine()) != null) {
                // include this line if line doesn't start with commentToken "#" or "PMX"
                if (!lineString.trim().startsWith("#") && !lineString.trim().startsWith("PMX")) {
                    //out.println("\n"+lineString);

                    lineFields = tabsPattern.split(lineString);

                    probeIDsSet.add(lineFields[mapSeqLabelColumn]+mapKeySeparator+lineFields[mapSeqIndexColumn]);
                }
            }

            out.println(" Done.");

            inputBuffer.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
        return(probeIDsSet);
    }

}
