/*
 * KmerCounts.java - Todd Riley
 *
 */



import org.biojava.bio.symbol.*;
import org.biojava.bio.dist.*;
import org.biojava.utils.*;
import org.biojava.bio.dp.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.db.*;
import org.biojava.bio.seq.io.*;
import org.biojava.bio.mydp.*;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.math.BigInteger;
import java.util.zip.GZIPInputStream;

import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.iterator.TLongIntIterator;
import gnu.trove.map.hash.TLongIntHashMap;
import gnu.trove.set.hash.TLongHashSet;
import gnu.trove.map.hash.TLongDoubleHashMap;
import gnu.trove.iterator.TLongDoubleIterator;

import org.apache.commons.lang.mutable.MutableInt;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.BinomialDistribution;


///////////////////////////////////////////////////////////////////////////////////////////////////////
// This class uses the Trove clases while encoding Kmers into Longs. Since Java uses signed longs,
// this class works for Kmer Length <= 31
//
// In Stored hashmaps, the length of the kmer keys is stored in (key = -1)
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

public class KmerCounts {

	private static String tabs = "\\t+";
	private static Pattern tabsPattern = Pattern.compile(tabs);
	FiniteAlphabet alphabet = DNATools.getDNA();
	ReversibleTranslationTable complementTable = DNATools.complementTable();

	String leftAdapterSeq = "";
	String rightAdapterSeq = "";

	public static void main(String[] args) throws IOException {
		KmerCounts builder = new KmerCounts();
		builder.parseCommandLine(args);
	}

	public void parseCommandLine(String[] args) throws IOException {

		try {

			String title = "KmerCounts";
			String version = "1.09";

			int numMandatoryArgs = 1;

			int kmerLength = 0;
			double countMin = 0;
			String combineRcString = "No";
			boolean combineRC = false;
			String multString = "Yes";
			boolean mult = true;
			String skipHeaderString = "Yes";
			boolean skipHeader = true;
			String probSkipHeaderString = "Yes";
			boolean probSkipHeader = true;

			// String overlapAdaptersString = "No";
			// boolean overlapAdapters = false;
			int overlapAdaptersLen = 0;
			int start = 0;
			int regionLength = -1;

			// String probOverlapAdaptersString = "No";
			// boolean probOverlapAdapters = false;
			int probOverlapAdaptersLen = 0;
			int probStart = 0;
			int probRegionLength = -1;

			String regexString = null;

			String sortBy = null;
			// String sortBy = "counts";

			String subKmerCountsFile = null;

			String ampKmerCountsFile = null;
			String ampSubKmerCountsFile = null;

			String ampR0KmerCountsFile = null;
			String ampR0SubKmerCountsFile = null;

			String observedCountsFile = null;
			String observedCountsFileFormat = null;

			String sortByFile = null;
			boolean expectedCounts = false;
			boolean filterObservedCounts = false;
			double filterPvalueThresh = .001;
			boolean countStats = false;

			String saveFileFormat = null;
			String saveFile = null;

			String proteinString = null;
			String labelString = null;
			String roundString = null;
			int round = -1;

			System.out.println("\n" + title	+ " v" + version + " - A Kmer-counts construction tool that reads in and parses sequences to generate kmer-tab-count tables using a sliding window approach. Uses compact bit representation and the gnu trove classes for compact, fast, and efficient hashmap memory management. Also can pass in K-mer counts tables to build an Nth Order Markov Model and calculate the probabilities of seeing the sequences.\n");
			System.out.println("\t\t    - written by Todd R. Riley");
			// System.out.print( "\n\tHi!");

			if ((!(args.length == 0 && numMandatoryArgs == 0))
                && (args.length < numMandatoryArgs || args[0].equalsIgnoreCase("-help") || args[0].equalsIgnoreCase("-?")))
            {
				System.out.println("\n Usage: " + title + " <Seq-tab-Count-File> -Option\n");
				System.out.println(" Options:");
				System.out.println("\n   -kmerLength <Int> (default = "	+ "regionLength" + ")");
				System.out.println();
				System.out.println("\n\t The length of the Kmers to count.");
				System.out.println();
				System.out.println("\n   -regionStart <Int> (default " + start + ")");
				System.out.println();
				System.out.println("\n\t The start position in the sequences to start counting kmers (0-based indexing).");
				System.out.println();
				System.out.println("\n   -regionLength <Int> (default " + "max"	+ ")");
				System.out.println();
				System.out.println("\n\t The length of the region to count kmers in.");
				System.out.println();
				System.out.println("\n   -countMin <Number> (default " + countMin + ")");
				System.out.println();
				System.out.println("\t Sets a minimum occurrence count in order to include sequence in the kmer counts.");
				System.out.println();
				System.out.println("\n   -match <regex>");
				System.out.println();
				System.out.println("\t A regex expression for a kmer within the subsequence (specified by -regionStart and -regionLength)");
				System.out.println("\t that must match in order to include the kmer in the generation of the kmer counts.");
				System.out.println();
				System.out.println("\n   -revComp <Yes/No> (default " + combineRcString + ")");
				System.out.println();
				System.out.println("\t Turns on/off combining reverse complement of the sequences.");
				System.out.println();
				System.out.println("\n   -mult <Yes/No> (default " + multString	+ ")");
				System.out.println();
				System.out.println("\t Turns on/off using the multiplicity of the sequences in the kmer counts.");
				System.out.println();
				System.out.println("\n   -skipHeader <Yes/No> (default " + skipHeaderString + ")");
				System.out.println();
				System.out.println("\t Turns on/off skipping the first line of the input file for counting.");
				System.out.println();
				System.out.println("\n   -probSkipHeader <Yes/No> (default " + probSkipHeaderString + ")");
				System.out.println();
				System.out.println("\t Turns on/off skipping the first line of the input file for probabilities.");
				System.out.println();
				System.out.println("\n   -overlapAdaptersLen <Int> (default "	+ overlapAdaptersLen + ")");
				System.out.println();
				System.out.println("\t Includes kmers that overlap at least N nucleotides with either side adapter.");
				System.out.println();
				System.out.println("\n   -sortBy <counts/kmers> (default NONE)");
				System.out.println();
				System.out.println("\t Flag to sort the kmer-tab-count list by kmers or by counts.");
				System.out.println();
				System.out.println("\n   -subKmers <Seq-tab-Count-File>");
				System.out.println();
				System.out.println("\t Loads a (k-1)mer counts file in order to calculate marginal frequencies.");
				System.out.println();
				System.out.println("\n   -ampKmers <Seq-tab-Count-File>");
				System.out.println();
				System.out.println("\t Loads a Kmer counts file in order to calculate 1 round amplification bias correction.");
				System.out.println();
				System.out.println("\n   -ampSubKmers <Seq-tab-Count-File>");
				System.out.println();
				System.out.println("\t Loads a (K-1)mer counts file in order to calculate 1 round amplification bias correction.");
				System.out.println();
				System.out.println("\n   -ampR0Kmers <Seq-tab-Count-File>");
				System.out.println();
				System.out.println("\t Loads a Kmer counts file in order to calculate 1 round amplification bias correction.");
				System.out.println();
				System.out.println("\n   -ampR0SubKmers <Seq-tab-Count-File>");
				System.out.println();
				System.out.println("\t Loads a (K-1)mer counts file in order to calculate 1 round amplification bias correction.");
				System.out.println();
				System.out.println("\n   -prob <largerKmer-tab-Count-File>");
				System.out.println();
				System.out.println("\t Parses a kmer counts file in order to calculate the probability of observing these larger kmer-counts.");
				System.out.println();
				System.out.println("\n   -probRegionStart <Int> (default "	+ probStart + ")");
				System.out.println();
				System.out.println("\n\t The start position in the sequences to start probability calc of kmers (0-based indexing).");
				System.out.println();
				System.out.println("\n   -probRegionLength <Int> (default "	+ "max" + ")");
				System.out.println();
				System.out.println("\n\t The length of the region to calc probability of kmers in.");
				System.out.println();
				System.out.println("\n   -probOverlapAdaptersLen <Int> (default " + probOverlapAdaptersLen + ")");
				System.out.println();
				System.out.println("\t include kmers that overlap at least N nucleotides with either side adapter in probability calc.");
				System.out.println();
				System.out.println("\n   -expected <format=txt/hash> <largerKmer-tab-Count-File>");
				System.out.println();
				System.out.println("\t Parses a kmer counts file in order to calculate the expected counts for these larger Kmers.");
				System.out.println();
				// System.out.println("\n   -load <table/hash/txt> <File-Name>");
				// System.out.println();
				System.out.println("\n   -save <table/hash/txt> <File-Name>");
				System.out.println();
				System.out.println("\n   -help, -? \t Displays this help message");
				System.out.println();
				System.exit(-1);
			}

			for (int a = 0; a < args.length; a++) {
				if (args[a].equalsIgnoreCase("-regionStart")) {
					start = Integer.parseInt(args[a + 1]);
				}
                else if (args[a].equalsIgnoreCase("-kmerLength")) {
					kmerLength = Integer.parseInt(args[a + 1]);
				}
                else if (args[a].equalsIgnoreCase("-regionLength")) {
					regionLength = Integer.parseInt(args[a + 1]);
				}
                else if (args[a].equalsIgnoreCase("-match")) {
					regexString = new String(args[a + 1]);
				}
                else if (args[a].equalsIgnoreCase("-revComp")) {
					combineRcString = new String(args[a + 1]);
					combineRC = StringTools.parseBoolean(combineRcString);
				}
                else if (args[a].equalsIgnoreCase("-countMin")) {
					countMin = Double.parseDouble(args[a + 1]);
				}
                else if (args[a].equalsIgnoreCase("-mult")) {
					multString = new String(args[a + 1]);
					mult = StringTools.parseBoolean(multString);
				}
                else if (args[a].equalsIgnoreCase("-skipHeader")) {
					skipHeaderString = new String(args[a + 1]);
					skipHeader = StringTools.parseBoolean(skipHeaderString);
				}
                else if (args[a].equalsIgnoreCase("-probSkipHeader")) {
					probSkipHeaderString = new String(args[a + 1]);
					probSkipHeader = StringTools.parseBoolean(probSkipHeaderString);
				}
                else if (args[a].equalsIgnoreCase("-overlapAdaptersLen")) {
					overlapAdaptersLen = Integer.parseInt(args[a + 1]);
					// overlapAdaptersString = new String(args[a+1]);
					// overlapAdapters =
					// StringTools.parseBoolean(overlapAdaptersString);
				}
                else if (args[a].equalsIgnoreCase("-sortBy")) {
					sortBy = new String(args[a + 1]);
				}
                else if (args[a].equalsIgnoreCase("-subKmers")) {
					subKmerCountsFile = new String(args[a + 1]);
				}
                else if (args[a].equalsIgnoreCase("-ampKmers")) {
					ampKmerCountsFile = new String(args[a + 1]);
				}
                else if (args[a].equalsIgnoreCase("-ampSubKmers")) {
					ampSubKmerCountsFile = new String(args[a + 1]);
				}
                else if (args[a].equalsIgnoreCase("-ampR0Kmers")) {
					ampR0KmerCountsFile = new String(args[a + 1]);
				}
                else if (args[a].equalsIgnoreCase("-ampR0SubKmers")) {
					ampR0SubKmerCountsFile = new String(args[a + 1]);
				}
                else if (args[a].equalsIgnoreCase("-prob")) {
					observedCountsFileFormat = new String(args[a + 1]);
					observedCountsFile = new String(args[a + 2]);
				}
                else if (args[a].equalsIgnoreCase("-expected")) {
					observedCountsFileFormat = new String(args[a + 1]);
					observedCountsFile = new String(args[a + 2]);
					expectedCounts = true;
				}
                else if (args[a].equalsIgnoreCase("-filter")) {
					filterPvalueThresh = Double.parseDouble(args[a + 1]);
					observedCountsFileFormat = new String(args[a + 2]);
					observedCountsFile = new String(args[a + 3]);
					filterObservedCounts = true;
				}
                else if (args[a].equalsIgnoreCase("-sortByFile")) {
					sortByFile = new String(args[a + 1]);
				}
                else if (args[a].equalsIgnoreCase("-probRegionStart")) {
					probStart = Integer.parseInt(args[a + 1]);
				}
                else if (args[a].equalsIgnoreCase("-probRegionLength")) {
					probRegionLength = Integer.parseInt(args[a + 1]);
				}
                else if (args[a].equalsIgnoreCase("-probOverlapAdaptersLen")) {
					probOverlapAdaptersLen = Integer.parseInt(args[a + 1]);
					// probOverlapAdaptersString = new String(args[a+1]);
					// probOverlapAdapters =
					// StringTools.parseBoolean(probOverlapAdaptersString);
				}
                else if (args[a].equalsIgnoreCase("-countStats")) {
					String countStatsString = new String(args[a + 1]);
					countStats = StringTools.parseBoolean(countStatsString);
				}
                else if (args[a].equalsIgnoreCase("-save")) {
					saveFileFormat = new String(args[a + 1]);
					saveFile = new String(args[a + 2]);
				}
				// else if (args[a].equalsIgnoreCase("-load")) {
				// loadFileFormat = new String(args[a+1]);
				// loadFile = new String(args[a+2]);
				// }
				// ////////////////////////////////////////////////////////////////////////////
				// ////////////////////////////////////////////////////////////////////////////
				else if (args[a].startsWith("-")) {
					System.out.println("\n\nError: Unknown argument " + args[a] + "\n");
					System.exit(-1);
				}
			}

			// System.err.println( "\nHi3!");
			proteinString = args[0];
			labelString = args[1];
			roundString = args[2];
			round = Integer.parseInt(roundString.substring(1, 2));

			String filePathName = args[3];

			// for (kmerLength <= 31)
			TLongIntHashMap kmerToIntCountMap = new TLongIntHashMap();
			TLongDoubleHashMap kmerToDoubleSumMap = new TLongDoubleHashMap();

			long kmerCountsTotal = updateCounts(kmerToIntCountMap,
                kmerToDoubleSumMap, filePathName, kmerLength, start,
                regionLength, regexString, countMin, mult,
                // overlapAdapters,
                overlapAdaptersLen, combineRC, skipHeader);
			// System.err.println("\nkmerCountsTotal = "+kmerCountsTotal);

			TLongIntHashMap subKmerToCountMap = null;
			long subKmerCountsTotal = 0;

			TLongIntHashMap ampKmerToCountMap = null;
			long ampKmerCountsTotal = 0;

			TLongIntHashMap ampSubKmerToCountMap = null;
			long ampSubKmerCountsTotal = 0;

			TLongIntHashMap ampR0KmerToCountMap = null;
			long ampR0KmerCountsTotal = 0;

			TLongIntHashMap ampR0SubKmerToCountMap = null;
			long ampR0SubKmerCountsTotal = 0;

			TLongDoubleHashMap marginalFreqs = null;

			if (ampKmerCountsFile != null) {
				ampKmerToCountMap = new TLongIntHashMap();

				ampKmerCountsTotal = updateCounts(ampKmerToCountMap, null,
                    ampKmerCountsFile, kmerLength, start, regionLength,
                    regexString, countMin, mult,
                    // overlapAdapters,
                    overlapAdaptersLen, combineRC, skipHeader);
			}

			if (ampSubKmerCountsFile != null) {
				ampSubKmerToCountMap = new TLongIntHashMap();

				ampSubKmerCountsTotal = updateCounts(ampSubKmerToCountMap,
                    null, ampSubKmerCountsFile, kmerLength - 1, start,
                    regionLength - 1, regexString, countMin, mult,
                    // overlapAdapters,
                    overlapAdaptersLen, combineRC, skipHeader);
			}

			if (ampR0KmerCountsFile != null) {
				ampR0KmerToCountMap = new TLongIntHashMap();

				ampR0KmerCountsTotal = updateCounts(ampR0KmerToCountMap, null,
                    ampR0KmerCountsFile, kmerLength, start, regionLength,
                    regexString, countMin, mult,
                    // overlapAdapters,
                    overlapAdaptersLen, combineRC, skipHeader);
			}

			if (ampR0SubKmerCountsFile != null) {
				ampR0SubKmerToCountMap = new TLongIntHashMap();

				ampR0SubKmerCountsTotal = updateCounts(ampR0SubKmerToCountMap,
                    null, ampR0SubKmerCountsFile, kmerLength - 1, start,
                    regionLength - 1, regexString, countMin, mult,
                    // overlapAdapters,
                    overlapAdaptersLen, combineRC, skipHeader);
			}

			if (subKmerCountsFile != null) {
				subKmerToCountMap = new TLongIntHashMap();

				subKmerCountsTotal = updateCounts(subKmerToCountMap, null,
                    subKmerCountsFile, kmerLength - 1, start,
                    regionLength - 1, regexString, countMin, mult,
                    // overlapAdapters,
                    overlapAdaptersLen, combineRC, skipHeader);

				// System.err.println("\nsubKmerCountsTotal = "+subKmerCountsTotal);

			}

			TLongIntHashMap observedCountMap = null;

			// expectedCountMap holds either probs or expected counts
			TLongDoubleHashMap expectedCountMap = null;
			long observedCountsTotal = 0;

			Table outputTable = null;

			// ///////////////////////////////////////////////////////////////////
			// calculate expected counts
			// ///////////////////////////////////////////////////////////////////
			if (observedCountsFile != null) {

				if (observedCountsFileFormat.equalsIgnoreCase("txt")) {
					// observedCountMap =
					// readTLongIntHashMap(observedCountsFile, "\t", 1, 0,
					// Integer.MAX_VALUE);
				}
                else if (observedCountsFileFormat.equalsIgnoreCase("hash")) {
					observedCountMap = (TLongIntHashMap) FileTools.readSerializedFile(observedCountsFile);

					if (!observedCountMap.containsKey(-1L)) {
						observedCountMap.put(-1L, 12);
					}

				}

				// String observedCountsBaseName =
				// FileTools.stripPath(observedCountsFile);
				// String[] pars = observedCountsBaseName.split("\\.");

				// for (int i=0; i < pars.length; i++) {
				// if ( pars[i].matches("[rR]\d+")) {
				// roundString = pars[i];
				// round = Integer.parseInt(roundString.substring(1,2));
				// break;
				// }
				// }

				// proteinString = pars[0];
				// labelString = pars[1];

				// // roundString = pars[2];
				// // round = Integer.parseInt(pars[2].substring(1));

				// if (pars[2].startsWith("r") || pars[2].startsWith("R")) {
				// round = Integer.parseInt(pars[2].substring(1,2));
				// roundString = pars[2];
				// }
				// else {
				// round = Integer.parseInt(pars[2].substring(0,1));
				// // //roundString = "R"+round;
				// roundString = "R"+pars[2];
				// // roundString = "R0";
				// }

				if (subKmerCountsFile != null) {

					marginalFreqs = marginalFreqs(kmerToIntCountMap,
                        kmerCountsTotal, subKmerToCountMap,
                        subKmerCountsTotal, ampKmerToCountMap,
                        ampKmerCountsTotal, ampSubKmerToCountMap,
                        ampSubKmerCountsTotal, ampR0KmerToCountMap,
                        ampR0KmerCountsTotal, ampR0SubKmerToCountMap,
                        ampR0SubKmerCountsTotal, kmerLength, round);
				}

				// Use the intermediate HashMap if we need to sort or if we are
				// filtering, or we need to get expected counts
				// if ((sortBy != null) || expectedCounts ||
				// filterObservedCounts) {
				if ((sortBy != null) || filterObservedCounts) {
					// populate the expectedCountMap (hashmap)
					expectedCountMap = new TLongDoubleHashMap();
				}
                else { // sortBy == null
                    // populate the outputTable
					outputTable = new Table();
					// System.err.println( "Hi1!");
				}

				MutableInt keyLength = new MutableInt(-1);

				if (observedCountMap == null) {
					// Parse the observedCountsFile one line at a time
					observedCountsTotal = probabilities(expectedCountMap,
                        outputTable, observedCountsFile, kmerToIntCountMap,
                        kmerCountsTotal, ampKmerToCountMap,
                        ampKmerCountsTotal, ampR0KmerToCountMap,
                        ampR0KmerCountsTotal, round, marginalFreqs,
                        kmerLength, keyLength, probStart, probRegionLength,
                        // probOverlapAdapters);
                        probOverlapAdaptersLen, probSkipHeader);
				}
                else {
					// Use the cached observedCountsMap
					observedCountsTotal = probabilities(expectedCountMap,
                        outputTable, observedCountsFile, observedCountMap,
                        kmerToIntCountMap, kmerCountsTotal,
                        ampKmerToCountMap, ampKmerCountsTotal,
                        ampR0KmerToCountMap, ampR0KmerCountsTotal, round,
                        marginalFreqs, kmerLength, keyLength, probStart,
                        probRegionLength,
                        // probOverlapAdapters);
                        probOverlapAdaptersLen);
				}

				// if we want expected counts then multiply all the
				// probabilities in expectedCountMap
				// by the observedCountsTotal
				// System.err.println("\nobservedCountsTotal = "+observedCountsTotal);

				if (expectedCounts) {
					// if (expectedCounts || filterObservedCounts) {
					getExpectedCounts(expectedCountMap, outputTable, observedCountsTotal);
				}

				// If output type is hash then don't create the sorted table or
				// dumpString
				// If output type is sorted table then don't create the
				// dumpString
				if ((!(saveFileFormat != null && saveFileFormat.equalsIgnoreCase("hash"))) && (outputTable == null)) {
					outputTable = getSortedKmerCounts(expectedCountMap, sortBy);
					// System.err.println( "Hi2!");
				}

				// new merLength is the length of the substring(observed seqs)
				kmerLength = keyLength.intValue();

			}

			// ///////////////////////////////////////////////////////////////////
			// DON'T calculate probabilities or expected counts
			// If output type is hash then don't create the sorted table or
			// dumpString
			// If output type is sorted table then don't create the dumpString
			// ///////////////////////////////////////////////////////////////////
			else {

				// System.err.println( "Hi3!");
				// System.err.println(
				// "\nkmerToIntCountMap.size() == "+kmerToIntCountMap.size());
				// System.err.println(
				// "\nkmerToDoubleSumMap.size() == "+kmerToDoubleSumMap.size());

				if (!(saveFileFormat != null && saveFileFormat.equalsIgnoreCase("hash"))) {
					// System.err.println( "Hi4!");

					// Perform proper sorting
					if (sortByFile == null) {
						// System.err.println( "Hi5!");
						if (kmerToIntCountMap.size() > 0) {
							outputTable = getSortedKmerCounts(kmerToIntCountMap, sortBy);
						}
                        else {
							// System.err.println( "Hi6!");
							outputTable = getSortedKmerCounts(
                                kmerToDoubleSumMap, sortBy);
						}
					}
                    else {
						outputTable = getSortedKmerCountsForFile(
                            kmerToIntCountMap, sortByFile);
					}

				}

				// String countsTableBaseName =
				// FileTools.stripPath(filePathName);
				// String[] pars = countsTableBaseName.split("\\.");
				// proteinString = pars[0];
				// labelString = pars[1];

				// // //round = Integer.parseInt(pars[2].substring(0,1));
				// // //roundString = "R"+round;
				// // roundString = "R"+pars[2];

				// roundString = pars[2];
			}

			// ///////////////////////////////////////////////////////////////////
			// perform
			// ///////////////////////////////////////////////////////////////////
			if (filterObservedCounts) {

				outputTable = filterKmers(filterPvalueThresh, observedCountMap,
                    expectedCountMap, observedCountsTotal);

			}
			// ///////////////////////////////////////////////////////////////////
			// save output
			// ///////////////////////////////////////////////////////////////////
			else if (saveFileFormat == null) { // no output file so dump to
                // stdout
				String dumpString = null;

				if (countStats) {
					dumpString = getCountStats(kmerToIntCountMap, kmerLength);
				}
                else {
					// System.err.println( "Hi4!");
					dumpString = toString(outputTable, kmerLength,
                        proteinString, labelString, roundString,
                        leftAdapterSeq, rightAdapterSeq);
					// String dumpString = toString(outputTable, kmerLength,
					// protein, roundString);
				}

				System.out.println(dumpString);

			}
            else if (saveFileFormat.equalsIgnoreCase("txt")) { // format ==
                // txt
				// save the table directly to file!!
				saveTableToText(saveFile, false, outputTable, kmerLength,
                    proteinString, labelString, roundString,
                    leftAdapterSeq, rightAdapterSeq);
			}
            else if (saveFileFormat.equalsIgnoreCase("table")) { // format ==
                // serialized
                // table
				FileTools.writeSerializedFile(outputTable, saveFile);
			}
            else if (saveFileFormat.equalsIgnoreCase("hash")) { // format ==
                // serialized
                // hash
				if (expectedCountMap != null) {
					FileTools.writeSerializedFile(expectedCountMap, saveFile);
				}
                else {
					FileTools.writeSerializedFile(kmerToIntCountMap, saveFile);
				}
			}
            else { // unknown save format
				System.out.println("Error: unknown save format " + saveFileFormat + " not recognized.");
			}

		}

		catch (Exception ex) {
			ex.printStackTrace();
		}

	}

	// assumes key=column0 and value=column1
	// 0-based indexing
	public static TLongIntHashMap readTLongIntHashMap(String aFilePathName,
        String delimiter, int numLinesSkip) {
		return (readTLongIntHashMap(aFilePathName, delimiter, numLinesSkip, 0,
				1));
	}

	public static TLongIntHashMap readTLongIntHashMap(String aFilePathName,
        String delimiter, int numLinesSkip, int keyColumn, int valueColumn) {

		Pattern delimiterPattern = Pattern.compile(delimiter);

		String keyValuePairs[] = FileTools.readStrings(aFilePathName,
            numLinesSkip);
		TLongIntHashMap tLongIntHM = new TLongIntHashMap(keyValuePairs.length);
		boolean columnsNumsChecked = false;
		boolean longIsLongStringChecked = false;
		boolean longIsLongString = false;

		for (String keyValuePair : keyValuePairs) {

			String[] keyValuePairArray = delimiterPattern.split(keyValuePair);

			if (keyValuePairArray.length >= 2) {

				if (!columnsNumsChecked) {
					if (keyColumn == Integer.MAX_VALUE) {
						keyColumn = keyValuePairArray.length - 1;
					}
					if (valueColumn == Integer.MAX_VALUE) {
						valueColumn = keyValuePairArray.length - 1;
					}
					columnsNumsChecked = true;
				}

				if (!longIsLongStringChecked) {
					try {
						Long.parseLong((String) keyValuePairArray[keyColumn]);
						longIsLongString = true;
					} catch (Exception e) {
						longIsLongString = false;

						// Set the Word Length!!!!!
						tLongIntHM.put(
                            -1L,
                            ((String) keyValuePairArray[keyColumn])
                            .length());

					}
					longIsLongStringChecked = true;
				}

				long key;
				if (longIsLongString) {
					key = (long) Long
                        .parseLong((String) keyValuePairArray[keyColumn]);
				}
                else {
					key = (long) getWordLong((String) keyValuePairArray[keyColumn]);
				}

				int retrievedValue = tLongIntHM.get(key);

				if (valueColumn < keyValuePairArray.length) {
					tLongIntHM
                        .put(key,
                            retrievedValue
                            + Integer
                            .parseInt(keyValuePairArray[valueColumn]));
				}
			}
		}
		return (tLongIntHM);
	}

	public TLongDoubleHashMap marginalFreqs(TLongIntHashMap kmerToIntCountMap,
        long kmerCountsTotal, TLongIntHashMap subKmerToCountMap,
        long subKmerCountsTotal, TLongIntHashMap ampKmerToCountMap,
        long ampKmerCountsTotal, TLongIntHashMap ampSubKmerToCountMap,
        long ampSubKmerCountsTotal, TLongIntHashMap ampR0KmerToCountMap,
        long ampR0KmerCountsTotal, TLongIntHashMap ampR0SubKmerToCountMap,
        long ampR0SubKmerCountsTotal, int kmerLength, int round) {

		TLongDoubleHashMap marginalFreqs = new TLongDoubleHashMap();
		for (TLongIntIterator iter = kmerToIntCountMap.iterator(); iter
                 .hasNext();) {
			iter.advance();

			// SKIP THE KEY THAT GIVES THE kmerLength
			if (iter.key() == -1L) {
				continue;
			}

			long kmerKey = iter.key();
			char[] kmer = getWordCharArray(kmerKey, kmerLength);
			char[] subKmer = Arrays.copyOf(kmer, kmerLength - 1);
			long subKmerKey = getWordLong(subKmer);

			int kmerCount = iter.value();
			int subKmerCount = subKmerToCountMap.get(subKmerKey);
			// System.err.print( "\nkmerCount = "+kmerCount);
			// System.err.print( "\nkmerCountsTotal = "+kmerCountsTotal);
			// System.err.print( "\nsubKmerCount = "+subKmerCount);
			// System.err.print( "\nsubKmerCountsTotal = "+subKmerCountsTotal);

			double probKmer = ((double) kmerCount) / kmerCountsTotal;
			double probSubKmer = ((double) subKmerCount) / subKmerCountsTotal;
			double marginalFreq = (probKmer) / (probSubKmer);

			if (ampKmerToCountMap != null) {
				int ampKmerCount = ampKmerToCountMap.get(kmerKey);
				int ampSubKmerCount = ampSubKmerToCountMap.get(subKmerKey);

				double probAmpKmer = ((double) ampKmerCount)
                    / ampKmerCountsTotal;
				double probAmpSubKmer = ((double) ampSubKmerCount)
                    / ampSubKmerCountsTotal;

				int ampR0KmerCount = ampR0KmerToCountMap.get(kmerKey);
				int ampR0SubKmerCount = ampR0SubKmerToCountMap.get(subKmerKey);

				double probAmpR0Kmer = ((double) ampR0KmerCount)
                    / ampR0KmerCountsTotal;
				double probAmpR0SubKmer = ((double) ampR0SubKmerCount)
                    / ampR0SubKmerCountsTotal;

				// double ampCorrection = ( (probAmpKmer/probAmpR0Kmer) /
				// (probAmpSubKmer/probAmpR0SubKmer) )^round;
				double ampCorrection = Math.pow((probAmpKmer / probAmpR0Kmer)
                    / (probAmpSubKmer / probAmpR0SubKmer), round);

				marginalFreq *= ampCorrection;
			}

			marginalFreqs.put(kmerKey, marginalFreq);
		}
		return (marginalFreqs);
	}

	// public double probability(
	// String observedCountsFile,
	// TLongIntHashMap kmerToIntCountMap,
	// int kmerCountsTotal,
	// TLongDoubleHashMap marginalFreqs,
	// int kmerLength) {

	// double probSum = 0;
	// int numSeqs = 0;
	// try {
	// BufferedReader bufferedReader = new BufferedReader(new
	// FileReader(observedCountsFile));

	// // Parse the table and calculate the average probability
	// String countsTableBaseName = FileTools.stripPath(observedCountsFile);
	// System.err.print(
	// "Parsing "+countsTableBaseName+" to determine the probabilities of generating the sequences....");
	// String lineString;
	// boolean firstLine = true;
	// while ( (lineString = bufferedReader.readLine()) != null ) {
	// //String lineEntries[] = lineString.split("\t");
	// String lineEntries[] = tabsPattern.split(lineString);
	// if (lineEntries.length > 1) {
	// String seq = lineEntries[0];
	// //String seqUpperCase = lineEntries[0].toUpperCase();
	// int seqReadCount = 0;
	// if (lineEntries.length == 2) {
	// seqReadCount = Integer.parseInt(lineEntries[1]);
	// }
	// else {
	// if (firstLine) { // skip the R colnames header
	// firstLine = false;
	// continue;
	// }
	// else {
	// seqReadCount = Integer.parseInt(lineEntries[4]);
	// }
	// }

	// // get probability of seeing this seq given the markov model
	// // markov model consists of kmerToIntCountMap, kmerCpuntsTotal,
	// marginalFreqs
	// double prob = probability(seq, kmerToIntCountMap, kmerCountsTotal,
	// marginalFreqs);

	// // Now need to perform binomial
	// prob *= seqReadCount; // not correct

	// numSeqs += seqReadCount;
	// probSum += prob;
	// }
	// }
	// System.err.println( "Done.");

	// }
	// catch (Exception ex) {
	// ex.printStackTrace();
	// }
	// //System.err.print( "\nprobSum = "+probSum);
	// //System.err.print( "\nnumSeqs = "+numSeqs);
	// return(probSum/numSeqs);

	// }

	public Table filterKmers(double filterPvalueThresh,
        TLongIntHashMap observedCountMap,
        TLongDoubleHashMap probabilityMap, long observedCountsTotal) {
		Table contaminationTable = new Table();

		try {
			int totalLinesCount = observedCountMap.size() - 1; // one element is
            // the
            // kmerLength
			double singleLinePercent = (1.0 / (totalLinesCount - 0)) * 100; // -1
																			// is
																			// for
																			// the
																			// header
																			// row
			int percentStepSize = Math.max(10,
                (int) Math.floor(singleLinePercent));
			int numberStepSize = (int) Math.floor(totalLinesCount
                / percentStepSize);
			int nextPercentThreshold = percentStepSize;
			int nextNumberThreshold = numberStepSize;

			System.err.print("Parsing all the Kmers to determine if any have revcomp count assymmetry greater than the contamination threshold....");

			int kmerLength = observedCountMap.get(-1L);

			int lineCounter = 0;
			for (TLongIntIterator iterator = observedCountMap.iterator(); iterator
                     .hasNext();) {

				iterator.advance();

				// Update parsing percentage
				lineCounter++;
				if ((lineCounter >= nextNumberThreshold)
                    && (nextPercentThreshold <= 100)) {
					System.err.print("\n\t" + nextPercentThreshold + "% complete...");
					nextPercentThreshold += percentStepSize;
					nextNumberThreshold += numberStepSize;
				}

				// SKIP THE KEY THAT GIVES THE kmerLength
				if (iterator.key() == -1L) {
					continue;
				}

				long kmerWord = iterator.key();
				int observedCount = iterator.value();
				double expectedCount = probabilityMap.get(kmerWord)
                    * observedCountsTotal;

				long revCompKmerWord = getRevComp(kmerWord, kmerLength);
				int observedRevCompCount = observedCountMap
                    .get(revCompKmerWord);
				double expectedRevCompCount = probabilityMap
                    .get(revCompKmerWord) * observedCountsTotal;

				double observedRatio = ((double) observedCount)
                    / (observedRevCompCount + observedCount);
				double expectedRatio = expectedCount
                    / (expectedRevCompCount + expectedCount);

				double ratiosRatio = (((double) observedCount) / observedRevCompCount)
                    / (expectedCount / expectedRevCompCount);

				if ((observedCount > 4)
                    && (observedRatio > (expectedRatio * 1.1))) {
					// if ((observedCount > 4) && (observedRatio >
					// (expectedRatio*1.1)) && ((observedRevCompCount*1.5) <
					// observedCount)) {
					// if ((observedCount > 4) && (observedRatio >
					// expectedRatio) && (observedRevCompCount < observedCount))
					// {
					// if ((observedCount > 4) && ((observedRevCompCount*2) <
					// observedCount)) {

					// /////////////////////////////////////////////////////////////////////////////////////
					// gross ratios of Obs / (Obs + Obs_rc)
					// /////////////////////////////////////////////////////////////////////////////////////
					double cumulativeProb1;
					// // Binomial Distribution
					BinomialDistribution binomialDist = new BinomialDistribution(
                        observedCount + observedRevCompCount, expectedRatio);
					cumulativeProb1 = binomialDist.cumulativeProbability(observedCount);

					// // Poisson Distribution - approximates the binomial for
					// large n, small p
					// PoissonDistribution poissonDist = new
					// PoissonDistribution(expectedRatio, .001);
					// double cumulativeProb =
					// poissonDist.cumulativeProbability(observedRatio);
					// double pValue = 1 - cumulativeProb;

					// Normal Distribution - approximates the poisson for lambda
					// > 10
					// Continuity correction for lambda < 1000 ???????
					// NormalDistribution normalDist = new
					// NormalDistribution(expectedRatio,
					// Math.sqrt(expectedRatio));
					// double cumulativeProb =
					// normalDist.cumulativeProbability(observedRatio);

					// /////////////////////////////////////////////////////////////////////////////////////
					// gross observations of word
					// /////////////////////////////////////////////////////////////////////////////////////
					double cumulativeProb2;
					// if (true) {
					// }
					// else if (true) {
					if (true) {
						// Poisson Distribution - the poisson is good
						// approximation of the binomial when (n >= 20 and (p <=
						// 0.05 or p>= 0.95)
						PoissonDistribution poissonDist = new PoissonDistribution(
                            expectedCount);
						// PoissonDistribution poissonDist = new
						// PoissonDistribution(expectedCount, .0000000001);
						cumulativeProb2 = poissonDist.cumulativeProbability(observedCount);
					}
                    else {
						// if n*p >= 100 (large enough) then we can well
						// approximate the poisson with the normal distribution

						// to approximate the binomial with the normal then also
						// p should not be close to 0 or 1 (closer to .5 the
						// better)
						// (n*p >= 10) and (n*(1-p) >= 10)
						NormalDistribution normalDist = new NormalDistribution(
                            expectedCount, Math.sqrt(expectedCount));
						// Continuity correction is the "half-correction"!!
						cumulativeProb2 = normalDist.cumulativeProbability(observedCount + 0.5);
					}

					double pValue1 = 1 - cumulativeProb1;
					if (pValue1 == 0) {
						pValue1 = 1e-16;
					}
					double pValue2 = 1 - cumulativeProb2;
					if (pValue2 == 0) {
						pValue2 = 1e-16;
					}

					// Bonferroni correction
					pValue1 *= (observedCountMap.size() - 1); // one element is
                    // the
                    // kmerLength
					pValue2 *= (observedCountMap.size() - 1); // one element is
                    // the
                    // kmerLength

					// if (true) {
					if ((pValue1 < filterPvalueThresh)
                        || (pValue2 < filterPvalueThresh)) {

						// System.out.print(
						// "\n"+getWordString(kmerWord, kmerLength)+"\t"+
						// expectedCount+"\t"+
						// observedCount+"\t"+
						// expectedRatio+"\t"+
						// observedRatio+"\t"+
						// pValue+"\t"
						// );

						contaminationTable.add(
                            Arrays.asList(
                                getWordString(kmerWord, kmerLength),
								new Double(expectedCount),
                                new Integer(observedCount),
                                new Double(expectedRevCompCount),
                                new Integer(observedRevCompCount),
                                new Double(expectedRatio),
                                new Double(observedRatio),
								new Double(ratiosRatio), new Double(pValue2),
								new Double(pValue1)));

					}
				}

				// if (lineCounter > 100) {
				// break;
				// }

			}
			System.err.println("\nDone.");

			// int[] sortColumns = {9, 7};
			int[] sortColumns = { 8, 7 };
			int[] sortDirections = { 1, -1 };
			contaminationTable.sort(sortColumns, sortDirections); // sort by the
            // p-values
			// contaminationTable.reverse();

			System.out.print("\nSequence\tExpected Count\tObserved Count\tExpected RevComp Count\tObserved RevComp Count\tExp / (Exp + Exp_rc)\tObs / (Obs + Obs_rc)\t(Obs/Obs_rc) / (Exp/Exp_rc)\tp-value1\tp-value2");
			System.out.println("\n" + contaminationTable.toString());

		} catch (Exception ex) {
			ex.printStackTrace();
		}
		return (contaminationTable);
	}

	// /////////////////////////////////////////////////////////
	// Parse the observedCountsFile one line at a time
	// /////////////////////////////////////////////////////////
	public long probabilities(TLongDoubleHashMap probsHashMap,
        Table probsTable, String observedCountsFile,
        TLongIntHashMap kmerToIntCountMap, long kmerCountsTotal,
        TLongIntHashMap ampKmerToCountMap, long ampKmerCountsTotal,
        TLongIntHashMap ampR0KmerToCountMap, long ampR0KmerCountsTotal,
        int round, TLongDoubleHashMap marginalFreqs, int kmerLength,
        MutableInt keyLength, int probStart, int probRegionLength,
        // boolean probOverlapAdapters
        int probOverlapAdaptersLen, boolean skipHeader) {

		// TLongDoubleHashMap probsHashMap = new TLongDoubleHashMap();
		long countsTotal = 0;
		try {
			int totalElementsCount = FileTools.getLineCount(observedCountsFile);
			double singleLinePercent = (1.0 / (totalElementsCount - 0)) * 100; // -1
            // is
            // for
            // the
            // header
            // row
			int percentStepSize = Math.max(10, (int) Math.floor(singleLinePercent));
			int numberStepSize = (int) Math.floor(totalElementsCount / percentStepSize);
			int nextPercentThreshold = percentStepSize;
			int nextNumberThreshold = numberStepSize;

			BufferedReader bufferedReader = null;
			if (observedCountsFile.endsWith(".gz")) {
				bufferedReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(observedCountsFile))));
			}
            else {
				bufferedReader = new BufferedReader(new FileReader(observedCountsFile));
			}

			// Parse the table and calculate the average probability
			String countsTableBaseName = FileTools.stripPath(observedCountsFile);
			System.err.print("Parsing " + countsTableBaseName + " to determine the probabilities of generating the sequences....");
			String lineString;
			boolean firstLine = true;
			boolean keyLengthSet = false;

			int lineCounter = 0;
			while ((lineString = bufferedReader.readLine()) != null) {

				// Update parsing percentage
				lineCounter++;
				if ((lineCounter >= nextNumberThreshold) && (nextPercentThreshold <= 100)) {
					System.err.print("\n\t" + nextPercentThreshold + "% complete...");
					// System.err.print( "\n\tHi!");
					nextPercentThreshold += percentStepSize;
					nextNumberThreshold += numberStepSize;
				}

				// skip the R colnames header
				if (skipHeader && firstLine) {
					firstLine = false;
					continue;
				}

				// String lineEntries[] = lineString.split("\t");
				String lineEntries[] = tabsPattern.split(lineString);
				if (lineEntries.length > 1) {
					// String seq = lineEntries[0];
					String seqUpperCase = lineEntries[0].toUpperCase();
					int seqReadCount = 0;
					if (lineEntries.length == 2) {
						seqReadCount = Integer.parseInt(lineEntries[1]);
					}
                    else {
						seqReadCount = Integer.parseInt(lineEntries[4]);
					}

					// get probability of seeing this seq given the markov model
					// markov model consists of kmerToIntCountMap,
					// kmerCountsTotal, marginalFreqs
					double prob = probability(
                        seqUpperCase,
                        kmerToIntCountMap,
                        kmerCountsTotal,
                        ampKmerToCountMap,
                        ampKmerCountsTotal,
                        ampR0KmerToCountMap,
                        ampR0KmerCountsTotal,
                        round,
                        marginalFreqs,
                        kmerLength,
                        probStart,
                        probRegionLength,
                        // probOverlapAdapters);
                        probOverlapAdaptersLen);

					long wordLong = -1;

					if (probRegionLength == -1) {
						wordLong = getWordLong(seqUpperCase);
					}
                    else {
						String substring = seqUpperCase.substring(probStart, (probStart + probRegionLength));
						wordLong = getWordLong(substring);
						// System.err.println("wordString="+substring);
						// System.err.println("wordLong="+wordLong);
					}

					if (!keyLengthSet) {
						if (probRegionLength == -1) {
							keyLength.setValue(seqUpperCase.length());
						}
                        else {
							keyLength.setValue(probRegionLength);
							// leftAdapterSeq = seqUpperCase.substring(0,
							// probStart);
							// rightAdapterSeq =
							// seqUpperCase.substring(probStart+probRegionLength,
							// seqUpperCase.length());
						}
						if (probsHashMap != null) {
							probsHashMap.put(-1L, keyLength.intValue());
						}

						keyLengthSet = true;
						// System.err.println("keyLength="+keyLength.intValue());

					}

					// Store in either the hashmap or the table
					if (probsHashMap != null) {
						probsHashMap.put(wordLong, prob);
					}
                    else {
						probsTable.add(Arrays.asList(new Long(wordLong),
								new Double(prob)));
					}

					countsTotal += seqReadCount;
				}

				// System.exit(-1);
			}
			System.err.println("Done.");

		} catch (Exception ex) {
			ex.printStackTrace();
		}
		return (countsTotal);
		// return(probsHashMap);
	}

	// /////////////////////////////////////////////////////////
	// Use the cached observedCountsMap
	// /////////////////////////////////////////////////////////
	public long probabilities(TLongDoubleHashMap probsHashMap,
        Table probsTable, String observedCountFile,
        TLongIntHashMap observedCountMap,
        TLongIntHashMap kmerToIntCountMap, long kmerCountsTotal,
        TLongIntHashMap ampKmerToCountMap, long ampKmerCountsTotal,
        TLongIntHashMap ampR0KmerToCountMap, long ampR0KmerCountsTotal,
        int round, TLongDoubleHashMap marginalFreqs, int kmerLength,
        MutableInt keyLength, int probStart, int probRegionLength,
        // boolean probOverlapAdapters
        int probOverlapAdaptersLen) {

		// TLongDoubleHashMap probsHashMap = new TLongDoubleHashMap();
		long countsTotal = 0;
		try {
			int totalElementsCount = observedCountMap.size();
			double singleLinePercent = (1.0 / (totalElementsCount - 0)) * 100; // -1
            // is
            // for
            // the
            // header
            // row
			int percentStepSize = Math.max(10,
                (int) Math.floor(singleLinePercent));
			int numberStepSize = (int) Math.floor(totalElementsCount
                / percentStepSize);
			int nextPercentThreshold = percentStepSize;
			int nextNumberThreshold = numberStepSize;

			// Parse the table and calculate the average probability
			String countsTableBaseName = FileTools.stripPath(observedCountFile);
			System.err.print("Parsing " + countsTableBaseName + " to determine the probabilities of generating the sequences....");

			// The kmer length is stored in (key = -1)
			int fullObservedSeqLength = observedCountMap.get(-1L);
			if (fullObservedSeqLength == 0) {
				fullObservedSeqLength = 12;
			}

			// use the stored kmer length of the full sequence as the keyLength
			// or use the length of the substring
			if (probRegionLength == -1) {
				// key length is the length of the sequences that are used to
				// calc probs
				keyLength.setValue(fullObservedSeqLength);
			}
            else {
				keyLength.setValue(probRegionLength);
			}
			if (probsHashMap != null) {
				probsHashMap.put(-1L, keyLength.intValue());
			}

			int lineCounter = 0;
			for (TLongIntIterator iterator = observedCountMap.iterator(); iterator
                     .hasNext();) {

				iterator.advance();

				// Update parsing percentage
				lineCounter++;
				if ((lineCounter >= nextNumberThreshold) && (nextPercentThreshold <= 100)) {
					System.err.print("\n\t" + nextPercentThreshold + "% complete...");
					// System.err.print( "\n\tHi!");
					nextPercentThreshold += percentStepSize;
					nextNumberThreshold += numberStepSize;
				}

				// SKIP THE KEY THAT GIVES THE kmerLength
				if (iterator.key() == -1L) {
					continue;
				}

				// System.err.println(
				// "\nfullObservedSeqLength = "+fullObservedSeqLength);
				String seqUpperCase = getWordString(iterator.key(), fullObservedSeqLength).toUpperCase();
				int seqReadCount = iterator.value();

				// get probability of seeing this seq given the markov model
				// markov model consists of kmerToIntCountMap, kmerCountsTotal,
				// marginalFreqs
				double prob = probability(
                    seqUpperCase,
                    kmerToIntCountMap,
                    kmerCountsTotal,
                    ampKmerToCountMap,
                    ampKmerCountsTotal,
                    ampR0KmerToCountMap,
                    ampR0KmerCountsTotal,
                    round,
                    marginalFreqs,
                    kmerLength,
                    probStart,
                    probRegionLength,
                    // probOverlapAdapters);
                    probOverlapAdaptersLen);

				// //////////////////////////////////////////////////////////////////
				// get the wordLong for this full seq or the desired substring
				// //////////////////////////////////////////////////////////////////
				long wordLong;
				if (probRegionLength == -1) {
					wordLong = iterator.key();
				}
                else {
					String substring = seqUpperCase.substring(probStart, (probStart + probRegionLength));
					wordLong = getWordLong(substring);
					// System.err.println("wordString="+substring);
					// System.err.println("wordLong="+wordLong);
				}

				// //////////////////////////////////////////////////////////////////
				// Store in either the hashmap or the table
				// //////////////////////////////////////////////////////////////////
				if (probsHashMap != null) {
					probsHashMap.put(wordLong, prob);
				}
                else {
					probsTable.add(Arrays.asList(new Long(wordLong), new Double(prob)));
				}

				countsTotal += seqReadCount;

				// System.exit(-1);
			}
			System.err.println("\nDone.");

		} catch (Exception ex) {
			ex.printStackTrace();
		}
		return (countsTotal);
		// return(probsHashMap);
	}

	// get probability of seeing this seq given the markov model
	// markov model consists of kmerToIntCountMap, kmerCountsTotal,
	// marginalFreqs
	public double probability(String seq, TLongIntHashMap kmerToIntCountMap,
        long kmerCountsTotal, TLongIntHashMap ampKmerToCountMap,
        long ampKmerCountsTotal, TLongIntHashMap ampR0KmerToCountMap,
        long ampR0KmerCountsTotal, int round,
        TLongDoubleHashMap marginalFreqs, int kmerLength, int probStart,
        int probRegionLength,
        // boolean probOverlapAdapters)
        int probOverlapAdaptersLen) {
		String seqUpperCase = seq.toUpperCase();
		boolean getMarginal = false;
		double prob = 1;

		int firstStartPosInc;
		int lastStartPosExc;

		// if (probOverlapAdapters) {
		// // We want to include all kmer windows that overlap with the Variable
		// Region (by atleast 1 nt)
		// firstStartPosInc = (probStart+1)-kmerLength;
		// lastStartPosExc = probStart+probRegionLength;
		// }
		// else {
		// // We only want kmer windows that are completely contained in the
		// variable region (no overlapping with adapters)
		// firstStartPosInc = probStart;

		// if (probRegionLength == -1) {
		// lastStartPosExc = seqUpperCase.length()-kmerLength+1;
		// }
		// else {
		// lastStartPosExc = probStart+probRegionLength-kmerLength+1;
		// }
		// }
		firstStartPosInc = probStart - probOverlapAdaptersLen;

		if (probRegionLength == -1) {
			// compute over the whole sequence
			lastStartPosExc = seqUpperCase.length() - kmerLength + 1
                + probOverlapAdaptersLen;
		}
        else {
			// compute over a substring
			lastStartPosExc = probStart + probRegionLength - kmerLength + 1
                + probOverlapAdaptersLen;
		}

		// System.out.println( "\nkmerCountsTotal = "+kmerCountsTotal);
		String fullSeq = seqUpperCase.substring(firstStartPosInc, probStart
            + probRegionLength + probOverlapAdaptersLen);
		// System.out.println("\n"+fullSeq);

		// for (int startPos = 0; startPos < seqUpperCase.length()-kmerLength+1;
		// startPos++) {
		for (int startPos = firstStartPosInc; startPos < lastStartPosExc; startPos++) {

			String currentKmer = seqUpperCase.substring(startPos,
                (startPos + kmerLength));
			// System.out.println( currentKmer );
			long wordLong = getWordLong(currentKmer);

			if (getMarginal == false) {
				// System.out.print(
				// "\nkmerToIntCountMap.get(wordLong) = "+kmerToIntCountMap.get(wordLong));
				int kmerCount = kmerToIntCountMap.get(wordLong);
				double probKmer = ((double) kmerCount) / kmerCountsTotal;

				prob *= probKmer;

				if (ampKmerToCountMap != null) {
					int ampKmerCount = ampKmerToCountMap.get(wordLong);
					double probAmpKmer = ((double) ampKmerCount) / ampKmerCountsTotal;

					int ampR0KmerCount = ampR0KmerToCountMap.get(wordLong);
					double probAmpR0Kmer = ((double) ampR0KmerCount) / ampR0KmerCountsTotal;

					// double ampCorrection = (probAmpKmer/probAmpR0Kmer)^round;
					double ampCorrection = Math.pow(probAmpKmer / probAmpR0Kmer, round);

					prob *= ampCorrection;
				}

				if (kmerLength > 1) {
					getMarginal = true;
				}
			}
            else { // getMarginal == true
				prob *= marginalFreqs.get(wordLong);
			}
		}
		// System.out.println("\nprob="+prob+"\n");
		// System.out.flush();
		// System.exit(0);
		return (prob);
	}

	public long updateCounts(TLongIntHashMap kmerToIntCountMap,
        TLongDoubleHashMap kmerToDoubleSumMap, String filePathName,
        int kmerLength, int start, int regionLength, String regexString,
        double countMin, boolean mult,
        // boolean overlapAdapters,
        int overlapAdaptersLen, boolean combineRC, boolean skipHeader) {

		// Java HashMap class
		// Map<String, Integer> kmerToIntCountMap = new HashMap();

		// Trove TIntIntHashMap is faster and has higher capacity
		// for (kmerLength <= 15)
		// TIntIntHashMap kmerToIntCountMap = new TIntIntHashMap();

		boolean kmerLengthSet = false;
		boolean useKmerToIntCountMap = true;

		long countsTotal = 0;
		try {

			if (kmerLength == 0) {
                // System.out.println("Hi");
				kmerLength = regionLength;
			}

			int totalLinesCount = FileTools.getLineCount(filePathName);
			double singleLinePercent = (1.0 / (totalLinesCount - 0)) * 100; // -1 is for the header row
			int percentStepSize = Math.max(10, (int) Math.floor(singleLinePercent));
			int numberStepSize = (int) Math.floor(totalLinesCount / percentStepSize);
			int nextPercentThreshold = percentStepSize;
			int nextNumberThreshold = numberStepSize;

			BufferedReader bufferedReader = null;
			if (filePathName.endsWith(".gz")) {
				bufferedReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filePathName))));
			}
            else {
				bufferedReader = new BufferedReader(new FileReader(filePathName));
			}

			int firstStartPosInc;
			int lastStartPosExc;

			// if (overlapAdapters) {
			// // We want to include all kmer windows that overlap with the
			// Variable Region (by atleast 1 nt)
			// firstStartPosInc = (start+1)-kmerLength;
			// lastStartPosExc = start+regionLength;
			// }
			// else {
			// // We only want kmer windows that are completely contained in the
			// variable region (no overlapping with adapters)
			// firstStartPosInc = start;
			// lastStartPosExc = start+regionLength-kmerLength+1;
			// }
			firstStartPosInc = start - overlapAdaptersLen;
			lastStartPosExc = start + regionLength - kmerLength + 1 + overlapAdaptersLen;

			SymbolTokenization symToken = alphabet.getTokenization("token");
			Pattern regexPattern = null;
			if (regexString != null) {
				regexPattern = Pattern.compile(regexString);
			}

			// Parse the table and populate the kmerToIntCountMap
			String countsTableBaseName = FileTools.stripPath(filePathName);
			System.err.print("Parsing " + countsTableBaseName + " for kmers of length " + kmerLength + " and populating the kmerToCount HashMap....");
			String lineString;
			boolean firstLine = true;

			int lineCounter = 0;
			while ((lineString = bufferedReader.readLine()) != null) {

				// Update parsing percentage
				lineCounter++;
				if ((lineCounter >= nextNumberThreshold) && (nextPercentThreshold <= 100)) {
					System.err.print("\n\t" + nextPercentThreshold + "% complete...");
					// System.err.print( "\n\tHi!");
					nextPercentThreshold += percentStepSize;
					nextNumberThreshold += numberStepSize;
				}

				// /////////////////////////////////////////////////////////////////////////////////////////////
				// /////////////////////////////////////////////////////////////////////////////////////////////
				// /////////////////////////////////////////////////////////////////////////////////////////////
				// /////////////////////////////////////////////////////////////////////////////////////////////
				// /////////////////////////////////////////////////////////////////////////////////////////////
				// /////////////////////////////////////////////////////////////////////////////////////////////

				// if (nextPercentThreshold == 60) {
				// break;
				// }

				// /////////////////////////////////////////////////////////////////////////////////////////////
				// /////////////////////////////////////////////////////////////////////////////////////////////
				// /////////////////////////////////////////////////////////////////////////////////////////////
				// /////////////////////////////////////////////////////////////////////////////////////////////
				// /////////////////////////////////////////////////////////////////////////////////////////////
				// /////////////////////////////////////////////////////////////////////////////////////////////

				// skip the R header
				if (skipHeader && firstLine) {
					firstLine = false;
					continue;
				}

				// String lineEntries[] = lineString.split("\t");
				String lineEntries[] = tabsPattern.split(lineString);
				if (lineEntries.length > 0) {

					String seqUpperCase = lineEntries[0].toUpperCase();

					Double seqReadCount = null;
					if (lineEntries.length == 2) {
						// seqReadCount = Integer.parseInt(lineEntries[1]);
						// seqReadCount =
						// Double.parseDouble(lineEntries[1]).intValue();
						seqReadCount = new Double(lineEntries[1]);
					}
                    else if (lineEntries.length >= 5) {
						// seqReadCount = Integer.parseInt(lineEntries[4]);
						// seqReadCount =
						// Double.parseDouble(lineEntries[1]).intValue();
						seqReadCount = new Double(lineEntries[4]);
					}
                    else {
                        seqReadCount = 1.0;
                    }

					if (!kmerLengthSet) {
						if (MathTools.isInteger(seqReadCount.doubleValue())) { // They
                            // are
                            // integers!
							useKmerToIntCountMap = true;

							// set the kmer length in (key = -1)
							if (!kmerToIntCountMap.containsKey(-1)) {
								kmerToIntCountMap.put(-1L, kmerLength);
							}
						}
                        else { // They are doubles!!!!
							useKmerToIntCountMap = false;

							// set the kmer length in (key = -1)
							if (!kmerToDoubleSumMap.containsKey(-1)) {
								kmerToDoubleSumMap.put(-1, kmerLength);
							}
						}
						kmerLengthSet = true;
					}

					if (seqReadCount.doubleValue() < 0) {
						System.err.println("\nERROR: seqReadCount = " + seqReadCount);
					}

					// System.err.println( "\nHi!");

					// /////////////////////////////////////////////////////////////////////////////////////////////
					// include each substring that matches the regex
					// /////////////////////////////////////////////////////////////////////////////////////////////
					if (regexPattern != null) {
						// System.err.println( "\nHi!");
						String origRegionString = seqUpperCase.substring(
                            firstStartPosInc, lastStartPosExc + kmerLength);
						String[] regionStrings = null;
						if (combineRC) {
							regionStrings = new String[2];
							SymbolList symList = DNATools.createDNA(origRegionString);
							SymbolList revCompSymList = SymbolListTools.reverseComplement(symList, complementTable);
							String revCompString = symToken.tokenizeSymbolList(revCompSymList);
							regionStrings[1] = revCompString.toUpperCase();
							// System.err.println("revCompString="+revCompString);
						}
                        else {
							regionStrings = new String[1];
						}
						regionStrings[0] = origRegionString;

						for (String regionString : regionStrings) {
							Matcher matcher = regexPattern.matcher(regionString);
							while (matcher.find()) {
								String matchedString = matcher.group();
								long wordLong = getWordLong(matchedString);
								long wordLongRC = getRevComp(wordLong, matchedString.length());

								if (useKmerToIntCountMap) {
									updateHashMap(kmerToIntCountMap,
                                        wordLong,
                                        wordLongRC,
                                        seqReadCount.intValue(),
                                        combineRC);
								}
                                else {
									updateHashMap(kmerToDoubleSumMap,
                                        wordLong,
                                        wordLongRC,
                                        seqReadCount.doubleValue(),
                                        combineRC);
								}

								countsTotal += seqReadCount.intValue();
							}
						}
					}
					// /////////////////////////////////////////////////////////////////////////////////////////////
					// include the kmer from every start pos
					// /////////////////////////////////////////////////////////////////////////////////////////////
					else {
						TLongHashSet kmerHashSet = null;
						if (mult == false) {
							kmerHashSet = new TLongHashSet();
						}
						// System.err.println( "\nHi5!");
						// System.err.print(
						// "\nfirstStartPosInc = "+firstStartPosInc);
						// System.err.print(
						// "\nlastStartPosExc = "+lastStartPosExc);

						// System.err.println( "\nHi!");
						for (int startPos = firstStartPosInc; startPos < lastStartPosExc; startPos++) {

							// System.err.println( "\nHi4!");
							// System.err.flush();

							String currentKmer = seqUpperCase.substring(
                                startPos, (startPos + kmerLength));
							// System.err.println( currentKmer );

							long wordLong = getWordLong(currentKmer);
							long wordLongRC = getRevComp(wordLong, kmerLength);

							// check the kmerHashSet to see if we've already
							// seen this kmer for this read
							if (mult == false) {
								if ((kmerHashSet.contains(wordLong))) {
									continue;
								}
								if (combineRC
                                    && kmerHashSet.contains(wordLongRC)) {
									continue;
								}
								kmerHashSet.add(wordLong);
							}

							if (useKmerToIntCountMap) {
								updateHashMap(kmerToIntCountMap, wordLong,
                                    wordLongRC, seqReadCount.intValue(),
                                    combineRC);
							}
                            else {
								updateHashMap(kmerToDoubleSumMap, wordLong,
                                    wordLongRC, seqReadCount.doubleValue(),
                                    combineRC);
							}

							countsTotal += seqReadCount.intValue();

							// System.err.println( "\nHi2!");

						}
					}
				}
			}
			System.err.println("\nDone.");

			// purge the kmerToIntCountMap if required
			if (countMin > 1) {
				System.err.print("Purging the kmerToCount HashMap of any entries with count < " + countMin + "...");

				if (useKmerToIntCountMap) {
					for (TLongIntIterator iter = kmerToIntCountMap.iterator(); iter.hasNext();) {
						iter.advance();

						// SKIP THE KEY THAT GIVES THE kmerLength
						if (iter.key() == -1L) {
							continue;
						}

						int count = iter.value();
						if (count < countMin) {
							// countsTotal -= count;
							iter.remove();// avoids
                            // ConcurrentModificationException
						}
					}
				}
                else {
					for (TLongDoubleIterator iter = kmerToDoubleSumMap.iterator(); iter.hasNext();) {
						iter.advance();

						// SKIP THE KEY THAT GIVES THE kmerLength
						if (iter.key() == -1L) {
							continue;
						}

						double count = iter.value();
						if (count < countMin) {
							// countsTotal -= count;
							iter.remove();// avoids
                            // ConcurrentModificationException
						}
					}
				}

				System.err.println("Done.");
			}

		} catch (Exception ex) {
			ex.printStackTrace();
		}
		// System.err.println("\ncountsTotal = "+countsTotal);
		return (countsTotal);
	}

	public void getExpectedCounts(TLongDoubleHashMap kmerToProbCountMap,
        Table outputTable, long totalCounts) {
		if (outputTable == null) {
			for (TLongDoubleIterator iter = kmerToProbCountMap.iterator(); iter.hasNext();) {
				iter.advance();

				// SKIP THE KEY THAT GIVES THE kmerLength
				if (iter.key() == -1L) {
					continue;
				}

				iter.setValue(iter.value() * totalCounts);
			}
		}
        else {
			int columns = outputTable.columns();
			outputTable.scaleColumnBy(columns - 1, totalCounts);
		}
	}

	public String getCountStats(TLongIntHashMap kmerToIntCountMap,
        int kmerLength) {
		int upperBoundKmer = multBy4ToPower(1, kmerLength);

		int ltoet0 = 0;
		int ltoet1 = 0;
		int ltoet2 = 0;
		int ltoet4 = 0;
		int ltoet6 = 0;
		int ltoet8 = 0;
		int ltoet10 = 0;
		int ltoet12 = 0;

		for (int kmer = 0; kmer < upperBoundKmer; kmer++) {
			int count = kmerToIntCountMap.get(kmer);

			if (count <= 12) {
				ltoet12++;
				if (count <= 10) {
					ltoet10++;
					if (count <= 8) {
						ltoet8++;
						if (count <= 6) {
							ltoet6++;
							if (count <= 4) {
								ltoet4++;
								if (count <= 2) {
									ltoet2++;
									if (count <= 1) {
										ltoet1++;
										if (count <= 0) {
											ltoet0++;
										}
									}
								}
							}
						}
					}
				}
			}
		}

		StringBuffer stringBuffer = new StringBuffer();

		stringBuffer.append("\n" + ltoet12 + " " + kmerLength + "-mers found on <= 12 probes the positive strand.");
		stringBuffer.append("\n" + ltoet10 + " " + kmerLength + "-mers found on <= 10 probes the positive strand.");
		stringBuffer.append("\n" + ltoet8 + " " + kmerLength + "-mers found on <= 8 probes the positive strand.");
		stringBuffer.append("\n" + ltoet6 + " " + kmerLength + "-mers found on <= 6 probes the positive strand.");
		stringBuffer.append("\n" + ltoet4 + " " + kmerLength + "-mers found on <= 4 probes the positive strand.");
		stringBuffer.append("\n" + ltoet2 + " " + kmerLength + "-mers found on <= 2 probes the positive strand.");
		stringBuffer.append("\n" + ltoet1 + " " + kmerLength + "-mers found on <= 1 probe the positive strand.");
		stringBuffer.append("\n" + ltoet0 + " " + kmerLength + "-mers not found anywhere on the positive strand.");

		return (stringBuffer.toString());
	}

	public Table getSortedKmerCountsForFile(TLongIntHashMap kmerToIntCountMap,
        String sortByFile) {
		// Table kmerCountTable = new Table(kmerToIntCountMap);
		Table kmerCountTable = new Table();

		try {
			String aLine;

			BufferedReader inBuffer = null;
			if (sortByFile.endsWith(".gz")) {
				inBuffer = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(sortByFile))));
			}
            else {
				inBuffer = new BufferedReader(new FileReader(sortByFile));
			}

			String commentToken = null;

			// Populate a Table with the HashMap entries
			System.err.println("size(kmerToIntCountMap) = " + kmerToIntCountMap.size());
			long[] keys = kmerToIntCountMap.keys();
			System.err.println("Kmers in kmerToIntCountMap are "
                + getWordString(keys[0], 16) + ", "
                + getWordString(keys[1], 16) + ", "
                + getWordString(keys[2], 16) + ", "
                + getWordString(keys[3], 16) + "....");

			System.err.print("Populating a Table with the entries found in "
                + sortByFile + "...");

			// skip the first line
			inBuffer.readLine();

			while (((aLine = inBuffer.readLine()) != null)) {
				// include this line if line doesn't start with commentToken
				if ((commentToken == null) || !aLine.trim().startsWith(commentToken)) {
					// returnBuffer.append(aLine+"\n");
					String lineEntries[] = tabsPattern.split(aLine);
					long wordLong = getWordLong(lineEntries[0]);
					int seqReadCount = kmerToIntCountMap.get(wordLong);
					kmerCountTable.add(Arrays.asList(new Long(wordLong), new Double((double) seqReadCount)));
				}
			}
			inBuffer.close();
			System.err.println("Done.");
		} catch (IOException e) {
			e.printStackTrace();
		}

		return (kmerCountTable);
	}

	public Table getSortedKmerCounts(TLongIntHashMap kmerToIntCountMap,
        String sortBy) {
		// Populate a Table with the HashMap entries
		System.err.print("Populating a Table with the entries left in the HashMap...");

		// Table kmerCountTable = new Table(kmerToIntCountMap);
		Table kmerCountTable = new Table();
		for (TLongIntIterator iter = kmerToIntCountMap.iterator(); iter.hasNext();) {
			iter.advance();
			kmerCountTable.add(Arrays.asList(new Long(iter.key()), new Double((double) iter.value())));
		}

		System.err.println("Done.");

		// Sort if necessary
		if (sortBy != null) {
			System.err.print("Sorting the Table by " + sortBy + "...");

			int sortColumn;
			if (sortBy.equalsIgnoreCase("counts")) {
				sortColumn = 1;
			}
            else {
				sortColumn = 0;
			}

			kmerCountTable.sort(sortColumn);
			kmerCountTable.reverse();
			System.err.println("Done.");
		}

		return (kmerCountTable);
	}

	public Table getSortedKmerCounts(TLongDoubleHashMap kmerToIntCountMap, String sortBy) {
		// Populate a Table with the HashMap entries
		System.err.print("Populating a Table with the entries left in the HashMap...");

		// Table kmerCountTable = new Table(kmerToIntCountMap);
		Table kmerCountTable = new Table();
		for (TLongDoubleIterator iter = kmerToIntCountMap.iterator(); iter.hasNext();) {
			iter.advance();
			kmerCountTable.add(Arrays.asList(new Long(iter.key()), new Double(iter.value())));
		}

		System.err.println("Done.");

		// Sort if necessary
		if (sortBy != null) {
			System.err.print("Sorting the Table by " + sortBy + "...");
			int sortColumn;
			if (sortBy.equalsIgnoreCase("counts")) {
				sortColumn = 1;
			}
            else {
				sortColumn = 0;
			}

			kmerCountTable.sort(sortColumn);
			kmerCountTable.reverse();
			System.err.println("Done.");
		}

		return (kmerCountTable);
	}

	public void saveTableToText(String aFilePathName, boolean append,
        Table kmerCountTable, int kmerLength, String protein, String label,
        String round, String aLeftAdapterSeq, String aRightAdapterSeq) {
		try {
			if (aLeftAdapterSeq == null) {
				aLeftAdapterSeq = "";
			}
			if (aRightAdapterSeq == null) {
				aRightAdapterSeq = "";
			}

			BufferedWriter outBuffer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(aFilePathName, append)));

			System.err.print("Saving the counts to a " + aFilePathName + "...");
			// System.err.println( "\n\tkmerLength="+kmerLength);

			if (protein != null) {
				outBuffer.write("Kmer" + "\t" + "Protein" + "\t" + "Label" + "\t" + "Round" + "\t" + "Count");
			}

			for (Iterator rowIter = kmerCountTable.iterator(); rowIter.hasNext();) {
				List currentRow = (List) rowIter.next();
				// String kmer = (String)currentRow.get(0);

				long kmerLong = ((Long) currentRow.get(0)).longValue();

				// SKIP THE KEY THAT GIVES THE kmerLength
				if (kmerLong == -1L) {
					continue;
				}

				String kmerString = getWordString(kmerLong, kmerLength);
				Double count = (Double) currentRow.get(1);

				// System.err.println("kmerString="+kmerString);
				// System.err.println("kmerLong="+kmerLong);

				String countString = null;
				if (MathTools.isInteger(count.doubleValue())) {
					countString = Integer.toString(count.intValue());
				}
                else {
					countString = count.toString();
				}
				// countString = count.toString();

				if (kmerLong < 0) {
					System.err.println("\nERROR: kmerString=" + kmerString + ", kmerLong=" + kmerLong);
				}

				if (protein != null) {
					outBuffer.write("\n" + aLeftAdapterSeq + kmerString + aRightAdapterSeq + "\t" + protein + "\t" + label + "\t" + round + "\t" + countString);
					// outBuffer.write("\n"+aLeftAdapterSeq+kmerString+aRightAdapterSeq
					// +"\t"+ countString);
					// outBuffer.write("\n"+kmerString+"\t"+ countString);
				}
                else {
					outBuffer.write("\n" + aLeftAdapterSeq + kmerString + aRightAdapterSeq + "\t" + countString);
					// outBuffer.write("\n"+kmerString+"\t"+ countString);
				}
			}
			System.err.println("Done.");
			outBuffer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public String toString(Table kmerCountTable, int kmerLength,
        String protein, String label, String round, String aLeftAdapterSeq,
        String aRightAdapterSeq) {
		if (aLeftAdapterSeq == null) {
			aLeftAdapterSeq = "";
		}
		if (aRightAdapterSeq == null) {
			aRightAdapterSeq = "";
		}

		StringBuffer stringBuffer = new StringBuffer();
		System.err.print("Printing the counts to a StringBuffer...");
		System.err.println("\n\tkmerLength=" + kmerLength);

		if (protein != null) {
			stringBuffer.append("Kmer" + "\t" + "Protein" + "\t" + "Label" + "\t" + "Round" + "\t" + "Count");
		}

		for (Iterator rowIter = kmerCountTable.iterator(); rowIter.hasNext();) {
			List currentRow = (List) rowIter.next();
			// String kmer = (String)currentRow.get(0);

			long kmerLong = ((Long) currentRow.get(0)).longValue();

			// SKIP THE KEY THAT GIVES THE kmerLength
			if (kmerLong == -1L) {
				continue;
			}

			String kmerString = getWordString(kmerLong, kmerLength);
			Double count = (Double) currentRow.get(1);

			// System.err.println("kmerString="+kmerString);
			// System.err.println("kmerLong="+kmerLong);

			String countString = null;
			if (MathTools.isInteger(count)) {
				countString = Integer.toString(count.intValue());
			}
            else {
				countString = count.toString();
			}

			if (kmerLong < 0) {
				System.err.println("\nERROR: kmerString=" + kmerString + ", kmerLong=" + kmerLong);
			}

			if (protein != null) {
				stringBuffer.append("\n" + aLeftAdapterSeq + kmerString + aRightAdapterSeq + "\t" + protein + "\t" + label + "\t" + round + "\t" + countString);
				// stringBuffer.append("\n"+aLeftAdapterSeq+kmerString+aRightAdapterSeq
				// +"\t"+ countString);
				// stringBuffer.append("\n"+kmerString+"\t"+ countString);
			}
            else {
				stringBuffer.append("\n" + aLeftAdapterSeq + kmerString + aRightAdapterSeq + "\t" + countString);
				// stringBuffer.append("\n"+kmerString+"\t"+ countString);
			}
		}
		System.err.println("Done.");
		return (stringBuffer.toString());
	}

	public String toString(Table kmerCountTable, int kmerLength,
        String protein, String round, String aLeftAdapterSeq,
        String aRightAdapterSeq) {
		StringBuffer stringBuffer = new StringBuffer();
		System.err.print("Printing the counts to a StringBuffer...");

		if (protein != null) {
			stringBuffer.append("Kmer" + "\t" + "Protein" + "\t" + "Round" + "\t" + "Count");
		}

		for (Iterator rowIter = kmerCountTable.iterator(); rowIter.hasNext();) {
			List currentRow = (List) rowIter.next();

			long kmerLong = ((Long) currentRow.get(0)).longValue();

			// SKIP THE KEY THAT GIVES THE kmerLength
			if (kmerLong == -1L) {
				continue;
			}

			String kmerString = getWordString(kmerLong, kmerLength);
			Double count = (Double) currentRow.get(1);

			String countString = null;
			if (MathTools.isInteger(count)) {
				countString = Integer.toString(count.intValue());
			}
            else {
				countString = count.toString();
			}

			if (protein != null) {
				stringBuffer.append("\n" + aLeftAdapterSeq + kmerString + aRightAdapterSeq + "\t" + protein + "\t" + round + "\t" + countString);
			}
            else {
				stringBuffer.append("\n" + aLeftAdapterSeq + kmerString + aRightAdapterSeq + "\t" + countString);
			}
		}
		System.err.println("Done.");
		return (stringBuffer.toString());
	}

	public void updateHashMap(TLongIntHashMap kmerToIntCountMap, long wordLong,
        long wordLongRC, int seqReadCount, boolean combineRC) {
		// store into HashMap
		// Integer oldCount = (Integer)kmerToIntCountMap.get(currentKmer);
		int oldCount = (int) kmerToIntCountMap.get(wordLong);

		if (oldCount == 0 && combineRC) {
			// int wordLongRC = getRevCompWord(wordLong, kmerLength);
			int oldCountRC = (int) kmerToIntCountMap.get(wordLongRC);
			if (oldCountRC == 0) {
				kmerToIntCountMap.put(wordLong, seqReadCount);
			}
            else {
				kmerToIntCountMap.put(wordLongRC, oldCountRC + seqReadCount);
			}
		}
        else {
			kmerToIntCountMap.put(wordLong, oldCount + seqReadCount);
		}
	}

	public void updateHashMap(TLongDoubleHashMap kmerToDoubleSumMap,
        long wordLong, long wordLongRC, double seqReadCount,
        boolean combineRC) {
		// store into HashMap
		// Integer oldCount = (Integer)kmerToDoubleSumMap.get(currentKmer);
		double oldCount = (double) kmerToDoubleSumMap.get(wordLong);

		if (oldCount == 0 && combineRC) {
			// int wordLongRC = getRevCompWord(wordLong, kmerLength);
			double oldCountRC = (double) kmerToDoubleSumMap.get(wordLongRC);
			if (oldCountRC == 0) {
				kmerToDoubleSumMap.put(wordLong, seqReadCount);
			}
            else {
				kmerToDoubleSumMap.put(wordLongRC, oldCountRC + seqReadCount);
			}
		}
        else {
			kmerToDoubleSumMap.put(wordLong, oldCount + seqReadCount);
		}
	}

	// word.length == arbitrary
	// "static final" is to strongly hint to the JVM to inline expand this
	// function
	public static final BigInteger multBy4ToPower(BigInteger scaler,
        int exponent) {
		// Here's the safer and maybe slower way if the javac compiler cannot
		// fully optimize
		// return(scaler * (int)Math.pow(4, exponent));

		// A potentially faster way using bit shifts
		// bit shifting left is equivalent to multiplication by powers of 2
		// bit shifting to the right is equivalent to division by powers of 2
		// x << y <===> x * (2^y)
		// x * (4^y) <==> x * (2^(2y)) <==> x << 2y
		return (scaler.shiftLeft(2 * exponent));
	}

	// word.length <= 31
	// "static final" is to strongly hint to the JVM to inline expand this
	// function
	public static final long multBy4ToPower(long scaler, int exponent) {
		// Here's the safer and maybe slower way if the javac compiler cannot
		// fully optimize
		// return(scaler * (int)Math.pow(4, exponent));

		// A potentially faster way using bit shifts
		// bit shifting left is equivalent to multiplication by powers of 2
		// bit shifting to the right is equivalent to division by powers of 2
		// x << y <===> x * (2^y)
		// x * (4^y) <==> x * (2^(2y)) <==> x << 2y
		return (scaler << (2 * exponent));
	}

	// word.length <= 15
	// "static final" is to strongly hint to the JVM to inline expand this
	// function
	public static final int multBy4ToPower(int scaler, int exponent) {
		// Here's the safer and maybe slower way if the javac compiler cannot
		// fully optimize
		// return(scaler * (int)Math.pow(4, exponent));

		// A potentially faster way using bit shifts
		// bit shifting left is equivalent to multiplication by powers of 2
		// bit shifting to the right is equivalent to division by powers of 2
		// x << y <===> x * (2^y)
		// x * (4^y) <==> x * (2^(2y)) <==> x << 2y
		return (scaler << (2 * exponent));
	}

    public static final char[] DNA_ALPHA = {'A','C','G','T'};

    public static final Symbol[] DNA_SYMBOL_ALPHA = {DNATools.a(), DNATools.c(), DNATools.g(), DNATools.t()};

    public static final int[] charToInt = {
        0, //  0  NUL (null)
        0, //  1  SOH (start of heading)
        0, //  2  STX (start of text)
        0, //  3  ETX (end of text)
        0, //  4  EOT (end of transmission)
        0, //  5  ENQ (enquiry)
        0, //  6  ACK (acknowledge)
        0, //  7  BEL (bell)
        0, //  8  BS  (backspace)
        0, //  9  TAB (horizontal tab)
        0, // 10  LF  (NL line feed, new line)
        0, // 11  VT  (vertical tab)
        0, // 12  FF  (NP form feed, new page)
        0, // 13  CR  (carriage return)
        0, // 14  SO  (shift out)
        0, // 15  SI  (shift in)
        0, // 16  DLE (data link escape)
        0, // 17  DC1 (device control 1)
        0, // 18  DC2 (device control 2)
        0, // 19  DC3 (device control 3)
        0, // 20  DC4 (device control 4)
        0, // 21  NAK (negative acknowledge)
        0, // 22  SYN (synchronous idle)
        0, // 23  ETB (end of trans. block)
        0, // 24  CAN (cancel)
        0, // 25  EM  (end of medium)
        0, // 26  SUB (substitute)
        0, // 27  ESC (escape)
        0, // 28  FS  (file separator)
        0, // 29  GS  (group separator)
        0, // 30  RS  (record separator)
        0, // 31  US  (unit separator)
        0, // 32  SPACE
        0, // 33  !
        0, // 34  "
        0, // 35  #
        0, // 36  $
        0, // 37  %
        0, // 38  &
        0, // 39  '
        0, // 40  (
        0, // 41  )
        0, // 42  *
        0, // 43  +
        0, // 44  ,
        0, // 45  -
        0, // 46  .
        0, // 47  /
        0, // 48  0
        0, // 49  1
        0, // 50  2
        0, // 51  3
        0, // 52  4
        0, // 53  5
        0, // 54  6
        0, // 55  7
        0, // 56  8
        0, // 57  9
        0, // 58  :
        0, // 59  ;
        0, // 60  <
        0, // 61  =
        0, // 62  >
        0, // 63  ?
        0, // 64  @
        0, // 65  A
        0, // 66  B
        1, // 67  C
        0, // 68  D
        0, // 69  E
        0, // 70  F
        2, // 71  G
        0, // 72  H
        0, // 73  I
        0, // 74  J
        0, // 75  K
        0, // 76  L
        0, // 77  M
        0, // 78  N
        0, // 79  O
        0, // 80  P
        0, // 81  Q
        0, // 82  R
        0, // 83  S
        3, // 84  T
        0, // 85  U
        0, // 86  V
        0, // 87  W
        0, // 88  X
        0, // 89  Y
        0, // 90  Z
        0, // 91  [
        0, // 92  \
        0, // 93  ]
        0, // 94  ^
        0, // 95  _
        0, // 96  `
        0, // 97  a
        0, // 98  b
        1, // 99  c
        0, // 100  d
        0, // 101  e
        0, // 102  f
        2, // 103  g
        0, // 104  h
        0, // 105  i
        0, // 106  j
        0, // 107  k
        0, // 108  l
        0, // 109  m
        0, // 110  n
        0, // 111  o
        0, // 112  p
        0, // 113  q
        0, // 114  r
        0, // 115  s
        3, // 116  t
        0, // 117  u
        0, // 118  v
        0, // 119  w
        0, // 120  x
        0, // 121  y
        0, // 122  z
        0, // 123  {
        0, // 124  |
        0, // 125  }
        0, // 126  ~
        0  // 127  DEL
    };

    public static final long[] charToLong = {
        0L, //  0  NUL (null)
        0L, //  1  SOH (start of heading)
        0L, //  2  STX (start of text)
        0L, //  3  ETX (end of text)
        0L, //  4  EOT (end of transmission)
        0L, //  5  ENQ (enquiry)
        0L, //  6  ACK (acknowledge)
        0L, //  7  BEL (bell)
        0L, //  8  BS  (backspace)
        0L, //  9  TAB (horizontal tab)
        0L, // 10  LF  (NL line feed, new line)
        0L, // 11  VT  (vertical tab)
        0L, // 12  FF  (NP form feed, new page)
        0L, // 13  CR  (carriage return)
        0L, // 14  SO  (shift out)
        0L, // 15  SI  (shift in)
        0L, // 16  DLE (data link escape)
        0L, // 17  DC1 (device control 1)
        0L, // 18  DC2 (device control 2)
        0L, // 19  DC3 (device control 3)
        0L, // 20  DC4 (device control 4)
        0L, // 21  NAK (negative acknowledge)
        0L, // 22  SYN (synchronous idle)
        0L, // 23  ETB (end of trans. block)
        0L, // 24  CAN (cancel)
        0L, // 25  EM  (end of medium)
        0L, // 26  SUB (substitute)
        0L, // 27  ESC (escape)
        0L, // 28  FS  (file separator)
        0L, // 29  GS  (group separator)
        0L, // 30  RS  (record separator)
        0L, // 31  US  (unit separator)
        0L, // 32  SPACE
        0L, // 33  !
        0L, // 34  "
        0L, // 35  #
        0L, // 36  $
        0L, // 37  %
        0L, // 38  &
        0L, // 39  '
        0L, // 40  (
        0L, // 41  )
        0L, // 42  *
        0L, // 43  +
        0L, // 44  ,
        0L, // 45  -
        0L, // 46  .
        0L, // 47  /
        0L, // 48  0
        0L, // 49  1
        0L, // 50  2
        0L, // 51  3
        0L, // 52  4
        0L, // 53  5
        0L, // 54  6
        0L, // 55  7
        0L, // 56  8
        0L, // 57  9
        0L, // 58  :
        0L, // 59  ;
        0L, // 60  <
        0L, // 61  =
        0L, // 62  >
        0L, // 63  ?
        0L, // 64  @
        0L, // 65  A
        0L, // 66  B
        1L, // 67  C
        0L, // 68  D
        0L, // 69  E
        0L, // 70  F
        2L, // 71  G
        0L, // 72  H
        0L, // 73  I
        0L, // 74  J
        0L, // 75  K
        0L, // 76  L
        0L, // 77  M
        0L, // 78  N
        0L, // 79  O
        0L, // 80  P
        0L, // 81  Q
        0L, // 82  R
        0L, // 83  S
        3L, // 84  T
        0L, // 85  U
        0L, // 86  V
        0L, // 87  W
        0L, // 88  X
        0L, // 89  Y
        0L, // 90  Z
        0L, // 91  [
        0L, // 92  \
        0L, // 93  ]
        0L, // 94  ^
        0L, // 95  _
        0L, // 96  `
        0L, // 97  a
        0L, // 98  b
        1L, // 99  c
        0L, // 100  d
        0L, // 101  e
        0L, // 102  f
        2L, // 103  g
        0L, // 104  h
        0L, // 105  i
        0L, // 106  j
        0L, // 107  k
        0L, // 108  l
        0L, // 109  m
        0L, // 110  n
        0L, // 111  o
        0L, // 112  p
        0L, // 113  q
        0L, // 114  r
        0L, // 115  s
        3L, // 116  t
        0L, // 117  u
        0L, // 118  v
        0L, // 119  w
        0L, // 120  x
        0L, // 121  y
        0L, // 122  z
        0L, // 123  {
        0L, // 124  |
        0L, // 125  }
        0L, // 126  ~
        0L  // 127  DEL
    };

	// word.length <= 15
	public static final int getWordInt(String word) {
		int k = 0;
		int len = word.length();

		for (int i = 0; i < len; i++) {


            // a lookup table is marginally faster than a switch (as long as the lookup table is in memory)
            k += multBy4ToPower(charToInt[word.charAt(i)], len - i - 1);


            // switch (word.charAt(i)) {
			// case 'A':
			// case 'a':
			// 	// k += 0;
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// case 'C':
			// case 'c':
			// 	// k += 1 * Math.pow(4,(len-i-1));
			// 	k += multBy4ToPower(1, len - i - 1);
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// case 'G':
			// case 'g':
			// 	// k += 2 * Math.pow(4,(len-i-1));
			// 	k += multBy4ToPower(2, len - i - 1);
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// case 'T':
			// case 't':
			// 	// k += 3 * Math.pow(4,(len-i-1));
			// 	k += multBy4ToPower(3, len - i - 1);
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// }

		}
		return (k);
	}

	// word.length <= 31
	public static final long getWordLong(String word) {
		long k = 0;
		int len = word.length();

		for (int i = 0; i < len; i++) {

            // a lookup table is marginally faster than a switch (as long as the lookup table is in memory)
            k += multBy4ToPower((long) (charToLong[word.charAt(i)]), len - i - 1);


            // switch (word.charAt(i)) {
			// case 'A':
			// case 'a':
			// 	// k += 0;
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// case 'C':
			// case 'c':
			// 	// k += 1 * Math.pow(4,(len-i-1));
			// 	k += multBy4ToPower(1L, len - i - 1);
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// case 'G':
			// case 'g':
			// 	// k += 2 * Math.pow(4,(len-i-1));
			// 	k += multBy4ToPower(2L, len - i - 1);
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// case 'T':
			// case 't':
			// 	// k += 3 * Math.pow(4,(len-i-1));
			// 	k += multBy4ToPower(3L, len - i - 1);
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// }

			// System.err.println("k="+k);
		}
		return (k);
	}

	// word.length <= 15
	public static final int getWordInt(char[] word) {
		int k = 0;
		int len = word.length;

		for (int i = 0; i < len; i++) {

            // a lookup table is marginally faster than a switch (as long as the lookup table is in memory)
            k += multBy4ToPower(charToInt[word[i]], len - i - 1);

            // switch (word[i]) {
			// case 'A':
			// case 'a':
			// 	// k += 0;
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// case 'C':
			// case 'c':
			// 	// k += 1 * Math.pow(4,(len-i-1));
			// 	k += multBy4ToPower(1, len - i - 1);
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// case 'G':
			// case 'g':
			// 	// k += 2 * Math.pow(4,(len-i-1));
			// 	k += multBy4ToPower(2, len - i - 1);
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// case 'T':
			// case 't':
			// 	// k += 3 * Math.pow(4,(len-i-1));
			// 	k += multBy4ToPower(3, len - i - 1);
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// }
		}
		return (k);
	}

	// word.length <= 31
	public static final long getWordLong(char[] word) {
		long k = 0;
		int len = word.length;

		for (int i = 0; i < len; i++) {

            // a lookup table is marginally faster than a switch (as long as the lookup table is in memory)
            k += multBy4ToPower((long) (charToLong[word[i]]), len - i - 1);


            // switch (word[i]) {
			// case 'A':
			// case 'a':
			// 	// k += 0;
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// case 'C':
			// case 'c':
			// 	// k += 1 * Math.pow(4,(len-i-1));
			// 	k += multBy4ToPower(1L, len - i - 1);
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// case 'G':
			// case 'g':
			// 	// k += 2 * Math.pow(4,(len-i-1));
			// 	k += multBy4ToPower(2L, len - i - 1);
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// case 'T':
			// case 't':
			// 	// k += 3 * Math.pow(4,(len-i-1));
			// 	k += multBy4ToPower(3L, len - i - 1);
			// 	// printf("maxmer[%d]=%c; k=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], k, maxWordSize);
			// 	break;
			// }
		}
		return (k);
	}

	// 0-based indexing, (start, end]
	public static final char getSubWord(char word, int len, int start, int end) {
		char reverseIndex = 0;

		for (int i = 0; i < len; i++) {
			char mod = (char) (word % 4);
			int pos = len - i - 1;

			if (pos < end) {
				reverseIndex += multBy4ToPower(mod, i - (len - end));
			}

			// fprintf(fout, "word=%d, mod=%d, reverseIndex=%d, i=%d\n", word,
			// mod, reverseIndex, i);
			word = (char) ((word - mod) / 4);

			if (pos <= start) {
				break;
			}
		}

		// fprintf(fout, "\n\n");
		return (reverseIndex);
	}

	// 0-based indexing, (start, end]
	public static final int getSubWord(int word, int len, int start, int end) {
		int reverseIndex = 0;

		for (int i = 0; i < len; i++) {
			int mod = (int) (word % 4);
			int pos = len - i - 1;

			if (pos < end) {
				reverseIndex += multBy4ToPower(mod, i - (len - end));
			}

			// fprintf(fout, "word=%d, mod=%d, reverseIndex=%d, i=%d\n", word,
			// mod, reverseIndex, i);
			word = (int) ((word - mod) / 4);

			if (pos <= start) {
				break;
			}
		}

		// fprintf(fout, "\n\n");
		return (reverseIndex);
	}

	// 0-based indexing, (start, end]
	public static final long getSubWord(long word, int len, int start, int end) {
		long reverseIndex = 0;

		for (int i = 0; i < len; i++) {
			long mod = (long) (word % 4L);
			int pos = len - i - 1;

			if (pos < end) {
				reverseIndex += multBy4ToPower(mod, i - (len - end));
			}

			// fprintf(fout, "word=%d, mod=%d, reverseIndex=%d, i=%d\n", word,
			// mod, reverseIndex, i);
			word = (long) ((word - mod) / 4L);

			if (pos <= start) {
				break;
			}
		}

		// fprintf(fout, "\n\n");
		return (reverseIndex);
	}

	// word.length <= 15
	public static final char[] getWordCharArray(int wordInt, int len) {
		char[] wordCharArray = new char[len];
		int temp = wordInt;

		for (int i = 0; i < len; i++) {
			int mod = temp % 4;

			// switch (mod) {
			// case 0:
			// 	wordCharArray[len - i - 1] = 'A';
			// 	break;
			// case 1:
			// 	wordCharArray[len - i - 1] = 'C';
			// 	break;
			// case 2:
			// 	wordCharArray[len - i - 1] = 'G';
			// 	break;
			// default:
			// 	wordCharArray[len - i - 1] = 'T';
			// 	break;
			// }

            // a lookup table is marginally faster than a switch (as long as the lookup table is in memory)
            wordCharArray[len - i - 1] = DNA_ALPHA[mod];

			temp = (temp - mod) / 4;
		}
		return (wordCharArray);
	}

	// word.length <= 31
	public static final Symbol[] getSymbolArray(long wordIndex, int len) {
		Symbol[] word = new Symbol[len];
		long temp = wordIndex;

		for (int i = 0; i < len; i++) {
			long modLong = temp % 4L;
			int mod = (int) modLong;

			// switch (mod) {
			// case 0:
			// 	word[len - i - 1] = DNATools.a();
			// 	break;
			// case 1:
			// 	word[len - i - 1] = DNATools.c();
			// 	break;
			// case 2:
			// 	word[len - i - 1] = DNATools.g();
			// 	break;
			// case 3:
			// 	word[len - i - 1] = DNATools.t();
			// 	break;
			// }

            // a lookup table is marginally faster than a switch (as long as the lookup table is in memory)
            word[len - i - 1] = DNA_SYMBOL_ALPHA[mod];

			temp = (temp - mod) / 4L;
		}
		return (word);
	}

	// word.length <= 31
	public static final char[] getWordCharArray(long wordLong, int len) {
		char[] wordCharArray = new char[len];
		long temp = wordLong;

		for (int i = 0; i < len; i++) {
			// BigInteger bigTemp = new BigInteger(Long.toString(temp));
			// mod = bigTemp.mod(new BigInteger(Long.toString(4L))).intValue();
			long modLong = (long) (temp % 4L);
			int mod = (int) modLong;

			// if ((mod < 0) || (mod > 3)) {
			// 	System.err.println("\nERROR: mod=" + mod + ", wordLong="
			// 			+ wordLong + ", len=" + len);
			// }

			// switch (mod) {
			// case 0:
			// 	wordCharArray[len - i - 1] = 'A';
			// 	break;
			// case 1:
			// 	wordCharArray[len - i - 1] = 'C';
			// 	break;
			// case 2:
			// 	wordCharArray[len - i - 1] = 'G';
			// 	break;
			// case 3:
			// 	wordCharArray[len - i - 1] = 'T';
			// 	break;
			// }

            // a lookup table is marginally faster than a switch (as long as the lookup table is in memory)
            wordCharArray[len - i - 1] = DNA_ALPHA[mod];

			temp = (temp - mod) / 4L;

			// if ((temp % 4L) != 0) {
			// System.err.println("\nERROR: temp="+temp+", temp%4L="+temp % 4L);
			// }
		}
		return (wordCharArray);
	}

	// word.length <= 15
	public static final String getWordString(int wordInt, int len) {
		return (new String(getWordCharArray(wordInt, len)));
	}

	// word.length <= 31
	public static final String getWordString(long wordLong, int len) {
		return (new String(getWordCharArray(wordLong, len)));
	}

	// word.length <= 7
	public static final char getRevComp(char word, int len) {
		char revCompIndex = 0;

		for (int i = 0; i < len; i++) {
			char mod = (char) (word % 4);

			// switch (mod) {
			// case 0: // A -> T
			// 	// revCompIndex += 3 * Math.pow(4,len-i-1);
			// 	revCompIndex += multBy4ToPower(3, len - i - 1);
			// 	break;
			// case 1: // C -> G
			// 	// revCompIndex += 2 * Math.pow(4,len-i-1);
			// 	revCompIndex += multBy4ToPower(2, len - i - 1);
			// 	break;
			// case 2: // G -> C
			// 	// revCompIndex += 1 * Math.pow(4,len-i-1);
			// 	revCompIndex += multBy4ToPower(1, len - i - 1);
			// 	break;
			// case 3: // T -> A
			// 	// revCompIndex += 0;
			// 	break;
			// }

            revCompIndex += multBy4ToPower(3 - mod, len - i - 1);

			// fprintf(fout, "word=%d, mod=%d, revCompIndex=%d, i=%d\n", word,
			// mod, revCompIndex, i);

			word = (char) ((word - mod) / 4);
		}

		// fprintf(fout, "\n\n");

		return (revCompIndex);
	}

	// word.length <= 15
	public static final int getRevComp(int word, int len) {
		int revCompIndex = 0;

		for (int i = 0; i < len; i++) {
			int mod = (int) (word % 4);

			// switch (mod) {
			// case 0: // A -> T
			// 	// revCompIndex += 3 * Math.pow(4,len-i-1);
			// 	revCompIndex += multBy4ToPower(3, len - i - 1);
			// 	break;
			// case 1: // C -> G
			// 	// revCompIndex += 2 * Math.pow(4,len-i-1);
			// 	revCompIndex += multBy4ToPower(2, len - i - 1);
			// 	break;
			// case 2: // G -> C
			// 	// revCompIndex += 1 * Math.pow(4,len-i-1);
			// 	revCompIndex += multBy4ToPower(1, len - i - 1);
			// 	break;
			// case 3: // T -> A
			// 	// revCompIndex += 0;
			// 	break;
			// }

            revCompIndex += multBy4ToPower(3 - mod, len - i - 1);

			// fprintf(fout, "word=%d, mod=%d, revCompIndex=%d, i=%d\n", word,
			// mod, revCompIndex, i);

			word = (int) ((word - mod) / 4);
		}

		// fprintf(fout, "\n\n");

		return (revCompIndex);
	}

	// word.length <= 31
	public static final long getRevComp(long word, int len) {
		long revCompIndex = 0;

		for (int i = 0; i < len; i++) {
			long modLong = (long) (word % 4L);

			// int mod = (int) modLong;
			// switch (mod) {
			// case 0: // A -> T
			// 	// revCompIndex += 3 * Math.pow(4,len-i-1);
			// 	revCompIndex += multBy4ToPower(3L, len - i - 1);
			// 	break;
			// case 1: // C -> G
			// 	// revCompIndex += 2 * Math.pow(4,len-i-1);
			// 	revCompIndex += multBy4ToPower(2L, len - i - 1);
			// 	break;
			// case 2: // G -> C
			// 	// revCompIndex += 1 * Math.pow(4,len-i-1);
			// 	revCompIndex += multBy4ToPower(1L, len - i - 1);
			// 	break;
			// case 3: // T -> A
			// 	// revCompIndex += 0;
			// 	break;
			// }

            revCompIndex += multBy4ToPower((long) (3L - modLong), len - i - 1);

			// fprintf(fout, "word=%d, mod=%d, revCompIndex=%d, i=%d\n", word,
			// mod, revCompIndex, i);

			// word = (word - mod) / 4L;
			word = (long) ((word - modLong) / 4L);
		}

		// fprintf(fout, "\n\n");

		return (revCompIndex);
	}

	// word.length <= 7
	public static final char getReverse(char word, int len) {
		char reverseIndex = 0;

		for (int i = 0; i < len; i++) {
			char mod = (char) (word % 4);

			// reverseIndex += mod * Math.pow(4,len-i-1);
			reverseIndex += multBy4ToPower(mod, len - i - 1);

			// fprintf(fout, "word=%d, mod=%d, reverseIndex=%d, i=%d\n", word,
			// mod, reverseIndex, i);

			word = (char) ((word - mod) / 4);
		}

		// fprintf(fout, "\n\n");

		return (reverseIndex);
	}

	// word.length <= 15
	public static final int getReverse(int word, int len) {
		int reverseIndex = 0;

		for (int i = 0; i < len; i++) {
			int mod = word % 4;

			// reverseIndex += mod * Math.pow(4,len-i-1);
			reverseIndex += multBy4ToPower(mod, len - i - 1);

			// fprintf(fout, "word=%d, mod=%d, reverseIndex=%d, i=%d\n", word,
			// mod, reverseIndex, i);

			word = (word - mod) / 4;
		}

		// fprintf(fout, "\n\n");

		return (reverseIndex);
	}

	// word.length <= 31
	public static final long getReverse(long word, int len) {
		long reverseIndex = 0;

		for (int i = 0; i < len; i++) {
			long mod = word % 4L;

			// reverseIndex += mod * Math.pow(4,len-i-1);
			reverseIndex += multBy4ToPower(mod, len - i - 1);

			// fprintf(fout, "word=%d, mod=%d, reverseIndex=%d, i=%d\n", word,
			// mod, reverseIndex, i);

			word = (word - mod) / 4L;
		}

		// fprintf(fout, "\n\n");

		return (reverseIndex);
	}

	// word.length <= 31
	public static final boolean isRepeat(long wordIndex, int len) {
		int mid = len / 2;

		long leftHalf = getSubWord(wordIndex, len, 0, mid);
		long rightHalf = getSubWord(wordIndex, len, mid, len);

		// System.out.println("\nword="+getWordString(wordIndex,len)+",  leftHalf="+getWordString(leftHalf,
		// mid)+",  rightHalf="+getWordString(rightHalf,mid));

		// check that they are the same length
		if ((mid - 0) != (len - mid)) {
			return (false);
		}
		return (leftHalf == rightHalf);
	}

	// word.length <= 31
	public static final boolean isSelfReverse(long wordIndex, int len) {
		long wordReverse = getReverse(wordIndex, len);
		if (wordReverse == wordIndex) {
			return (true);
		}
		return (false);
	}

	// word.length <= 31
	public static final boolean isSelfReverseComplement(long wordIndex, int len) {
		// int wordRevComp = revCompMatrix[(int)wordIndex];
		long wordRevComp = getRevComp(wordIndex, len);
		if (wordRevComp == wordIndex) {
			return (true);
		}
		return (false);
	}

	// word.length <= 31
	public static final int getContent(long word, int len, char nucleotide) {
		int nuclInt = 0;

        // a lookup table is marginally faster than a switch (as long as the lookup table is in memory)
        nuclInt = charToInt[nucleotide];

		// switch (nucleotide) {
		// case 'A':
		// case 'a':
		// 	nuclInt = 0;
		// 	break;
		// case 'C':
		// case 'c':
		// 	nuclInt = 1;
		// 	break;
		// case 'G':
		// case 'g':
		// 	nuclInt = 2;
		// 	break;
		// case 'T':
		// case 't':
		// 	nuclInt = 3;
		// 	break;
		// }

		return (getContent(word, len, nuclInt));
	}

	// word.length <= 15
	public static final int getContent(int word, int len, char nucleotide) {
		int nuclInt = 0;

        // a lookup table is marginally faster than a switch (as long as the lookup table is in memory)
        nuclInt = charToInt[nucleotide];

		// switch (nucleotide) {
		// case 'A':
		// case 'a':
		// 	nuclInt = 0;
		// 	break;
		// case 'C':
		// case 'c':
		// 	nuclInt = 1;
		// 	break;
		// case 'G':
		// case 'g':
		// 	nuclInt = 2;
		// 	break;
		// case 'T':
		// case 't':
		// 	nuclInt = 3;
		// 	break;
		// }

		return (getContent(word, len, nuclInt));
	}

	// word.length <= 15
	// A = 0
	// C = 1
	// G = 2
	// T = 3
	public static final int getContent(int word, int len, int nucleotide) {
		int mod;
		int content = 0;

		int i = 0;
		while (i < len) {
			mod = word % 4;
			if (mod == nucleotide) {
				content++;
			}

			word = (word - mod) / 4;
			i++;
		}
		return (content);
	}

	// word.length <= 31
	// A = 0
	// C = 1
	// G = 2
	// T = 3
	public static final int getContent(long word, int len, int nucleotide) {
		int content = 0;

		int i = 0;
		while (i < len) {
			// BigInteger bigTemp = new BigInteger(Long.toString(word));
			// mod = bigTemp.mod(new BigInteger(Long.toString(4L))).intValue();
			long modLong = (long) (word % 4L);
			int mod = (int) modLong;

			if (mod == nucleotide) {
				content++;
			}
			word = (long) ((word - modLong) / 4L);
			i++;
		}
		return (content);
	}

}
