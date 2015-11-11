

import java.util.Set;
import java.util.Arrays;
import java.io.*;
import java.lang.Math;
import java.lang.String;
import java.lang.StringBuffer;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.BioException;
import org.biojava.bio.mydp.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.io.SeqIOTools;
import java.util.Iterator;
import org.biojava.bio.seq.io.SymbolTokenization;
import java.util.NoSuchElementException;
import java.util.regex.Pattern;
import java.math.BigInteger;


import gnu.trove.set.hash.TLongHashSet;
import gnu.trove.set.hash.TCharHashSet;
import gnu.trove.map.hash.TObjectCharHashMap;
import gnu.trove.list.array.TCharArrayList;

import org.apache.commons.lang.mutable.MutableDouble;


public class KmerMatrix implements Serializable {

	private static final long serialVersionUID = 7923870404512774823L;

	public int kmerLength = 0;
	public int numSeqIDs = 0;
	public int numKmerWindows = 0;

	// //////////////////////////////////////////////////////////////////////////////
	// Even though we wish to use only one of the following structs at any given
	// time,
	// Java does not support "unions".
	//
	// The goal here is to have a viable data struct available for most
	// scenarios
	// //////////////////////////////////////////////////////////////////////////////

	// for (kmerLength <= 8)
	// char[seqIDNum][startPos] = Kmer
	char[][] charKmerMatrix = null;

	// for (kmerLength <= 15)
	// int[seqIDNum][startPos] = Kmer
	int[][] intKmerMatrix = null;

	// for (kmerLength <= 31)
	// long[seqIDNum][startPos] = Kmer
	long[][] longKmerMatrix = null;

	// for (kmerLength >= 32)
	// BigInteger[seqIDNum][startPos] = Kmer
	BigInteger[][] bigIntegerKmerMatrix = null;

	// //////////////////////////////////////////////////////////////////////////////
	// //////////////////////////////////////////////////////////////////////////////

	// If kmer Length <= 8 then create char lookup tables
	char[] revCompMatrix = null;
	TObjectCharHashMap symListToWordMap = null;
	SymbolList[] wordToSymListMatrix = null;
	char[][] wordToProbesMatrix = null;

	TLongHashSet allWordsSet = null;
	Table allSeqIdIndexes = null;

	private static String tabs = "\\t+";
	private static Pattern tabsPattern = Pattern.compile(tabs);

	private static String alphabetName = "DNA";
	private Alphabet alphabet;
	private SymbolTokenization symbolTokenization = null;

	// //////////////////////////////////////////////////////////////////////////////
	// //////////////////////////////////////////////////////////////////////////////

	public KmerMatrix(int aNumSeqIDs, int aKmerLength, int seqLengths) {
		this(aNumSeqIDs, aKmerLength, seqLengths, null, null, null);
	}

	public KmerMatrix(int aNumSeqIDs, int aKmerLength, int seqLengths,
			String anAlphabetName, Alphabet anAlphabet,
			SymbolTokenization aSymTok) {
		if (anAlphabetName == null) {
			setAlphabet(alphabetName);
		} else {
			setAlphabet(anAlphabetName, anAlphabet, aSymTok);
		}

		int numKmerWindows = seqLengths - aKmerLength + 1;

		this.numSeqIDs = aNumSeqIDs;
		this.kmerLength = aKmerLength;
		this.numKmerWindows = numKmerWindows;

		// ****KmerMatrix[seqIDNum][startPos] = Kmer;

		if (aKmerLength <= 8) {
			this.charKmerMatrix = new char[aNumSeqIDs][numKmerWindows];
		} else if (aKmerLength <= 15) {
			this.intKmerMatrix = new int[aNumSeqIDs][numKmerWindows];
		} else if (aKmerLength <= 31) {
			this.longKmerMatrix = new long[aNumSeqIDs][numKmerWindows];
		} else {
			// need to add BigInteger support
			System.err.println("Error : KmerMatrix.java - kmerLength "
					+ aKmerLength + " is to long.");
		}

		if (kmerLength <= 8) {
			System.out.print("Creating reverse complement lookup table...");
			revCompMatrix = makeRevCompMatrix(kmerLength);
			System.out.println("Done.");

			System.out.print("Creating Symbol List lookup table...");
			wordToSymListMatrix = makeWordToSymListMatrix(kmerLength);
			System.out.println("Done.");

			System.out.print("Creating Symbol List to Kmer hashtable...");
			symListToWordMap = makeSymListToWordMap(kmerLength);
			System.out.println("Done.");

			// System.out.print("Creating Kmer to Probes lookup table...");
			// wordToProbesMatrix = makeKmerToProbesMatrix(kmerLength);
			// System.out.println("Done.");
		}

	}

	// //////////////////////////////////////////////////////////////////////////////
	// //////////////////////////////////////////////////////////////////////////////

	// double[] getKmerToAffinityMatrix() {
	// return(kmerToAffinityMatrix);
	// }

	public char[] getRevCompMatrix() {
		return (revCompMatrix);
	}

	public TObjectCharHashMap getSymListToWordMap() {
		return (symListToWordMap);
	}

	public void setAlphabet(String anAlphabetName, Alphabet anAlphabet,
			SymbolTokenization aSymTok) {
		KmerMatrix.alphabetName = anAlphabetName;
		this.alphabet = anAlphabet;
		this.symbolTokenization = aSymTok;
	}

	public void setAlphabet(String anAlphabetName) {
		try {
			KmerMatrix.alphabetName = anAlphabetName;
			try {
				this.alphabet = AlphabetManager.alphabetForName(alphabetName);
			} catch (NoSuchElementException ex) {
				// try it upper case
				this.alphabet = AlphabetManager.alphabetForName(alphabetName
						.toUpperCase());
			}
			this.symbolTokenization = alphabet.getTokenization("token");
		} catch (Exception ex) {
			ex.printStackTrace(System.err);
		}
	}

	public void setOrderedSubset(String[] sourceSeqIDs, String[] orderedSeqIDs) {

		this.numSeqIDs = orderedSeqIDs.length;

		// Table has rows of (sedID, index)
		Table seqIDsTable = new Table(sourceSeqIDs);

		// sort by seqID
		seqIDsTable.sort(0);

		if (this.charKmerMatrix != null) {
			char[][] new_charKmerMatrix = new char[orderedSeqIDs.length][];
			for (int i = 0; i < orderedSeqIDs.length; i++) {
				String aSeqID = orderedSeqIDs[i];
				int sortedIndex = seqIDsTable.binarySearch(0, aSeqID);
				int oldIndex = (Integer) seqIDsTable
						.getRow(sortedIndex).get(1);
				new_charKmerMatrix[i] = this.charKmerMatrix[oldIndex];
			}
			this.charKmerMatrix = new_charKmerMatrix;
		} else if (this.intKmerMatrix != null) {
			int[][] new_intKmerMatrix = new int[orderedSeqIDs.length][];
			for (int i = 0; i < orderedSeqIDs.length; i++) {
				String aSeqID = orderedSeqIDs[i];
				int sortedIndex = seqIDsTable.binarySearch(0, aSeqID);
				int oldIndex = (Integer) seqIDsTable
						.getRow(sortedIndex).get(1);
				new_intKmerMatrix[i] = this.intKmerMatrix[oldIndex];
			}
			this.intKmerMatrix = new_intKmerMatrix;
		} else if (this.longKmerMatrix != null) {
			long[][] new_longKmerMatrix = new long[orderedSeqIDs.length][];
			for (int i = 0; i < orderedSeqIDs.length; i++) {
				String aSeqID = orderedSeqIDs[i];
				int sortedIndex = seqIDsTable.binarySearch(0, aSeqID);
				int oldIndex = (Integer) seqIDsTable
						.getRow(sortedIndex).get(1);
				new_longKmerMatrix[i] = this.longKmerMatrix[oldIndex];
			}
			this.longKmerMatrix = new_longKmerMatrix;
		} else if (this.bigIntegerKmerMatrix != null) {

		}

	}

	public long[] getAllWordIndexes() {
		if (allWordsSet == null) {
			allWordsSet = new TLongHashSet(getWordIndexes(this.numSeqIDs));
		}
		return (allWordsSet.toArray());
	}

	// return the wordIndexes of all kmers found in the first numberOfSequences
	// The sequences are in intensity-descending order
	public long[] getWordIndexes(int numberOfSequences) {
		TLongHashSet wordsSet = null;
		int numKmerWindows = this.numKmerWindows;

		if (this.charKmerMatrix != null) {
			int entriesTotal = Math.min(multBy4ToPower(1, this.kmerLength),
					numberOfSequences * numKmerWindows);
			allWordsSet = new TLongHashSet(entriesTotal);

			for (int aSeqIdIndex = 0; aSeqIdIndex < numberOfSequences; aSeqIdIndex++) {
				allWordsSet.addAll(ArrayTools
						.toLongArray(this.charKmerMatrix[aSeqIdIndex]));
			}
		} else if (this.intKmerMatrix != null) {
			int entriesTotal = Math.min(multBy4ToPower(1, this.kmerLength),
					numberOfSequences * numKmerWindows);
			allWordsSet = new TLongHashSet(entriesTotal);

			for (int aSeqIdIndex = 0; aSeqIdIndex < numberOfSequences; aSeqIdIndex++) {
				allWordsSet.addAll(ArrayTools
						.toLongArray(this.intKmerMatrix[aSeqIdIndex]));
			}
		} else if (this.longKmerMatrix != null) {
			// long entriesTotal = Math.min(multBy4ToPower(1L, this.kmerLength),
			// numberOfSequences*numKmerWindows);
			int entriesTotal = numberOfSequences * numKmerWindows; // can't
																	// support
																	// all
																	// possible
																	// long
																	// words!!!
			allWordsSet = new TLongHashSet(entriesTotal);

			for (int aSeqIdIndex = 0; aSeqIdIndex < numberOfSequences; aSeqIdIndex++) {
				allWordsSet.addAll(this.longKmerMatrix[aSeqIdIndex]);
			}
		}

		// BigIntegers cannot be converted to longs!

		return (allWordsSet.toArray());
	}

	public int getCount(int aSeqIdIndex, long aWordIndex,
			WeightMatrixTools.BindingStrand strand) {
		int count = 0;

		long revCompIndex = getRevComp(aWordIndex, this.kmerLength);

		if ((strand == WeightMatrixTools.BindingStrand.POS)
				|| (strand == WeightMatrixTools.BindingStrand.BOTH)) {

			if (this.charKmerMatrix != null) {
				for (int startPos = 0; startPos < this.charKmerMatrix[aSeqIdIndex].length; startPos++) {
					if (this.charKmerMatrix[aSeqIdIndex][startPos] == aWordIndex) {
						count++;
					}
				}
			} else if (this.intKmerMatrix != null) {
				for (int startPos = 0; startPos < this.intKmerMatrix[aSeqIdIndex].length; startPos++) {
					if (this.intKmerMatrix[aSeqIdIndex][startPos] == aWordIndex) {
						count++;
					}
				}
			} else if (this.longKmerMatrix != null) {
				for (int startPos = 0; startPos < this.longKmerMatrix[aSeqIdIndex].length; startPos++) {
					if (this.longKmerMatrix[aSeqIdIndex][startPos] == aWordIndex) {
						count++;
					}
				}
			}

		}

		// DON'T COUNT SELF-REV-COMPS TWICE!!
		if ((strand == WeightMatrixTools.BindingStrand.NEG)
				|| ((strand == WeightMatrixTools.BindingStrand.BOTH) && (revCompIndex != aWordIndex))) {

			if (this.charKmerMatrix != null) {
				for (int startPos = 0; startPos < this.charKmerMatrix[aSeqIdIndex].length; startPos++) {
					if (this.charKmerMatrix[aSeqIdIndex][startPos] == revCompIndex) {
						count++;
					}
				}
			} else if (this.intKmerMatrix != null) {
				for (int startPos = 0; startPos < this.intKmerMatrix[aSeqIdIndex].length; startPos++) {
					if (this.intKmerMatrix[aSeqIdIndex][startPos] == revCompIndex) {
						count++;
					}
				}
			} else if (this.longKmerMatrix != null) {
				for (int startPos = 0; startPos < this.longKmerMatrix[aSeqIdIndex].length; startPos++) {
					if (this.longKmerMatrix[aSeqIdIndex][startPos] == revCompIndex) {
						count++;
					}
				}
			}

		}
		return (count);
	}

	public int getCount(int aSeqIdIndex, long aWordIndex, int startPos,
			WeightMatrixTools.BindingStrand strand) {
		int count = 0;

		long revCompIndex = getRevComp(aWordIndex, this.kmerLength);

		if ((strand == WeightMatrixTools.BindingStrand.POS)
				|| (strand == WeightMatrixTools.BindingStrand.BOTH)) {

			if (this.charKmerMatrix != null) {
				if (this.charKmerMatrix[aSeqIdIndex][startPos] == aWordIndex) {
					count++;
				}
			} else if (this.intKmerMatrix != null) {
				if (this.intKmerMatrix[aSeqIdIndex][startPos] == aWordIndex) {
					count++;
				}
			} else if (this.longKmerMatrix != null) {
				if (this.longKmerMatrix[aSeqIdIndex][startPos] == aWordIndex) {
					count++;
				}
			}

		}

		// DON'T COUNT SELF-REV-COMPS TWICE!!
		if ((strand == WeightMatrixTools.BindingStrand.NEG)
				|| ((strand == WeightMatrixTools.BindingStrand.BOTH) && (revCompIndex != aWordIndex))) {

			if (this.charKmerMatrix != null) {
				if (this.charKmerMatrix[aSeqIdIndex][startPos] == revCompIndex) {
					count++;
				}
			} else if (this.intKmerMatrix != null) {
				if (this.intKmerMatrix[aSeqIdIndex][startPos] == revCompIndex) {
					count++;
				}
			} else if (this.longKmerMatrix != null) {
				if (this.longKmerMatrix[aSeqIdIndex][startPos] == revCompIndex) {
					count++;
				}
			}

		}
		return (count);
	}

	// weights[0] = positive strand weights
	// weights[1] = negative strand weights
	public double getCount(int aSeqIdIndex, long aWordIndex,
			double[][] weights, WeightMatrixTools.BindingStrand strand) {
		double count = 0;

		long revCompIndex = getRevComp(aWordIndex, this.kmerLength);

		if ((strand == WeightMatrixTools.BindingStrand.POS)
				|| (strand == WeightMatrixTools.BindingStrand.BOTH)) {

			if (this.charKmerMatrix != null) {
				for (int startPos = 0; startPos < this.charKmerMatrix[aSeqIdIndex].length; startPos++) {
					if (this.charKmerMatrix[aSeqIdIndex][startPos] == aWordIndex) {
						count += weights[0][startPos];
					}
				}
			} else if (this.intKmerMatrix != null) {
				for (int startPos = 0; startPos < this.intKmerMatrix[aSeqIdIndex].length; startPos++) {
					if (this.intKmerMatrix[aSeqIdIndex][startPos] == aWordIndex) {
						count += weights[0][startPos];
					}
				}
			} else if (this.longKmerMatrix != null) {
				for (int startPos = 0; startPos < this.longKmerMatrix[aSeqIdIndex].length; startPos++) {
					if (this.longKmerMatrix[aSeqIdIndex][startPos] == aWordIndex) {
						count += weights[0][startPos];
					}
				}
			}

		}

		if ((strand == WeightMatrixTools.BindingStrand.NEG)
				|| (strand == WeightMatrixTools.BindingStrand.BOTH)) {
			// ((strand == WeightMatrixTools.BindingStrand.BOTH) &&
			// (revCompIndex != aWordIndex))) {

			if (this.charKmerMatrix != null) {
				for (int startPos = 0; startPos < this.charKmerMatrix[aSeqIdIndex].length; startPos++) {
					if (this.charKmerMatrix[aSeqIdIndex][startPos] == revCompIndex) {
						count += weights[1][startPos];
					}
				}
			} else if (this.intKmerMatrix != null) {
				for (int startPos = 0; startPos < this.intKmerMatrix[aSeqIdIndex].length; startPos++) {
					if (this.intKmerMatrix[aSeqIdIndex][startPos] == revCompIndex) {
						count += weights[1][startPos];
					}
				}
			} else if (this.longKmerMatrix != null) {
				for (int startPos = 0; startPos < this.longKmerMatrix[aSeqIdIndex].length; startPos++) {
					if (this.longKmerMatrix[aSeqIdIndex][startPos] == revCompIndex) {
						count += weights[1][startPos];
					}
				}
			}

		}

		// DON'T COUNT SELF-REV-COMPS TWICE!!
		// self-reverse complements words have been counted twice IF SCORED BOTH
		// STRANDS!
		if ((revCompIndex == aWordIndex)
				&& (strand == WeightMatrixTools.BindingStrand.BOTH)) {
			count = count / 2.0;
		}

		return (count);
	}

	// weights[0] = positive strand weights
	// weights[1] = negative strand weights
	public double getCount(int aSeqIdIndex, char specialWord,
			double[][] weights,
			// WeightMatrixTools.BindingStrand strand,
			double[] oldKmerToAffinityMatrix, MutableDouble strandAffinityMD,
			double initAffinity, boolean includeRevComps) {
		double weightedCount = 0;
		double strandAffinity = 0;

		char specialWordRevComp = revCompMatrix[specialWord];

		boolean isSelfRevComp = false;
		if (specialWordRevComp == specialWord) {
			isSelfRevComp = true;
		}

		boolean oldAffinitiesIsNull = false;
		if (oldKmerToAffinityMatrix == null) {
			oldAffinitiesIsNull = true;
		}

		for (int startPos = 0; startPos < this.charKmerMatrix[aSeqIdIndex].length; startPos++) {
			char aWordIndex = this.charKmerMatrix[aSeqIdIndex][startPos];
			double affinity = 0;

			// Positive Strand
			if (aWordIndex == specialWord) {
				weightedCount += weights[0][startPos];
			} else if (!oldAffinitiesIsNull) {
				affinity += weights[0][startPos]
						* oldKmerToAffinityMatrix[aWordIndex];
			} else {
				affinity += initAffinity;
			}

			char aRevCompIndex = revCompMatrix[aWordIndex];
			if (includeRevComps) {
				// Negative Strand
				if (aRevCompIndex == specialWord) {
					// if (aWordIndex == specialWordRevComp) {
					weightedCount += weights[1][startPos];
				} else if (!oldAffinitiesIsNull) {
					// strandAffinity += weights[1][startPos] *
					// oldKmerToAffinityMatrix[aWordIndex];
					affinity += weights[1][startPos]
							* oldKmerToAffinityMatrix[aRevCompIndex];
				} else {
					affinity += initAffinity;
				}
			}

			// increase strand affinity
			if (affinity != 0) {

				// self-reverse-complement words have been counted twice
				if ((includeRevComps) && (aRevCompIndex == aWordIndex)) {
					affinity = affinity / 2.0;
				}

				strandAffinity += affinity;
			}

		}

		// self-reverse-complement words have been counted twice
		if (includeRevComps && isSelfRevComp) {
			weightedCount = weightedCount / 2.0;
		}

		strandAffinityMD.setValue(strandAffinity);
		return (weightedCount);
	}

	public boolean contains(int aSeqIdIndex, char aWordIndex,
			WeightMatrixTools.BindingStrand strand) {

		char revCompIndex = revCompMatrix[aWordIndex];

		for (int startPos = 0; startPos < this.charKmerMatrix[aSeqIdIndex].length; startPos++) {
			if (this.charKmerMatrix[aSeqIdIndex][startPos] == aWordIndex) {
				return (true);
			}
		}

		if (revCompIndex != aWordIndex) {
			for (int startPos = 0; startPos < this.charKmerMatrix[aSeqIdIndex].length; startPos++) {
				if (this.charKmerMatrix[aSeqIdIndex][startPos] == revCompIndex) {
					return (true);
				}
			}
		}

		return (false);
	}

	public double[] getCountsArray(int seqIdIndex, SymbolList[] someWords,
			WeightMatrixTools.BindingStrand strand) {
		double[] counts = new double[someWords.length];

		for (int i = 0; i < someWords.length; i++) {
			if (someWords[i].length() != this.kmerLength) {
				System.out
						.println("Error : SymbolList length is not the same as the words counts length.\n");
				return (null);
			}

			long aWordIndex = getWordIndex(someWords[i], symListToWordMap);
			counts[i] = getCount(seqIdIndex, aWordIndex, strand);
		}
		return (counts);
	}

	public double[] getCountsArray(SymbolList aWord,
			WeightMatrixTools.BindingStrand strand) {
		if (aWord.length() != this.kmerLength) {
			System.out
					.println("Error : SymbolList length is not the same as the words counts length.\n");
			return (null);
		}

		double[] counts = new double[this.numSeqIDs];
		long aWordIndex = getWordIndex(aWord, symListToWordMap);

		for (int aSeqIdIndex = 0; aSeqIdIndex < this.numSeqIDs; aSeqIdIndex++) {
			counts[aSeqIdIndex] = getCount(aSeqIdIndex, aWordIndex, strand);
		}
		return (counts);
	}

	public double[] getCountsArray(SymbolList aWord, double[][] weights,
			WeightMatrixTools.BindingStrand strand) {
		if (aWord.length() != this.kmerLength) {
			System.out
					.println("Error : SymbolList length is not the same as the words counts length.\n");
			return (null);
		}

		double[] counts = new double[this.numSeqIDs];
		long aWordIndex = getWordIndex(aWord, symListToWordMap);

		for (int aSeqIdIndex = 0; aSeqIdIndex < this.numSeqIDs; aSeqIdIndex++) {
			counts[aSeqIdIndex] = getCount(aSeqIdIndex, aWordIndex, weights,
					strand);
		}
		return (counts);
	}

	public double[] getCountsArray(SymbolList aWord, int startPos,
			WeightMatrixTools.BindingStrand strand) {
		if (aWord.length() != this.kmerLength) {
			System.out
					.println("Error : SymbolList length is not the same as the words counts length.\n");
			return (null);
		}

		double[] counts = new double[this.numSeqIDs];
		long aWordIndex = getWordIndex(aWord, symListToWordMap);

		for (int aSeqIdIndex = 0; aSeqIdIndex < this.numSeqIDs; aSeqIdIndex++) {
			counts[aSeqIdIndex] = getCount(aSeqIdIndex, aWordIndex, startPos,
					strand);
		}
		return (counts);
	}

	public double[] getCountsArray(SymbolList aWord, Set<Integer> someProbeIDs,
			WeightMatrixTools.BindingStrand strand) {
		if (aWord.length() != this.kmerLength) {
			System.out
					.println("Error : SymbolList length is not the same as the words counts length.\n");
			return (null);
		}

		double[] counts = new double[someProbeIDs.size()];
		long aWordIndex = getWordIndex(aWord, symListToWordMap);

		int i = 0;
		for (int aSeqID : someProbeIDs) {
			counts[i] = getCount(aSeqID, aWordIndex, strand);
			i++;
		}
		return (counts);
	}

	public double[] getCountsArray(SymbolList aWord, int[] someProbeIDs,
			WeightMatrixTools.BindingStrand strand) {
		if (aWord.length() != this.kmerLength) {
			System.out
					.println("Error : SymbolList length is not the same as the words counts length.\n");
			return (null);
		}

		double[] counts = new double[someProbeIDs.length];
		long aWordIndex = getWordIndex(aWord, symListToWordMap);

		for (int i = 0; i < someProbeIDs.length; i++) {
			counts[i] = getCount(someProbeIDs[i], aWordIndex, strand);
		}
		return (counts);
	}

	public void writeSerializedFile(String aFilePathName) {
		FileTools.writeSerializedFile(this, aFilePathName);
	}

	public static KmerMatrix readSerializedFile(String filePathName) {
		return ((KmerMatrix) FileTools.readSerializedFile(filePathName));
	}

	public void addSymbolLists(String fileName, int kmerLength,
			String alphabetName, boolean windowed)
			throws IllegalAlphabetException, IllegalSymbolException,
			BioException {
		setAlphabet(alphabetName);
		addSymbolLists(fileName, kmerLength, windowed);
	}

	public void addSymbolLists(String fileName, int kmerLength, boolean windowed)
			throws IllegalAlphabetException, IllegalSymbolException,
			BioException {

		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(fileName));
		} catch (Exception ex) {
			ex.printStackTrace(System.err);
		}

		// Replaced sr with a symbollist (or a sequence)
		SequenceIterator seqIter = null;
		Sequence aSequence = null;

		try {
			seqIter = SeqIOTools.readFasta(br, symbolTokenization);
			// aSequence = seqIter.nextSequence();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		addSymbolLists(seqIter, kmerLength, windowed);
	}

	public void addSymbolLists(SequenceIterator iter)
			throws IllegalAlphabetException, IllegalSymbolException,
			BioException {
		addSymbolLists(this.asIterator(iter), 1, false);
	}

	public void addSymbolLists(SequenceIterator iter, int kmerLength,
			boolean windowed) throws IllegalAlphabetException,
			IllegalSymbolException, BioException {
		addSymbolLists(this.asIterator(iter), kmerLength, windowed);
	}

	public void addSymbolLists(Iterator<Sequence> seqIter, int kmerLength,
			boolean windowed) throws IllegalAlphabetException,
			IllegalSymbolException, BioException {

		int seqNameIndex = 0;
		while (seqIter.hasNext()) {
			Sequence source = seqIter.next();
			this.numSeqIDs++;

			addSymbolLists(seqNameIndex, source, kmerLength, windowed);
			seqNameIndex++;
		}
	}

	public void addSymbolLists(int aSeqIdIndex, SymbolList aSymList,
			int kmerLength, boolean windowed) throws IllegalAlphabetException,
			IllegalSymbolException, BioException {

		// for (int startPos = 1; startPos <= aSymList.length() - kmerLength +
		// 1; startPos++) {
		for (int startPos = 1; startPos <= this.numKmerWindows; startPos++) {
			SymbolList aKmerSymList = aSymList.subList(startPos, startPos
					+ kmerLength - 1);
			long aKmerIndex = getWordIndex(aKmerSymList, symListToWordMap);

			if (this.charKmerMatrix != null) {
				this.charKmerMatrix[aSeqIdIndex][startPos - 1] = (char) aKmerIndex;
			} else if (this.intKmerMatrix != null) {
				this.intKmerMatrix[aSeqIdIndex][startPos - 1] = (int) aKmerIndex;
			} else if (this.longKmerMatrix != null) {
				this.longKmerMatrix[aSeqIdIndex][startPos - 1] = aKmerIndex;
			} else if (this.bigIntegerKmerMatrix != null) {

			}

		}
	}

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

    public static final char[] DNA_CHAR_ALPHA = {'A','C','G','T'};

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

	// word.length <= 31
	public static final long getWordIndex(SymbolList aWord) {
		return (getWordIndex(aWord, null));
	}

	// word.length <= 31
	public static final long getWordIndex(SymbolList aWord,
			TObjectCharHashMap aSymListToWordMap) {

		if (aSymListToWordMap != null) {
			return (aSymListToWordMap.get(aWord));
		}

		long wordIndex = 0;
		int len = aWord.length();

		for (int i = 1; i <= len; i++) {
			Symbol symbol = aWord.symbolAt(i);

			// if symbol==A then do nothing!!

			if (symbol == DNATools.c()) {
				// wordIndex += 1 * Math.pow(4,(len-i));
				wordIndex += multBy4ToPower(1L, len - i);
				// printf("maxmer[%d]=%c; wordIndex=%d ; maxWordSize=%d | ", i,
				// maxmer[i], wordIndex, maxWordSize);
			} else if (symbol == DNATools.g()) {
				// wordIndex += 2 * Math.pow(4,(len-i));
				wordIndex += multBy4ToPower(2L, len - i);
				// printf("maxmer[%d]=%c; wordIndex=%d ; maxWordSize=%d | ", i,
				// maxmer[i], wordIndex, maxWordSize);
			} else if (symbol == DNATools.t()) {
				// wordIndex += 3 * Math.pow(4,(len-i));
				wordIndex += multBy4ToPower(3L, len - i);
				// printf("maxmer[%d]=%c; wordIndex=%d ; maxWordSize=%d | ", i,
				// maxmer[i], wordIndex, maxWordSize);
			}
		}
		return (wordIndex);
	}

	// word.length <= 31
	// 1-based indexing
	// start_incl, end_incl
	public static final long getWordIndex(SymbolList aWord, int start_incl,
			int end_incl) {

		long wordIndex = 0;
		int len = end_incl - start_incl + 1;

		for (int i = 0; i < len; i++) {
			Symbol symbol = aWord.symbolAt(start_incl + i);

			// if symbol==A then do nothing!!

			if (symbol == DNATools.c()) {
				// wordIndex += 1 * Math.pow(4,(len-i-1));
				wordIndex += multBy4ToPower(1L, len - i - 1);
				// printf("maxmer[%d]=%c; wordIndex=%d ; maxWordSize=%d | ", i,
				// maxmer[i], wordIndex, maxWordSize);
			} else if (symbol == DNATools.g()) {
				// wordIndex += 2 * Math.pow(4,(len-i-1));
				wordIndex += multBy4ToPower(2L, len - i - 1);
				// printf("maxmer[%d]=%c; wordIndex=%d ; maxWordSize=%d | ", i,
				// maxmer[i], wordIndex, maxWordSize);
			} else if (symbol == DNATools.t()) {
				// wordIndex += 3 * Math.pow(4,(len-i-1));
				wordIndex += multBy4ToPower(3L, len - i - 1);
				// printf("maxmer[%d]=%c; wordIndex=%d ; maxWordSize=%d | ", i,
				// maxmer[i], wordIndex, maxWordSize);
			}
		}
		return (wordIndex);
	}

	// word.length <= 31
	public static final long getWordIndex(String word) {
		long wordIndex = 0;
		int len = word.length();

		for (int i = 0; i < len; i++) {

            // a lookup table is marginally faster than a switch (as long as the lookup table is in memory)
            wordIndex += multBy4ToPower((long) (charToLong[word.charAt(i)]), len - i - 1);

            // switch (word.charAt(i)) {
			// case 'A':
			// case 'a':
			// 	// wordIndex += 0;
			// 	// printf("maxmer[%d]=%c; wordIndex=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], wordIndex, maxWordSize);
			// 	break;
			// case 'C':
			// case 'c':
			// 	// wordIndex += 1 * Math.pow(4,(len-i-1));
			// 	wordIndex += multBy4ToPower(1L, len - i - 1);
			// 	// printf("maxmer[%d]=%c; wordIndex=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], wordIndex, maxWordSize);
			// 	break;
			// case 'G':
			// case 'g':
			// 	// wordIndex += 2 * Math.pow(4,(len-i-1));
			// 	wordIndex += multBy4ToPower(2L, len - i - 1);
			// 	// printf("maxmer[%d]=%c; wordIndex=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], wordIndex, maxWordSize);
			// 	break;
			// case 'T':
			// case 't':
			// 	// wordIndex += 3 * Math.pow(4,(len-i-1));
			// 	wordIndex += multBy4ToPower(3L, len - i - 1);
			// 	// printf("maxmer[%d]=%c; wordIndex=%d ; maxWordSize=%d | ", i,
			// 	// maxmer[i], wordIndex, maxWordSize);
			// 	break;
			// }
		}
		return (wordIndex);
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

	// word.length <= 7
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

	// word.length <= 15
	// 0-based indexing, (start, end]
	public static final int getSubWord(int word, int len, int start, int end) {
		int reverseIndex = 0;

		for (int i = 0; i < len; i++) {
			int mod = word % 4;
			int pos = len - i - 1;

			if (pos < end) {
				reverseIndex += multBy4ToPower(mod, i - (len - end));
			}

			// fprintf(fout, "word=%d, mod=%d, reverseIndex=%d, i=%d\n", word,
			// mod, reverseIndex, i);
			word = (word - mod) / 4;

			if (pos <= start) {
				break;
			}
		}

		// fprintf(fout, "\n\n");
		return (reverseIndex);
	}

	// word.length <= 31
	// 0-based indexing, (start, end]
	public static final long getSubWord(long word, int len, int start, int end) {
		long reverseIndex = 0;

		for (int i = 0; i < len; i++) {
			long mod = word % 4;
			int pos = len - i - 1;

			if (pos < end) {
				reverseIndex += multBy4ToPower(mod, i - (len - end));
			}

			// fprintf(fout, "word=%d, mod=%d, reverseIndex=%d, i=%d\n", word,
			// mod, reverseIndex, i);
			word = (word - mod) / 4;

			if (pos <= start) {
				break;
			}
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
	public SymbolList getSymbolList(long wordIndex, int len) {
		return (getSymbolList(wordIndex, len, null, alphabet));
	}

	// word.length <= 31
	public static final SymbolList getSymbolList(long wordIndex, int len,
			SymbolList[] aWordToSymListMatrix, Alphabet anAlphabet) {
		if (aWordToSymListMatrix != null) {
			return (aWordToSymListMatrix[(int) wordIndex]);
		} else {
			return (new SimpleSymbolList(getSymbolArray(wordIndex, len), len,
					anAlphabet));
		}
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
            wordCharArray[len - i - 1] = DNA_CHAR_ALPHA[mod];

			temp = (temp - mod) / 4L;

			// if ((temp % 4L) != 0) {
			// System.err.println("\nERROR: temp="+temp+", temp%4L="+temp % 4L);
			// }
		}
		return (wordCharArray);
	}

	// word.length <= 31
	public static final String getWordString(long wordIndex, int len) {
		return (getWordString(wordIndex, len, null));
	}

	// word.length <= 31
	public static final String getWordString(long wordIndex, int len,
			SymbolList[] aWordToSymListMatrix) {
		if (aWordToSymListMatrix != null) {
			return (aWordToSymListMatrix[(int) wordIndex].seqString());
		}

		return (new String(getWordCharArray(wordIndex, len)));
	}

	/**
	 * Makes a <code>SequenceIterator</code> look like an
	 * <code>Iterator {@code <Sequence>}</code>
	 * 
	 * @param iter
	 *            The <CODE>SequenceIterator</CODE>
	 * @return An <CODE>Iterator</CODE> that returns only <CODE>Sequence</CODE>
	 *         objects. <B>You cannot call <code>remove()</code> on this
	 *         iterator!</B>
	 */
	public Iterator<Sequence> asIterator(SequenceIterator iter) {
		final SequenceIterator it = iter;

		return new Iterator<Sequence>() {
			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public Sequence next() {
				try {
					return it.nextSequence();
				} catch (BioException e) {
					NoSuchElementException ex = new NoSuchElementException();
					ex.initCause(e);
					throw ex;
				}
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}

	// max length = 8
	SymbolList[] makeWordToSymListMatrix(int kmerLength) {
		// int maxWord = (int)Math.pow(4, kmerLength);
		int maxWord = multBy4ToPower(1, kmerLength);

		SymbolList[] aWordToSymListMatrix = new SymbolList[maxWord];

		// be sure to use "int word" and not "char word" or infinite loop for
		// 8mers
		for (int word = 0; word < maxWord; word++) {
			SymbolList symList = getSymbolList(word, kmerLength);
			aWordToSymListMatrix[word] = symList;
		}

		return (aWordToSymListMatrix);
	}

	// max length = 8
	TObjectCharHashMap makeSymListToWordMap(int kmerLength) {
		// int maxWord = (int)Math.pow(4, kmerLength);
		int maxWord = multBy4ToPower(1, kmerLength);

		TObjectCharHashMap aSymListToWordMap = new TObjectCharHashMap();

		// be sure to use "int word" and not "char word" or infinite loop for
		// 8mers
		for (int word = 0; word < maxWord; word++) {
			SymbolList symList = getSymbolList(word, kmerLength);
			aSymListToWordMap.put(symList, (char) word);
		}

		return (aSymListToWordMap);
	}

	// max length = 8
	char[] makeRevCompMatrix(int kmerLength) {
		// int maxWord = (int)Math.pow(4, kmerLength);
		int maxWord = multBy4ToPower(1, kmerLength);

		char[] aRevCompMatrix = new char[maxWord];

		// be sure to use "int word" and not "char word" or infinite loop for
		// 8mers
		for (int word = 0; word < maxWord; word++) {
			char revComp = getRevComp((char) word, kmerLength);
			aRevCompMatrix[word] = revComp;

			// if((((int)word) % 1000) == 0 ) {
			// System.out.println("\n\tAssociated "+((int)word)+" out of "+maxWord+" reverse complements so far...");
			// }

		}

		return (aRevCompMatrix);
	}

	// double getKmerAffinity(char wordIndex, int wordLen, double[] intensities,
	// double[][] weights) {
	// char revCompIndex = revCompMatrix[wordIndex];
	// ArrayList<Double> affinities = new ArrayList<Double>();

	// for (int aSeqIdIndex=0; aSeqIdIndex < this.numSeqIDs; aSeqIdIndex++) {
	// double weightedCount = getCount(aSeqIdIndex, wordIndex, weights,
	// WeightMatrixTools.BindingStrand.BOTH);
	// if (weightedCount > 0) {
	// affinities.add(intensities[aSeqIdIndex] / weightedCount);
	// }
	// }

	// //System.out.println("\tNumber of probes that contain this "+kmerLength+"-mer or its revComp is "+affinities.size());

	// // total drop fraction = .20
	// //double affinity =
	// MathTools.trimmedMean(ArrayTools.toDoubleArray(affinities), .20);
	// double affinity =
	// MathTools.trimmedMean(ArrayTools.toDoubleArray(affinities), .30);

	// return(affinity);
	// }

	double getKmerAffinity(char wordIndex, int wordLen, double[] intensities,
			double[][] weights, double fractionDrop) {
		char revCompIndex = revCompMatrix[wordIndex];
		ArrayList<Double> affinities = new ArrayList<Double>();

		for (int aSeqIdIndex = 0; aSeqIdIndex < this.numSeqIDs; aSeqIdIndex++) {
			double weightedCount = getCount(aSeqIdIndex, wordIndex, weights,
					WeightMatrixTools.BindingStrand.POS);
			if (weightedCount > 0) {
				affinities.add(intensities[aSeqIdIndex] / weightedCount);
			}
		}

		// System.out.println("\tNumber of probes that contain this "+kmerLength+"-mer or its revComp is "+affinities.size());

		// total drop fraction = .30
		// double affinity = MathTools.trimmedMean(affinities, 0, count,
		// fractionDrop);
		// double affinity = MathTools.trimmedMean(affinities, 0, count,
		// fractionDrop);
		// double affinity = MathTools.trimmedMean(affinities, 0, count, .30);

		double affinity = MathTools.trimmedMean(
				ArrayTools.toDoubleArray(affinities), fractionDrop);
		// double affinity =
		// MathTools.trimmedMean(ArrayTools.toDoubleArray(affinities), .20);
		// double affinity =
		// MathTools.trimmedMean(ArrayTools.toDoubleArray(affinities), .30);

		return (affinity);
	}

	double getKmerAffinity(
			int index,
			char[] wordsArray, // partial array of kmers
			int wordLen,
			double[] intensities,
			double[][] weights,
			char[][] kmerToProbesMatrix, // partial array of kmer keys, can
											// include seqs with revComps or not
			double[] oldKmerToAffinityMatrix, // full array for all kmers
			double fractionDrop, boolean nonNegativeFlag, double percentStep,
			double initValue, int probesPerKmerStandard,
			double probesPerKmerMultFactor, boolean includeRevComps) {
		char word = wordsArray[index];
		char revCompWord = revCompMatrix[word];
		char[] kmerProbeIDs = kmerToProbesMatrix[index];

		ArrayList<Double> affinitiesArrayList = new ArrayList<Double>();
		MutableDouble strandAffinity = new MutableDouble();

		if (nonNegativeFlag && (initValue < 0)) {
			initValue = 0;
		}

		// if (kmerProbeIDs.length < 9) {
		// if (oldKmerToAffinityMatrix != null) {
		// return(oldKmerToAffinityMatrix[revCompWord]);
		// //return(0);
		// }
		// else {
		// return(initValue);
		// //return(0); // needed by sparse matrix to get rid of these kmers
		// }
		// }

		// Some poly-G K-mers don't exist on any or only a few probes!!!!
		// set them to rev-comp affinity as a guess
		if (kmerProbeIDs.length == 0) {
			return (0);
		}

		// if (kmerProbeIDs.length == 0) {

		// //if (kmerProbeIDs.length < 10) {
		// //System.out.println("PROBLEM: affinitiesArrayList.size() == 0, kmer="+getWordString(word,
		// wordLen));
		// if (oldKmerToAffinityMatrix != null) {
		// double affinity = oldKmerToAffinityMatrix[revCompWord];

		// // scale the affinity by the weights ratio
		// double posMean = MathTools.mean(weights[0], 2, 7);
		// double negMean = MathTools.mean(weights[1], 2, 7);
		// // double posMean = MathTools.mean(weights[0]);
		// // double negMean = MathTools.mean(weights[1]);
		// if ((negMean != 0) && (posMean != 0)) {
		// affinity *= posMean/negMean;
		// }

		// return(affinity);
		// }
		// else {
		// return(initValue);
		// }
		// }

		for (int i = 0; i < kmerProbeIDs.length; i++) {
			int aSeqIdIndex = kmerProbeIDs[i];

			// ProbeIDs also include revcomps
			// if not include revComps then check to make sure it contains Kmer
			// if (!includeRevComps) {

			// int count = getCount(
			// aSeqIdIndex,
			// word,
			// WeightMatrixTools.BindingStrand.POS);

			// if (count == 0) {
			// continue;
			// }
			// }

			double weightedCount = getCount(aSeqIdIndex, word, weights,
					// WeightMatrixTools.BindingStrand.BOTH,
					oldKmerToAffinityMatrix, strandAffinity, initValue,
					includeRevComps);

			// affinities[i] = intensities[aSeqIdIndex] / weightedCount;

			double oldAffinity = initValue;
			if (oldKmerToAffinityMatrix != null) {
				oldAffinity = oldKmerToAffinityMatrix[word]; // oldKmerToAffinityMatrix
																// is a full
																// Kmers array
			}

			double newAffinity = (intensities[aSeqIdIndex] - strandAffinity
					.doubleValue()) / weightedCount;
			// double newAffinity = intensities[aSeqIdIndex] / weightedCount;

			// enforce non-negativity
			if ((nonNegativeFlag) && (newAffinity < 0)) {
				newAffinity = 0;
			} else if (weightedCount <= 0) {
				newAffinity = 0;
			}

			if (true) {
				// if (oldKmerToAffinityMatrix != null) {
				double affinityDelta = newAffinity - oldAffinity;
				affinitiesArrayList.add(oldAffinity
						+ (affinityDelta * percentStep));
			} else { // NEED THIS TO BE ABLE TO CALL METHOD TO GET TRIMMEDMEANS
				affinitiesArrayList.add(newAffinity);
			}

		}

		// add step to 0 !!!
		// includes for revComps and not revComps
		for (int i = 0; i < probesPerKmerStandard * probesPerKmerMultFactor; i++) {
			double oldAffinity = initValue;
			if (oldKmerToAffinityMatrix != null) {
				oldAffinity = oldKmerToAffinityMatrix[word]; // oldKmerToAffinityMatrix
																// is a full
																// Kmers array
			}

			double affinityDelta = 0 - oldAffinity;
			affinitiesArrayList
					.add(oldAffinity + (affinityDelta * percentStep));
		}

		// DON'T USE REVCOMPS
		// double[] affinities = ArrayTools.toDoubleArray(affinitiesArrayList);
		// double affinity = MathTools.trimmedMean(affinities, fractionDrop);

		// WEIGHTED AVERAGE
		// double affinity = 0;
		// double[] affinities = ArrayTools.toDoubleArray(affinitiesArrayList);
		// double kmerAffinity = MathTools.trimmedMean(affinities,
		// fractionDrop);

		// if (kmerProbeIDs.length >= kmerToProbesMatrix[revCompWord].length) {
		// affinity = kmerAffinity;
		// }
		// else {
		// double rcKmerAffinity = 0;
		// if (oldKmerToAffinityMatrix != null) {
		// rcKmerAffinity = oldKmerToAffinityMatrix[revCompWord];
		// }
		// else {
		// rcKmerAffinity = initValue;
		// }

		// double fraction = (2 * kmerProbeIDs.length) / (kmerProbeIDs.length +
		// kmerToProbesMatrix[revCompWord].length);
		// affinity = (fraction * kmerAffinity) + ((1 - fraction) *
		// rcKmerAffinity);
		// }

		// if (affinitiesArrayList.size() < 9) {
		// return(0);
		// }

		// USE AFFINITY(REVCOMP) WHEN PROBECOUNTS < THRESH
		double affinity = 0;
		// Some poly-G K-mers don't exist on any of the probes!!!!
		// set them to rev-comp affinity as a guess
		if (false) {
			// if (affinitiesArrayList.size() == 0) {
			// if ((affinitiesArrayList.size() < 10)
			// && (kmerProbeIDs.length <
			// kmerToProbesMatrix[revCompWord].length)){

			// System.out.println("PROBLEM: affinitiesArrayList.size() == 0, kmer="+getWordString(word,
			// wordLen));
			if (oldKmerToAffinityMatrix != null) {
				affinity = oldKmerToAffinityMatrix[revCompWord];

				// scale the affinity by the weights ratio
				double posMean = MathTools.mean(weights[0], 2, 7);
				double negMean = MathTools.mean(weights[1], 2, 7);
				// double posMean = MathTools.mean(weights[0]);
				// double negMean = MathTools.mean(weights[1]);
				if ((negMean != 0) && (posMean != 0)) {
					affinity *= posMean / negMean;
				}
			} else {
				affinity = initValue;
			}
		} else {
			double[] affinities = ArrayTools.toDoubleArray(affinitiesArrayList);
			affinity = MathTools.trimmedMean(affinities, fractionDrop);
		}

		if (Double.isNaN(affinity)) {
			// System.out.println("\nError: affinity is NaN!!!!");
			// affinity = initValue;
			affinity = 0;
		}

		// System.out.println("\tNumber of probes that contain this "+kmerLength+"-mer or its revComp is "+affinities.size());

		// total drop fraction = .30

		// double affinity = MathTools.trimmedMean(affinities, 0, count, .20);
		// double affinity = MathTools.trimmedMean(affinities, 0, count, .30);

		// double affinity =
		// MathTools.trimmedMean(ArrayTools.toDoubleArray(affinities), .20);
		// double affinity =
		// MathTools.trimmedMean(ArrayTools.toDoubleArray(affinities), .30);

		return (affinity);
	}

	double getAvgProbesPerKmer(char[][] aWordToProbesMatrix) {
		double sum = 0;
		int nullCount = 0;

		for (int i = 0; i < aWordToProbesMatrix.length; i++) {
			if (aWordToProbesMatrix[i] != null) {
				sum += aWordToProbesMatrix[i].length;
			} else {
				nullCount++;
			}
		}

		if (nullCount > 0) {
			System.out
					.println("\nWarning!!!!: nullCount in getAvgProbesPerKmer is "
							+ nullCount + "!!");
			System.out.println("\naWordToProbesMatrix.length = "
					+ aWordToProbesMatrix.length + ".\n");
		}

		return (sum / aWordToProbesMatrix.length);
	}

	char[][] makeKmerToProbesMatrix(int kmerLength, char[] wordsArray,
			boolean includeRevComps) {
		int maxWord = multBy4ToPower(1, kmerLength);
		TCharArrayList[] kmerToProbesArrayList = new TCharArrayList[maxWord];

		for (int i = 0; i < kmerToProbesArrayList.length; i++) {
			kmerToProbesArrayList[i] = new TCharArrayList();
		}

		// go through all the probe sequences
		for (int aSeqIdIndex = 0; aSeqIdIndex < this.numSeqIDs; aSeqIdIndex++) {

			for (int startPos = 0; startPos < this.charKmerMatrix[aSeqIdIndex].length; startPos++) {
				char aWordIndex = this.charKmerMatrix[aSeqIdIndex][startPos];
				kmerToProbesArrayList[aWordIndex].add((char) aSeqIdIndex);

				// Also add revComp if not self-reverse-complement
				if (includeRevComps) {
					// if (true) {
					char revCompIndex = revCompMatrix[aWordIndex];
					if (revCompIndex != aWordIndex) {
						kmerToProbesArrayList[revCompIndex]
								.add((char) aSeqIdIndex);
					}
				}

			}
		}

		// go through all the kmers
		int missingKmersCount = 0;
		char[][] kmerToProbesMatrix = new char[wordsArray.length][];
		Table probeCounts = new Table();

		for (int index = 0; index < wordsArray.length; index++) {
			// kmerToProbesMatrix[wordsArray[index]] = new
			// char[kmerToProbesArrayList[wordsArray[index]].size()];
			kmerToProbesMatrix[index] = kmerToProbesArrayList[wordsArray[index]]
					.toArray();

			if (kmerToProbesMatrix[index].length == 0) {
				missingKmersCount++;
			}

			probeCounts.add(Arrays.asList(
					new String(getWordString(wordsArray[index], kmerLength)),
					new Integer(kmerToProbesMatrix[index].length)));
		}

		probeCounts.sort(1);
		// FileTools.write(probeCounts.toString(0, 1000, "\t"), kmerLength
		// 		+ "mers.probeCounts.posStrand.bottom1000.txt", false);

		if (missingKmersCount > 0) {
			System.out
					.println("\n\nWarning!!!!: "
							+ missingKmersCount
							+ " "
							+ kmerLength
							+ "-mers are not found anywhere in any of the probe sequences on the positive strand!!\n");
			int lessThanCount = probeCounts.lessThanCount(1, 1);
			List<Object> missingKmers = probeCounts.getColumn(0,
					0, lessThanCount);
			// System.out.println("Missing "+kmerLength+"-mers are :"+StringTools.toString(missingKmers));
			System.out.println("The missing " + kmerLength + "-mers are "
					+ missingKmers.toString());
		}

		System.out.println("\nNumber of " + kmerLength
				+ "-mers found on <= 1 probe (positive strand) is "
				+ probeCounts.lessThanCount(1, 2) + ".");
		System.out.println("Number of " + kmerLength
				+ "-mers found on <= 2 probes (positive strand) is "
				+ probeCounts.lessThanCount(1, 3) + ".");
		System.out.println("Number of " + kmerLength
				+ "-mers found on <= 4 probes (positive strand) is "
				+ probeCounts.lessThanCount(1, 5) + ".");
		System.out.println("Number of " + kmerLength
				+ "-mers found on <= 8 probes (positive strand) is "
				+ probeCounts.lessThanCount(1, 9) + ".");
		System.out.println("Number of " + kmerLength
				+ "-mers found on <= 12 probes (positive strand) is "
				+ probeCounts.lessThanCount(1, 13) + ".\n");

		return (kmerToProbesMatrix);

	}

	// TCharDoubleHashMap makeKmerToAffinityMap(int kmerLength, double[]
	// intensities, double[][] positionalBias) {
	// TCharDoubleHashMap kmerToAffinityMap = new TCharDoubleHashMap();

	// int maxWord = (int)Math.pow(4, kmerLength);

	// // be sure to use "int word" and not "char word" or infinite loop for
	// 8mers
	// for (int word=0; word < maxWord; word++) {
	// double affinity = getKmerAffinity((char)word, kmerLength, intensities,
	// positionalBias);
	// kmerToAffinityMap.put((char)word, affinity);
	// }

	// return(kmerToAffinityMap);
	// }

	double[] makeKmerToTrimmedMeanMatrix(int kmerLength, double[] intensities,
			double[][] positionalBias, char[] wordsArray, // partial Kmers array
			char[][] kmerToProbesMatrix, // partial Kmers array
			double fractionDrop, boolean includeRevComps) {
		double[] kmerToTrimmedMeanMatrix = new double[wordsArray.length];

		for (int index = 0; index < wordsArray.length; index++) {

			double trimmedMean = getKmerAffinity(index,
					wordsArray, // partial Kmers array
					kmerLength, intensities, positionalBias,
					kmerToProbesMatrix, // partial Kmers keys
					null, // oldKmerToAffinityMatrix
					fractionDrop, false, // nonNegativeFlag
					1.0, // percentStep
					0.0, // initValue
					0, // probesPerKmerStandard
					0, // probesPerKmerMultFactor
					includeRevComps);

			kmerToTrimmedMeanMatrix[index] = trimmedMean;

			if (((index) % 5000) == 0) {
				System.out.println("\n\tCreated affinities for " + index
						+ " out of " + wordsArray.length + " " + kmerLength
						+ "-mers so far...");
			}
		}

		return (kmerToTrimmedMeanMatrix);
	}

	double[] makeKmerToAffinityMatrix(int kmerLength, double[] intensities,
			double[][] positionalBias, double fractionDrop) {
		// int maxWord = (int)Math.pow(4, kmerLength);
		int maxWord = multBy4ToPower(1, kmerLength);

		double[] kmerToAffinityMatrix = new double[maxWord];

		// be sure to use "int word" and not "char word" or infinite loop for
		// 8mers
		for (int word = 0; word < maxWord; word++) {
			double affinity = getKmerAffinity((char) word, kmerLength,
					intensities, positionalBias, fractionDrop);
			kmerToAffinityMatrix[word] = affinity;

			if (((word) % 5000) == 0) {
				System.out.println("\n\tCreated affinities for " + (word)
						+ " out of " + maxWord + " " + kmerLength
						+ "-mers so far...");
			}
		}

		return (kmerToAffinityMatrix);
	}

	// kmerLength <= 8
	// public TCharHashSet makeRevCompEquivKmers(int kmerLength) {
	public char[] makeRevCompEquivKmers(int kmerLength) {
		TCharHashSet words = new TCharHashSet();
		int maxWord = multBy4ToPower(1, kmerLength);

		// be sure to use "int word" and not "char word" or infinite loop for
		// 8mers
		for (int word = 0; word < maxWord; word++) {
			char revComp = revCompMatrix[word];
			if (!words.contains(revComp)) {
				words.add((char) word);
			}
		}
		// return(words);
		return (words.toArray());
	}

	// kmerLength <= 8
	public static char[] makeAllKmers(int kmerLength) {
		int maxWord = multBy4ToPower(1, kmerLength);
		char[] wordsArray = new char[maxWord];

		// be sure to use "int word" and not "char word" or infinite loop for
		// 8mers
		for (int word = 0; word < maxWord; word++) {
			wordsArray[word] = (char) word;
		}
		return (wordsArray);
	}

	double[][] makeSparseKmerMatrix(int kmerLength,
			char[] wordsArray, // partial array of kmers
			char[][] kmerToProbesMatrix, // partial array of kmer keys
			double[][] positionalBiases, double[] intensities,
			boolean includeRevComps, boolean intercept) {
		Table sparseMatrixTable = new Table();
		// TIntHashSet probeIDs = new TIntHashSet();
		// TIntHashSet words = new TIntHashSet();

		// loop through all the kmers
		for (int kmerIndex = 0; kmerIndex < wordsArray.length; kmerIndex++) {

			int kmer = wordsArray[kmerIndex];
			char[] kmerToProbesArray = kmerToProbesMatrix[kmerIndex];

			// for this kmer loop through all the probes
			for (int probeIndex = 0; probeIndex < kmerToProbesArray.length; probeIndex++) {
				int probe = kmerToProbesArray[probeIndex];

				// traverse the probe and get the kmer count for this probe
				// int count = getCount(probe, kmer,
				// WeightMatrixTools.BindingStrand.POS);
				// if (count == 0) {
				// //System.out.println("\nError: K-mer not found on probe!");
				// continue;
				// }

				double weightedCount;
				if (!includeRevComps) {
					weightedCount = getCount(probe, kmer, positionalBiases,
							WeightMatrixTools.BindingStrand.POS);
				} else {
					weightedCount = getCount(probe, kmer, positionalBiases,
							WeightMatrixTools.BindingStrand.BOTH);
				}

				// update the sparseKmerMatrix with this weightedCount
				sparseMatrixTable.add(Arrays.asList(new Double(probe + 1),
						new Double(kmerIndex + 1), new Double(weightedCount)));

				// probeIDs.add(probe+1);
				// words.add(kmerIndex+1);
			}
		}

		// if add an intercept, then add them as last column
		if (intercept) {
			// go through all the probe sequences
			// for (int aSeqIdIndex=0; aSeqIdIndex < newIntensities.size();
			// aSeqIdIndex++) {
			for (int aSeqIdIndex = 0; aSeqIdIndex < this.numSeqIDs; aSeqIdIndex++) {
				// add a column of all ones as the last column

				sparseMatrixTable.add(Arrays.asList(
						new Double(aSeqIdIndex + 1), // row
						new Double(wordsArray.length + 1), // intercept column
						new Double(1.0))); // value

			}
		}

		// create the array of intensities that contain these K-mers
		double[][] sparseKmerMatrix = new double[4][];

		sparseKmerMatrix[0] = ArrayTools
				.toDoubleArray((ArrayList) sparseMatrixTable.getColumn(0));
		sparseKmerMatrix[1] = ArrayTools
				.toDoubleArray((ArrayList) sparseMatrixTable.getColumn(1));
		sparseKmerMatrix[2] = ArrayTools
				.toDoubleArray((ArrayList) sparseMatrixTable.getColumn(2));
		sparseKmerMatrix[3] = intensities;

		// System.out.println("\n\tmin(rows)="+MathTools.min(sparseKmerMatrix[0])+"; max(rows)="+MathTools.max(sparseKmerMatrix[0]));
		// System.out.println("\n\tmin(kmers)="+MathTools.min(sparseKmerMatrix[1])+"; max(kmers)="+MathTools.max(sparseKmerMatrix[1]));
		// System.out.println("\n\tmin(values)="+MathTools.min(sparseKmerMatrix[2])+"; max(values)="+MathTools.max(sparseKmerMatrix[2]));
		// System.out.println("\n\tmin(intensities)="+MathTools.min(sparseKmerMatrix[3])+"; max(intensities)="+MathTools.max(sparseKmerMatrix[3]));

		return (sparseKmerMatrix);
	}

	// double[][] makeSparseKmerMatrix(
	// int kmerLength,
	// char[] wordsArray,
	// char[][] kmerToProbesMatrix,
	// double[][] positionalBiases,
	// double[] intensities,
	// boolean intercept)
	// {
	// Table sparseMatrixTable = new Table();
	// TIntHashSet probeIDs = new TIntHashSet();
	// TIntHashSet words = new TIntHashSet();

	// // loop through all the kmers
	// for (int kmerIndex=0; kmerIndex < wordsArray.length; kmerIndex++) {

	// int kmer = wordsArray[kmerIndex];
	// char[] kmerToProbesArray = kmerToProbesMatrix[kmer];

	// // for this kmer loop through all the probes
	// for (int probeIndex=0; probeIndex < kmerToProbesArray.length;
	// probeIndex++) {
	// int probe = kmerToProbesArray[probeIndex];

	// // traverse the probe and get the kmer count for this probe
	// int count = getCount(probe, kmer, WeightMatrixTools.BindingStrand.POS);
	// if (count == 0) {
	// //System.out.println("\nError: K-mer not found on probe!");
	// continue;
	// }

	// double weightedCount = getCount(probe, kmer, positionalBiases,
	// WeightMatrixTools.BindingStrand.POS);

	// // update the sparseKmerMatrix with this weightedCount
	// sparseMatrixTable.add(
	// Arrays.asList(
	// new Double(probe+1),
	// new Double(kmerIndex+1),
	// new Double(weightedCount)));

	// probeIDs.add(probe+1);
	// words.add(kmerIndex+1);
	// }
	// }

	// // sort the sparseMatrixTable by ascending probeID
	// sparseMatrixTable.sort(0);

	// // compress the table by getting rid of empty rows
	// int currentNewRow = 1;
	// int currentOldRow = ((Double)sparseMatrixTable.getElement(0,
	// 0)).intValue();

	// int skipped=0;
	// for (int i=0; i < sparseMatrixTable.rows(); i++) {
	// if (((Double)sparseMatrixTable.getElement(i, 0)).intValue() !=
	// currentOldRow) {

	// skipped += ((Double)sparseMatrixTable.getElement(i, 0)).intValue() -
	// currentOldRow -1;

	// // update CNR
	// currentNewRow++;

	// // update COR
	// currentOldRow = ((Double)sparseMatrixTable.getElement(i, 0)).intValue();
	// }

	// // set the row to the currentNewRow
	// sparseMatrixTable.setElement(i, 0, new Double(currentNewRow));
	// }
	// System.out.println("\nRemoved "+skipped+" rows when compressing the sparse matrix.");

	// // compress the intensities by getting rid of empty rows
	// ArrayList<Double> newIntensities = new ArrayList<Double>();
	// ArrayList<Double> skippedIntensities = new ArrayList<Double>();
	// skipped = 0;
	// double skippedMean = 0;
	// for (int i=0; i < intensities.length; i++) {
	// if (probeIDs.contains(i+1)) {
	// newIntensities.add(intensities[i]);
	// }
	// else {
	// skipped++;
	// skippedIntensities.add(intensities[i]);
	// }
	// }
	// if (skipped > 0) {
	// skippedMean = MathTools.mean(skippedIntensities);
	// }
	// System.out.println("\nRemoved "+skipped+" rows in new residuals Array. Mean(skippedIntensities) = "+skippedMean);

	// //
	// System.out.println("\n\tnewIntensities.size()="+newIntensities.size()+"; probeIDs.size()="+probeIDs.size());
	// //
	// System.out.println("\n\twords.size()="+words.size()+"; wordsArray.length="+wordsArray.length);

	// // if add an intercept, then add them as last column
	// if (intercept) {
	// // go through all the probe sequences
	// for (int aSeqIdIndex=0; aSeqIdIndex < newIntensities.size();
	// aSeqIdIndex++) {
	// //for (int aSeqIdIndex=0; aSeqIdIndex < this.numSeqIDs; aSeqIdIndex++) {
	// //add a column of all ones as the last column

	// sparseMatrixTable.add(
	// Arrays.asList(
	// new Double(aSeqIdIndex+1), // row
	// new Double(words.size()+1), // intercept column
	// new Double(1.0))); // value

	// }
	// }

	// if (skipped > 0) {
	// // add one row that has only the intercept
	// sparseMatrixTable.add(
	// Arrays.asList(
	// new Double(newIntensities.size()+1), // new row
	// new Double(words.size()+1), // intercept column
	// new Double(1.0))); // value
	// newIntensities.add(skippedMean);
	// }

	// // create the array of intensities that contain these K-mers
	// double[][] sparseKmerMatrix = new double[4][];

	// sparseKmerMatrix[0] =
	// ArrayTools.toDoubleArray((ArrayList)sparseMatrixTable.getColumn(0));
	// sparseKmerMatrix[1] =
	// ArrayTools.toDoubleArray((ArrayList)sparseMatrixTable.getColumn(1));
	// sparseKmerMatrix[2] =
	// ArrayTools.toDoubleArray((ArrayList)sparseMatrixTable.getColumn(2));
	// sparseKmerMatrix[3] = ArrayTools.toDoubleArray(newIntensities);

	// //
	// System.out.println("\n\tmin(rows)="+MathTools.min(sparseKmerMatrix[0])+"; max(rows)="+MathTools.max(sparseKmerMatrix[0]));
	// //
	// System.out.println("\n\tmin(kmers)="+MathTools.min(sparseKmerMatrix[1])+"; max(kmers)="+MathTools.max(sparseKmerMatrix[1]));
	// //
	// System.out.println("\n\tmin(values)="+MathTools.min(sparseKmerMatrix[2])+"; max(values)="+MathTools.max(sparseKmerMatrix[2]));
	// //
	// System.out.println("\n\tmin(intensities)="+MathTools.min(sparseKmerMatrix[3])+"; max(intensities)="+MathTools.max(sparseKmerMatrix[3]));

	// return(sparseKmerMatrix);
	// }

	// ///////////////////////////////////////////////////////////////////////////////////////////////////////
	// Now performs in-place replacement of affinity values and returns back the
	// oldKmerToAffinity Matrix
	// if oldKmerToAffinityMatrix != null
	// ///////////////////////////////////////////////////////////////////////////////////////////////////////
	double[] makeKmerToAffinityMatrix(
			int kmerLength,
			char[] kmersArray, // partial array of kmers
			double[] intensities,
			double[][] positionalBias,
			char[][] kmerToProbesMatrix, // partial array of kmer keys
			double[] oldKmerToAffinityMatrix, // full array for all kmers
			double fractionDrop, boolean nonNegativeFlag, double percentStep,
			double initValue, int probesPerKmerStandard,
			double probesPerKmerMultFactor, boolean includeRevComps) {

		double[] newKmerToAffinityMatrix = null;
		if (oldKmerToAffinityMatrix == null) {
			int maxWord = multBy4ToPower(1, kmerLength);
			newKmerToAffinityMatrix = new double[maxWord];
			Arrays.fill(newKmerToAffinityMatrix, initValue);
		} else {
			newKmerToAffinityMatrix = oldKmerToAffinityMatrix;
		}

		for (int index = 0; index < kmersArray.length; index++) {

			char word = kmersArray[index];

			double affinity = getKmerAffinity(index, kmersArray, kmerLength,
					intensities, positionalBias, kmerToProbesMatrix,
					oldKmerToAffinityMatrix, fractionDrop, nonNegativeFlag,
					percentStep, initValue, probesPerKmerStandard,
					probesPerKmerMultFactor, includeRevComps);

			newKmerToAffinityMatrix[word] = affinity;

			if (includeRevComps) {
				char revComp = revCompMatrix[word];
				newKmerToAffinityMatrix[revComp] = affinity;
			}

			// if((((int)word) % 5000) == 0 ) {
			// System.out.println("\n\tCreated affinities for "+((int)word)+" out of "+kmersArray.length+" "+kmerLength+"-mers so far...");
			// }
		}

		return (newKmerToAffinityMatrix);
	}

	public double[] getProbeAffinities(char kmerLength,
			double[] kmerToAffinityMatrix, double[][] weights,
			boolean includeRevComps) {
		double[] probeAffinities = new double[this.numSeqIDs];
		// int maxWord = (int)Math.pow(4, kmerLength);
		int maxWord = multBy4ToPower(1, kmerLength);

		for (int aSeqIdIndex = 0; aSeqIdIndex < this.numSeqIDs; aSeqIdIndex++) {

			// add the intecept if it exists
			if (kmerToAffinityMatrix.length > maxWord) {
				probeAffinities[aSeqIdIndex] += kmerToAffinityMatrix[maxWord];
			}

			// sum over each kmer window
			for (int startPos = 0; startPos < this.charKmerMatrix[aSeqIdIndex].length; startPos++) {
				char aWordIndex = this.charKmerMatrix[aSeqIdIndex][startPos];
				double affinity = 0;

				// Positive Strand
				affinity += weights[0][startPos]
						* kmerToAffinityMatrix[aWordIndex];

				// Negative Strand
				char aRevCompIndex = revCompMatrix[aWordIndex];
				if (includeRevComps) {
					affinity += weights[1][startPos]
							* kmerToAffinityMatrix[aRevCompIndex];

					// self-reverse-complement words have been counted twice
					if (aRevCompIndex == aWordIndex) {
						affinity = affinity / 2.0;
					}
				}

				probeAffinities[aSeqIdIndex] += affinity;

				// boolean isSelfRevComp = false;
				// if (aRevCompIndex == aWordIndex) {
				// isSelfRevComp = true;
				// }

			}
		}
		return (probeAffinities);
	}

	public static Table getRevCompPalindromes(int kmerLength, Table aTable,
			int cutoff) {
		Table palindromesTable = new Table();

		int count = 0;
		for (int i = 0; i < aTable.rows(); i++) {

			List<Object> row = aTable.getRow(i);
			Long word = ((Long) row.get(0)).longValue();
			long revCompIndex = getRevComp(word, kmerLength);
			if (word == revCompIndex) {
				palindromesTable.add(row);
				count++;
				if (count >= cutoff) {
					break;
				}
			}
		}
		return (palindromesTable);
	}

	public static String toString(int kmerLength, Table aTable) {
		return (toString(kmerLength, aTable, "\t", false, 0, aTable.rows()));
	}

	public static String toString(int kmerLength, Table aTable,
			String delimiter, boolean rowLabels, int start, int end) {

		if (start < 0 || start > end) {
			System.err.println("Invalid Indexes: start=" + start + " end="
					+ end + ".");
			return (null);
		}

		StringBuffer output = new StringBuffer();
		for (int i = start; i < end; i++) {

			if (i < aTable.rows()) {

				List<Object> row = aTable.getRow(i);
				for (int j = 0; j < row.size(); j++) {
					if (j > 0) {
						output.append(delimiter);
						output.append(row.get(j));
					} else {
						if (rowLabels) {
							output.append((i + 1) + "" + delimiter);
						}
						String wordString = getWordString(
								((Long) row.get(j)).longValue(), kmerLength);
						output.append(wordString);
					}
				}
				output.append("\n");

			}
		}
		return (output.toString());
	}

	public static String getSortedTableString(int kmerLength,
			double[] kmerToAffinityMatrix) {
		Table sortedTable = getSortedTable(kmerLength, kmerToAffinityMatrix);
		return (toString(kmerLength, sortedTable));
	}

	public static Table getSortedTable(int kmerLength,
			double[] kmerToAffinityMatrix) {
		// int maxWord = (int)Math.pow(4, kmerLength);

		int maxWord = multBy4ToPower(1, kmerLength);

		Table kmerAffinitiesTable = new Table();

		// be sure to use "int word" and not "char word" or infinite loop for
		// 8mers
		for (int word = 0; word < maxWord; word++) {
			kmerAffinitiesTable.add(Arrays.asList(new Long(word), new Double(
					kmerToAffinityMatrix[word])));
		}

		kmerAffinitiesTable.sort(1); // sort by ascending relativeAffinity
										// values
		kmerAffinitiesTable.reverse(); // get descending order

		return (kmerAffinitiesTable);
	}

}
