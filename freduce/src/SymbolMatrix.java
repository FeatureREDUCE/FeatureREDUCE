/*
 * Todd Riley
 */



import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import java.util.Set;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.LinkedHashMap;
import java.io.*;
import java.lang.Math;
import java.lang.String;
import java.lang.StringBuffer;
import java.util.ArrayList;
import java.util.List;

import org.biojava.utils.*;
import org.biojava.bio.*;
import org.biojava.bio.gui.sequence.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.symbol.*;
import org.biojava.bio.gui.sequence.FeatureLabelRenderer.TypeLabelMaker;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dist.Distribution;
import java.util.Collection;
import java.util.Iterator;
import java.util.Collections;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dist.DistributionFactory.DefaultDistributionFactory;

import org.biojava.bio.seq.io.SymbolTokenization;
import java.util.NoSuchElementException;
import java.text.NumberFormat;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

import org.biojava.bio.dp.WeightMatrix;

import org.biojava.bio.mydp.StringTools;
import org.biojava.bio.mydp.FileTools;
import org.biojava.bio.mydp.MathTools;
import org.biojava.bio.mydp.WeightMatrixTools;
import org.biojava.bio.mydp.FeatureKey;

import gnu.trove.map.hash.TLongShortHashMap;
import gnu.trove.map.hash.TIntShortHashMap;
import gnu.trove.set.hash.TLongHashSet;


public class SymbolMatrix implements Serializable {

    SymbolMatrixData symbolData = null;
    Table allSeqIdIndexes = null;

    private static String tabs = "\\t+";
    private static Pattern tabsPattern = Pattern.compile(tabs);

    private static String alphabetName = "DNA";
    private Alphabet alphabet;
    private SymbolTokenization symbolTokenization = null;
    private Packing packing;

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    public SymbolMatrix(String[] aSeqIDs){
        setAlphabet(alphabetName);

        initSymbolData(aSeqIDs);
    }

    public SymbolMatrix(String[] aSeqIDs, String anAlphabetName, Alphabet anAlphabet, SymbolTokenization aSymTok, Packing aPacking) {
        setAlphabet(alphabetName, anAlphabet, aSymTok, aPacking);

        initSymbolData(aSeqIDs);
    }

    // init from a saved SymbolMatrix File
    public SymbolMatrix(String anAlphabetName, Alphabet anAlphabet, SymbolTokenization aSymTok, String fileName, Packing aPacking) {
        setAlphabet(alphabetName, anAlphabet, aSymTok, aPacking);

        symbolData = readSerializedFile(fileName);

//         System.out.println("\noriginalSeqIDs.length() = "+symbolData.originalSeqIDs.length);
//         System.out.println("seqIDs.length() = "+symbolData.seqIDs.length);
//         System.out.println("paskedSymbolLists.length() = "+symbolData.packedSymbolLists.length);
    }

    private void initSymbolData(String[] aSeqIDs) {

        symbolData = new SymbolMatrixData();

        if (aSeqIDs != null) {
            symbolData.seqIDs = aSeqIDs;
            symbolData.packedSymbolLists = new PackedSymbolList[aSeqIDs.length];
        }

    }

    public void setAlphabet(String anAlphabetName, Alphabet anAlphabet, SymbolTokenization aSymTok, Packing aPacking) {
        this.alphabetName = anAlphabetName;
        this.alphabet = anAlphabet;
        this.symbolTokenization = aSymTok;
        this.packing = aPacking;
    }

    public void setAlphabet(String anAlphabetName) {
        try {
            this.alphabetName = anAlphabetName;
            try{
                this.alphabet = AlphabetManager.alphabetForName(alphabetName);
            }
            catch(NoSuchElementException ex){
                //try it upper case
                this.alphabet = AlphabetManager.alphabetForName(alphabetName.toUpperCase());
            }
            this.symbolTokenization = alphabet.getTokenization("token");
            this.packing = PackingFactory.getPacking((FiniteAlphabet)this.alphabet, false);
        }
        catch(Exception ex) {
            ex.printStackTrace(System.err);
        }
    }

    public int getSeqLength() {
        if ((symbolData.packedSymbolLists != null) && (symbolData.packedSymbolLists[0] != null)) {
            return(symbolData.packedSymbolLists[0].length());
        }
        return(0);
    }

    public int getNumSeqIDs() {
        return(symbolData.seqIDs.length);
    }

    public String getSeqID(int index) {
        return(symbolData.seqIDs[index]);
    }

    public String[] getSeqIDs() {
        return(symbolData.seqIDs);
    }

    public String[] getOriginalSeqIDs() {
        return(symbolData.originalSeqIDs);
    }

    public SymbolList getSymbolList(int index) {
        return(symbolData.packedSymbolLists[index]);
    }

    public void setOrderedSubset(String[] orderedSeqIDs) {

        SymbolMatrixData newSymbolData = new SymbolMatrixData();
        newSymbolData.seqIDs = orderedSeqIDs;
        newSymbolData.originalSeqIDs = symbolData.seqIDs;

        // System.out.print("\norderedSeqIDs[0]="+orderedSeqIDs[0]+"; orderedSeqIDs.length="+orderedSeqIDs.length);
        // System.out.print("\nseqIDs[0]="+symbolData.seqIDs[0]+"; seqIDs.length="+symbolData.seqIDs.length);

        // Table has rows of (sedID, index)
        Table seqIDsTable = new Table(symbolData.seqIDs);

        // sort by seqID
        seqIDsTable.sort(0);

        newSymbolData.packedSymbolLists = new PackedSymbolList[orderedSeqIDs.length];

        for (int i=0; i < orderedSeqIDs.length; i++) {
            String aSeqID = orderedSeqIDs[i];
            int sortedIndex = seqIDsTable.binarySearch(0, aSeqID);
            int oldIndex = (Integer)((List<Object>)seqIDsTable.getRow(sortedIndex)).get(1);
            newSymbolData.packedSymbolLists[i] = symbolData.packedSymbolLists[oldIndex];

            // newSymbolData.packedSymbolLists[i] = symbolData.packedSymbolLists[i];
        }

        symbolData = newSymbolData;
    }

    // For a WeightMatrix
    // returns the array of probability sums by scoring every kmer in the probe sequence, for each probe sequence
    public double[] getProbSums(
        WeightMatrix aPosStrandPWM,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double eToMu,
        double nonSpecKa) {

        try {
            WeightMatrix aNegStrandPWM = null;
            double wmDistance = 1.0;
            if (strand == WeightMatrixTools.BindingStrand.NEG || strand == WeightMatrixTools.BindingStrand.BOTH) {
                aNegStrandPWM = WeightMatrixTools.reverseComplement(aPosStrandPWM, DNATools.complementTable());

                // If the PSAM is an reverse-complement palindrome, then we should look at only the positve
                // or negative strand (BUT NOT BOTH)
//                 boolean isRevCompPalindrome = WeightMatrixTools.areEmissionSpectraEqual(aPosStrandPWM, aNegStrandPWM);
//                 if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.BOTH)) {
//                     //out.println("\n\nSince the motif is self reverse-complement, we should look only at the positive strand.");
//                     strand = WeightMatrixTools.BindingStrand.POS;
//                 }

            }
            if (strand == WeightMatrixTools.BindingStrand.BOTH) {
                wmDistance = WeightMatrixTools.similarity(aPosStrandPWM, aNegStrandPWM);
            }

            double[] probSums = new double[symbolData.seqIDs.length];

            for (int aSeqIdIndex=0; aSeqIdIndex < symbolData.seqIDs.length; aSeqIdIndex++) {

                probSums[aSeqIdIndex] = WeightMatrixTools.score(
                    aPosStrandPWM,
                    aNegStrandPWM,
                    symbolData.packedSymbolLists[aSeqIdIndex],
                    strand,
                    calc,
                    eToMu,
                    nonSpecKa,
                    wmDistance);

            }
            return(probSums);
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(null);
    }

    // For a WeightMatrix
    // includes weights
    // returns the array of probability sums by scoring every kmer in the probe sequence, for each probe sequence
    // then it scales each probability sum by the weight
    public double[] getProbSums(
        WeightMatrix aPosStrandPWM,
        double[][] weights,
        boolean[][] mandatoryColumns,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double eToMu,
        double nonSpecKa) {

        try {
            WeightMatrix aNegStrandPWM = null;
            double wmDistance = 1.0;
            if (strand == WeightMatrixTools.BindingStrand.NEG || strand == WeightMatrixTools.BindingStrand.BOTH) {
                aNegStrandPWM = WeightMatrixTools.reverseComplement(aPosStrandPWM, DNATools.complementTable());

                // If the PSAM is an reverse-complement palindrome, then we should look at only the positve
                // or negative strand (BUT NOT BOTH)
//                 boolean isRevCompPalindrome = WeightMatrixTools.areEmissionSpectraEqual(aPosStrandPWM, aNegStrandPWM);
//                 if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.BOTH)) {
//                     //out.println("\n\nSince the motif is self reverse-complement, we should look only at the positive strand.");
//                     strand = WeightMatrixTools.BindingStrand.POS;
//                 }

            }
            if (strand == WeightMatrixTools.BindingStrand.BOTH) {
                wmDistance = WeightMatrixTools.similarity(aPosStrandPWM, aNegStrandPWM);
            }

            double[] probSums = new double[symbolData.seqIDs.length];

            for (int aSeqIdIndex=0; aSeqIdIndex < symbolData.seqIDs.length; aSeqIdIndex++) {

                probSums[aSeqIdIndex] = WeightMatrixTools.score(
                    aPosStrandPWM,
                    aNegStrandPWM,
                    symbolData.packedSymbolLists[aSeqIdIndex],
                    weights,
                    mandatoryColumns,
                    strand,
                    calc,
                    eToMu,
                    nonSpecKa,
                    wmDistance);

            }
            return(probSums);
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(null);
    }

    // For a WeightMatrix
    // includes start position
    // returns the array of probabilities by scoring the kmer in the probe sequence that starts at position startPos, for each probe sequence
    public double[] getProbs(
        WeightMatrix aPosStrandPWM,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double eToMu,
        double nonSpecKa,
        int startPos) {

        try {
            WeightMatrix aNegStrandPWM = null;
            double wmDistance = 1.0;
            if (strand == WeightMatrixTools.BindingStrand.NEG || strand == WeightMatrixTools.BindingStrand.BOTH) {
                aNegStrandPWM = WeightMatrixTools.reverseComplement(aPosStrandPWM, DNATools.complementTable());

                // If the PSAM is an reverse-complement palindrome, then we should look at only the positve
                // or negative strand (BUT NOT BOTH)
//                 boolean isRevCompPalindrome = WeightMatrixTools.areEmissionSpectraEqual(aPosStrandPWM, aNegStrandPWM);
//                 if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.BOTH)) {
//                     //out.println("\n\nSince the motif is self reverse-complement, we should look only at the positive strand.");
//                     strand = WeightMatrixTools.BindingStrand.POS;
//                 }

            }
            if (strand == WeightMatrixTools.BindingStrand.BOTH) {
                wmDistance = WeightMatrixTools.similarity(aPosStrandPWM, aNegStrandPWM);
            }

            double[] probs = new double[symbolData.seqIDs.length];

            for (int aSeqIdIndex=0; aSeqIdIndex < symbolData.seqIDs.length; aSeqIdIndex++) {

                probs[aSeqIdIndex] = WeightMatrixTools.score(
                    aPosStrandPWM,
                    aNegStrandPWM,
                    symbolData.packedSymbolLists[aSeqIdIndex],
                    strand,
                    calc,
                    eToMu,
                    nonSpecKa,
                    wmDistance,
                    startPos);

            }
            return(probs);
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(null);
    }

    // For a WeightMatrix
    // returns the matrix of probabilities by scoring the kmer in the probe sequence that starts at position startPos, for each probe sequence
    public double[][] getProbs(
        WeightMatrix aPosStrandPWM,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double eToMu,
        double nonSpecKa) {

        try {
            WeightMatrix aNegStrandPWM = null;
            double wmDistance = 1.0;
            if (strand == WeightMatrixTools.BindingStrand.NEG || strand == WeightMatrixTools.BindingStrand.BOTH) {
                aNegStrandPWM = WeightMatrixTools.reverseComplement(aPosStrandPWM, DNATools.complementTable());

                // If the PSAM is an reverse-complement palindrome, then we should look at only the positve
                // or negative strand (BUT NOT BOTH)
//                 boolean isRevCompPalindrome = WeightMatrixTools.areEmissionSpectraEqual(aPosStrandPWM, aNegStrandPWM);
//                 if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.BOTH)) {
//                     //out.println("\n\nSince the motif is self reverse-complement, we should look only at the positive strand.");
//                     strand = WeightMatrixTools.BindingStrand.POS;
//                 }

            }
            if (strand == WeightMatrixTools.BindingStrand.BOTH) {
                wmDistance = WeightMatrixTools.similarity(aPosStrandPWM, aNegStrandPWM);
            }

            int seqLengths = symbolData.packedSymbolLists[0].length();
            int motifLength = aPosStrandPWM.columns();
            int numMotifWindows = seqLengths - motifLength + 1;
            double[][] probs = new double[symbolData.seqIDs.length][numMotifWindows];

            for (int aSeqIdIndex=0; aSeqIdIndex < symbolData.seqIDs.length; aSeqIdIndex++) {

                for (int startPos=0; startPos < numMotifWindows; startPos++) {

                    probs[aSeqIdIndex][startPos] = WeightMatrixTools.score(
                        aPosStrandPWM,
                        aNegStrandPWM,
                        symbolData.packedSymbolLists[aSeqIdIndex],
                        strand,
                        calc,
                        eToMu,
                        nonSpecKa,
                        wmDistance,
                        startPos);
                }
            }
            return(probs);
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(null);
    }

    // For a WeightMatrix
    // includes weights
    // returns the matrix of probabilities by scoring the kmer in the probe sequence that starts at position startPos, for each probe sequence
    public double[][] getProbs(
        WeightMatrix aPosStrandPWM,
        double[][] weights,
        boolean[][] mandatoryColumns,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double eToMu,
        double nonSpecKa) {

        try {
            WeightMatrix aNegStrandPWM = null;
            double wmDistance = 1.0;
            if (strand == WeightMatrixTools.BindingStrand.NEG || strand == WeightMatrixTools.BindingStrand.BOTH) {
                aNegStrandPWM = WeightMatrixTools.reverseComplement(aPosStrandPWM, DNATools.complementTable());

                // If the PSAM is an reverse-complement palindrome, then we should look at only the positve
                // or negative strand (BUT NOT BOTH)
//                 boolean isRevCompPalindrome = WeightMatrixTools.areEmissionSpectraEqual(aPosStrandPWM, aNegStrandPWM);
//                 if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.BOTH)) {
//                     //out.println("\n\nSince the motif is self reverse-complement, we should look only at the positive strand.");
//                     strand = WeightMatrixTools.BindingStrand.POS;
//                 }

            }
            if (strand == WeightMatrixTools.BindingStrand.BOTH) {
                wmDistance = WeightMatrixTools.similarity(aPosStrandPWM, aNegStrandPWM);
            }

            int seqLengths = symbolData.packedSymbolLists[0].length();
            int motifLength = aPosStrandPWM.columns();
            int numMotifWindows = seqLengths + motifLength - 1;
            double[][] probs = new double[symbolData.seqIDs.length][numMotifWindows];

            for (int aSeqIdIndex=0; aSeqIdIndex < symbolData.seqIDs.length; aSeqIdIndex++) {

                for (int startPos = -1 * (motifLength - 1); startPos <= seqLengths - 1; startPos++) {

                    probs[aSeqIdIndex][startPos] = WeightMatrixTools.score(
                        aPosStrandPWM,
                        aNegStrandPWM,
                        symbolData.packedSymbolLists[aSeqIdIndex],
                        weights,
                        mandatoryColumns,
                        strand,
                        calc,
                        eToMu,
                        nonSpecKa,
                        wmDistance,
                        startPos);
                }
            }
            return(probs);
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(null);
    }

    // For a FeaturedWeightMatrix
    // returns the array of probability sums by scoring every kmer in the probe sequence, for each probe sequence
    public double[] getProbSums(
        FeaturedWeightMatrix anFsam,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc aCalc,
        double eToMu,
        double nonSpecKa,
        double threshold) {

        double[] probSums = new double[symbolData.seqIDs.length];

        for (int aSeqIdIndex=0; aSeqIdIndex < symbolData.seqIDs.length; aSeqIdIndex++) {

            anFsam.setStrand(strand);
            anFsam.setCalc(aCalc);
            anFsam.setEToMu(eToMu);
            //anFsam.setNonSpecKa(nonSpecKa);

            probSums[aSeqIdIndex] = anFsam.score(
                symbolData.packedSymbolLists[aSeqIdIndex],
                (FeatureKey)null);
                //, threshold);

        }
        return(probSums);
    }

    // For a FeaturedWeightMatrix
    // includes weights
    // returns the array of probability sums by scoring every kmer in the probe sequence, for each probe sequence
    // then it scales each probability sum by the weight
    public double[] getProbSums(
        FeaturedWeightMatrix anFsam,
        double[][] weights,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc aCalc,
        double eToMu,
        double nonSpecKa,
        double threshold) {

        double[] probSums = new double[symbolData.seqIDs.length];

        for (int aSeqIdIndex=0; aSeqIdIndex < symbolData.seqIDs.length; aSeqIdIndex++) {

            anFsam.setStrand(strand);
            anFsam.setCalc(aCalc);
            anFsam.setEToMu(eToMu);
            //anFsam.setNonSpecKa(nonSpecKa);

            probSums[aSeqIdIndex] = anFsam.score(
                symbolData.packedSymbolLists[aSeqIdIndex],
                weights,
                (FeatureKey)null);
                //, threshold);

        }
        return(probSums);
    }

    // For a FeaturedWeightMatrix
    // includes start position
    // returns the array of probabilities by scoring the kmer in the probe sequence that starts at position startPos, for each probe sequence
    public double[] getProbs(
        FeaturedWeightMatrix anFsam,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc aCalc,
        double eToMu,
        double nonSpecKa,
        double threshold,
        int startPos) {

        double[] probs = new double[symbolData.seqIDs.length];

        for (int aSeqIdIndex=0; aSeqIdIndex < symbolData.seqIDs.length; aSeqIdIndex++) {

            anFsam.setStrand(strand);
            anFsam.setCalc(aCalc);
            anFsam.setEToMu(eToMu);
            //anFsam.setNonSpecKa(nonSpecKa);

            probs[aSeqIdIndex] = anFsam.score(
                symbolData.packedSymbolLists[aSeqIdIndex],
                (FeatureKey)null,
//                 threshold,
                startPos);

        }
        return(probs);
    }

    public void writeSerializedFile(String aFilePathName) {
        FileTools.writeSerializedFile(symbolData, aFilePathName);
    }


    // seqIdStrings must have same iteration order as the matrix rows!!
    public void writeTable(String aFilePathName) {
        try {
            BufferedWriter outBuffer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(aFilePathName)));

            for (int aSeqIdIndex=0; aSeqIdIndex < symbolData.seqIDs.length; aSeqIdIndex++) {
                String seqString = symbolData.packedSymbolLists[aSeqIdIndex].seqString();
                outBuffer.write(symbolData.seqIDs[aSeqIdIndex] + "\t" + seqString + "\n");
            }
            outBuffer.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public SymbolMatrixData readSerializedFile(String filePathName) {
        return((SymbolMatrixData)FileTools.readSerializedFile(filePathName));
    }

    public void addSymbolLists(String fileName, String alphabetName)
        throws IllegalAlphabetException, IllegalSymbolException, BioException {
        setAlphabet(alphabetName);
        addSymbolLists(fileName);
    }

    public void addSymbolLists(String fileName)
        throws IllegalAlphabetException, IllegalSymbolException, BioException{

        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(fileName));
        }
        catch(Exception ex) {
            ex.printStackTrace(System.err);
        }

        // Replaced sr with a symbollist (or a sequence)
        SequenceIterator seqIter = null;
        Sequence aSequence = null;

        try {
            seqIter = SeqIOTools.readFasta(br, symbolTokenization);
            //aSequence = seqIter.nextSequence();
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        addSymbolLists(seqIter);
    }

    public void addSymbolLists(SequenceIterator iter)
        throws IllegalAlphabetException, IllegalSymbolException, BioException{
        addSymbolLists(this.asIterator(iter));
    }

    public void addSymbolLists(Iterator<Sequence> seqIter)
        throws IllegalAlphabetException, IllegalSymbolException, BioException{

        int seqNameIndex = 0;
        while (seqIter.hasNext()) {
            Sequence sequence = seqIter.next();
            addSymbolList(seqNameIndex, sequence);
            seqNameIndex++;
        }
    }

    public void addSymbolList(int aSeqIDIndex, Sequence aSeq)
        throws IllegalAlphabetException, IllegalSymbolException, BioException{
            String seqName = aSeq.getName();
            symbolData.seqIDs[aSeqIDIndex] = seqName;
            symbolData.packedSymbolLists[aSeqIDIndex] = new PackedSymbolList(this.packing, aSeq);
    }

    /**
     * Makes a <code>SequenceIterator</code> look like an
     * <code>Iterator {@code <Sequence>}</code>
     * @param iter The <CODE>SequenceIterator</CODE>
     * @return An <CODE>Iterator</CODE> that returns only <CODE>Sequence</CODE>
     * objects. <B>You cannot call <code>remove()</code> on this iterator!</B>
     */
    public Iterator<Sequence> asIterator(SequenceIterator iter){
        final SequenceIterator it = iter;

        return new Iterator<Sequence>(){
            public boolean hasNext(){
                return it.hasNext();
            }
            public Sequence next() {
                try{
                    return it.nextSequence();
                }
                catch(BioException e) {
                    NoSuchElementException ex = new NoSuchElementException();
                    ex.initCause(e);
                    throw ex;
                }
            }
            public void remove(){
                throw new UnsupportedOperationException();
            }
        };
    }



}
