import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

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
import org.biojava.bio.mydp.GridLayoutFrame;
import org.biojava.bio.mydp.GraphicsTools;

//
// This uses a "sparse array" implementation and will be inefficient for very large sequences
//

public class WordCounts {

    private static String tabs = "\\t+";
    private static Pattern tabsPattern = Pattern.compile(tabs);

    public class MutableInt{
        int mInt;

        public MutableInt(int i) {
            mInt = i;
        }

        public String toString() {
            return(Integer.toString(mInt));
        }
    }

    private static String alphabetName = "DNA";
    private Alphabet alphabet;
    private SymbolTokenization symbolTokenization = null;
    private int kmerSize = -1;

    // make sure all entries extend AbstractSymbolList
    LinkedHashMap<SymbolList, MutableInt> wordCountsMap = null;

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    public WordCounts(int aKmerSize, int initSize){
        setAlphabet(alphabetName);
        this.kmerSize = aKmerSize;
        wordCountsMap = new LinkedHashMap<SymbolList, MutableInt>(initSize);
    }

    public WordCounts(int aKmerSize, int initSize, String anAlphabetName, Alphabet anAlphabet, SymbolTokenization aSymTok) {
        setAlphabet(alphabetName, anAlphabet, aSymTok);
        this.kmerSize = aKmerSize;
        wordCountsMap = new LinkedHashMap<SymbolList, MutableInt>(initSize);
    }

    // init from a saved WordCounts File
    public WordCounts(String fileName){
        setAlphabet(alphabetName);
        getCounts(fileName);
    }

    // init from a saved WordCounts File
    public WordCounts(String fileName, String anAlphabetName, Alphabet anAlphabet, SymbolTokenization aSymTok) {
        setAlphabet(alphabetName, anAlphabet, aSymTok);
        getCounts(fileName);
    }

    // init and calculate the word counts of a file
    public WordCounts(String fileName, int aKmerSize, int initSize, boolean windowed){
        this(aKmerSize, initSize);

        try {
            getCounts(fileName, aKmerSize, windowed);
        }
        catch(Exception ex) {
            ex.printStackTrace(System.err);
        }
    }

    // init and calculate the word counts of a file
    public WordCounts(String fileName, int aKmerSize, int initSize, boolean windowed, String anAlphabetName, Alphabet anAlphabet, SymbolTokenization aSymTok) {
        this(aKmerSize, initSize, anAlphabetName, anAlphabet, aSymTok);

        try {
            getCounts(fileName, aKmerSize, windowed);
        }
        catch(Exception ex) {
            ex.printStackTrace(System.err);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////


    public void setAlphabet(String anAlphabetName, Alphabet anAlphabet, SymbolTokenization aSymTok) {
        this.alphabetName = anAlphabetName;
        this.alphabet = anAlphabet;
        this.symbolTokenization = aSymTok;
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
        }
        catch(Exception ex) {
            ex.printStackTrace(System.err);
        }
    }

    public void write(String fileName) {
        FileTools.write(this.wordCountsMap, fileName, "\t", "");
    }

    public int getCount(SymbolList aSymbolList) {
        MutableInt mutableInt = wordCountsMap.get(aSymbolList);
        if (mutableInt == null) {
            return(0);
        }
        return(mutableInt.mInt);
    }

    public void putCount(SymbolList aSymbolList, int count) {
        //wordCountsMap.put(aSymbolList, new MutableInt(count));
        MutableInt mutableInt = wordCountsMap.get(aSymbolList);
        if (mutableInt == null) {
            wordCountsMap.put(aSymbolList, new MutableInt(count));
        }
        else {
            mutableInt.mInt = count;
        }
    }

    public void incrementCount(SymbolList aSymbolList, int count) {
        MutableInt mutableInt = wordCountsMap.get(aSymbolList);
        if (mutableInt == null) {
            wordCountsMap.put(aSymbolList, new MutableInt(count));
        }
        else {
            mutableInt.mInt += count;
        }
    }

    // getProbSum is a very fast calulation of the sum of all the probability scores of a weight matrix over the
    // entire sequence that is represented by this WordCounts object
    public double getProbSum(WeightMatrix aWeightMatrix) {
        double overallScore = 0;
        double aScore;
        int numColumns = aWeightMatrix.columns();

        try {

            if (numColumns != this.kmerSize) {
                System.out.println("Error : weight matrix length is not the same as the words counts length.\n");
                return(-1);
            }

            // score each word
            for (SymbolList aWord : wordCountsMap.keySet()) {
                MutableInt count = wordCountsMap.get(aWord);

                // score each numColumns-bp window on the sequence
                aScore = 0;
                for (int c = 0; c < aWeightMatrix.columns(); c++) {
                    aScore += Math.log(aWeightMatrix.getColumn(c).getWeight(aWord.symbolAt(c)));
                }
                aScore = Math.exp(aScore);

                // multiply by count to recover original sequence
                overallScore += count.mInt * aScore;
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }

        return(overallScore);
    }

    public int[] getCountsArray(SymbolList[] allWords) {
        int[] allCounts = new int[allWords.length];
        for (int i=0; i < allWords.length; i++) {
            allCounts[i] = getCount(allWords[i]);
        }
        return(allCounts);
    }

    public void getCounts(String aFilePathName) {
        try {
            String keyValuePairs[] = FileTools.readStrings(aFilePathName);

            if (wordCountsMap == null) {
                int initSize = keyValuePairs.length*2;
                LinkedHashMap<SymbolList, MutableInt> wordCountsMap = new LinkedHashMap<SymbolList, MutableInt>(initSize);
            }

            if (this.kmerSize==-1 && keyValuePairs!=null) {
                String[] keyValuePairArray = tabsPattern.split(keyValuePairs[0]);
                this.kmerSize = ((String)keyValuePairArray[0]).length();
            }

            for (String keyValuePair : keyValuePairs) {
                String[] keyValuePairArray = tabsPattern.split(keyValuePair);
                SymbolList aWord = new SimpleSymbolList(symbolTokenization, keyValuePairArray[0]);
                incrementCount(aWord, Integer.parseInt(keyValuePairArray[1]));
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getCounts(String fileName, int kmerSize, String alphabetName, boolean windowed)
        throws IllegalAlphabetException, IllegalSymbolException, BioException {
        setAlphabet(alphabetName);
        getCounts(fileName, kmerSize, windowed);
    }

    public void getCounts(String fileName, int kmerSize, boolean windowed)
        throws IllegalAlphabetException, IllegalSymbolException, BioException{

        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(fileName));
        }
        catch(Exception ex) {
            ex.printStackTrace(System.err);
        }

        //         RichStreamReader sr = new  RichStreamReader(
        //                 br, format, toke,
        //                 RichSequenceBuilderFactory.THRESHOLD,
        //                 RichObjectFactory.getDefaultNamespace());

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
        getCounts(seqIter, kmerSize, windowed);
    }

    public void getCounts(SequenceIterator iter)
        throws IllegalAlphabetException, IllegalSymbolException, BioException{
        getCounts(this.asIterator(iter), 1, false);
    }

    public void getCounts(SequenceIterator iter, int kmerSize, boolean windowed)
        throws IllegalAlphabetException, IllegalSymbolException, BioException{
        getCounts(this.asIterator(iter), kmerSize, windowed);
    }

    public void getCounts(Iterator<? extends SymbolList> iter, int kmerSize, boolean windowed)
        throws IllegalAlphabetException, IllegalSymbolException, BioException{

        if(kmerSize > 1){
            iter = this.nmerView(iter, kmerSize, windowed);
        }

        while(iter.hasNext()){

            SymbolList symList = iter.next();

            //System.out.print("\nkeyString="+keyString+" wordIndex="+wordIndex);
            incrementCount(symList, 1);
        }

    }

    /**
     * Takes each <CODE>SymbolList</CODE> from the <CODE>Iterator</CODE> and applies
     * a view to it. The view can be windowed (eg codons) or
     * sliding (eg overlapping dimers)
     * @param iter The input iterator
     * @param nmerSize The size of the window eg 3 for codons.
     * If the size is less than 2 then you get back
     * the original <CODE>Iterator</CODE>
     * @param windowed true if you want non-overlapping nmers (eg codons),
     * false if you want them to overlap.
     * @return An <CODE>Iterator</CODE> over <CODE>SymbolLists</CODE> with the
     * desired view applied. <B>You cannot call <code>remove()</code> on this iterator!</B>
     */
    public Iterator<SymbolList> nmerView(Iterator<? extends SymbolList> iter, int nmerSize, boolean windowed){

        if(nmerSize < 2) return (Iterator<SymbolList>)iter;

        final Iterator<? extends SymbolList> it = iter;
        final int size = nmerSize;
        final boolean w = windowed;

        return new Iterator<SymbolList>(){

            public boolean hasNext(){
                return it.hasNext();
            }

            public SymbolList next() {
                try{
                    SymbolList source = it.next();
                    if(w) {
                        // extends AbstractSymbolList
                        return SymbolListViews.windowedSymbolList(source, size);
                    }
                    else {
                        // extends AbstractSymbolList
                        return SymbolListViews.orderNSymbolList(source, size);
                    }
                }
                catch(BioException e){
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
