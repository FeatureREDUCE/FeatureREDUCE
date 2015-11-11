import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import java.util.Hashtable;
import java.util.Arrays;
import java.util.ArrayList;
import java.io.*;
import java.lang.Math;
import java.lang.String;
import java.lang.StringBuffer;
import java.util.TreeSet;

import org.biojava.utils.*;
import org.biojava.bio.*;
import org.biojava.bio.gui.sequence.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.symbol.*;
import org.biojava.bio.gui.sequence.FeatureLabelRenderer.TypeLabelMaker;
import org.biojava.bio.seq.io.SeqIOTools;

public class RocData{

    String dataName;

    //TreeSet positiveSet;
    //TreeSet negativeSet;
    double[] positiveSet;
    double[] negativeSet;

    double minScore;
    double maxScore;
    double affineThresholdOffset;
    double affineThresholdScalar;
    double thresholdStepSize;
    double datasetSize;
    int bins;

    double[][] rocArray;
    double auc = 0;

    public RocData(String aDataName, String aModelFile, int aBins) {
        String[] dataLines = getStringsFromFile(aModelFile);
        double[][] scores = convertLinesToDoubleMatrix(dataLines);
        init(aDataName, scores, aBins);
    }

    public RocData(String aDataName, String aModelFile, double aThresholdStepSize) {
        String[] dataLines = getStringsFromFile(aModelFile);
        double[][] scores = convertLinesToDoubleMatrix(dataLines);
        int aBins = (int) Math.round((double) 1.0 / aThresholdStepSize);
        init(aDataName, scores, aBins);
    }

    public RocData(String aDataName, double[] positiveSetScores, double[] negativeSetScores, int aBins) {
        init(aDataName, positiveSetScores, negativeSetScores, aBins);
    }

    public void init(String aDataName, double[][] scores, int aBins) {

        ArrayList positiveSetScores = new ArrayList();
        ArrayList negativeSetScores = new ArrayList();

        for (int i=0; i < scores[0].length; i++) {
            if (scores[0][i] == 0) {
                negativeSetScores.add(new Double(scores[1][i]));
            }
            else {
                positiveSetScores.add(new Double(scores[1][i]));
            }
        }

        init(aDataName, doubleObjsToPrimitives((Double[])positiveSetScores.toArray(new Double[0])), doubleObjsToPrimitives((Double[])negativeSetScores.toArray(new Double[0])), aBins);

    }


    public void init(String aDataName, double[] positiveSetScores, double[] negativeSetScores, int aBins) {
        dataName = aDataName;
        bins = aBins;
        thresholdStepSize = (double) 1.0 / aBins;

        //positiveSet = new TreeSet(positiveSetScores);
        //negativeSet = new TreeSet(negativeSetScores);
        positiveSet = positiveSetScores;
        negativeSet = negativeSetScores;
        Arrays.sort(positiveSet);
        Arrays.sort(negativeSet);

        //minScore = Math.min(positiveSet.first(), negativeSet.first());
        //maxScore = Math.max(positiveSet.last(), negativeSet.last());
        minScore = Math.min(positiveSet[0], negativeSet[0]);
        maxScore = Math.max(positiveSet[positiveSet.length-1], negativeSet[negativeSet.length-1]);

        affineThresholdOffset = minScore;
        affineThresholdScalar = maxScore - minScore;

        //datasetSize = positiveSet.size() + negativeSet.size();
        datasetSize = positiveSet.length + negativeSet.length;

        rocArray = new double[2][bins+1];


        int currentPosIndex;
        int currentNegIndex;
        int lastPosIndex = -1;
        int lastNegIndex = -1;
        double truePosFract;
        double falsePosFract;
        int binCounter = 0;

        //With thresold go backwards from 1 to 0, so that the (FPF, TPF) values go from (0,0) to (1,1)
        for (double thresh=1.0; thresh >= -0.00009; thresh -= thresholdStepSize) {
            double modelThresh = (thresh*affineThresholdScalar) + affineThresholdOffset;

            //double truePosFract = positiveSet.tailSet(modelThresh).size() / positiveSet.length;
            //double falsePosFract = negativeSet.tailSet(modelThresh).size() / negativeSet.length;


            // Only perform binary Search if the index must change!!
            if (lastPosIndex != -1 && lastPosIndex != 0 && positiveSet[lastPosIndex-1] < modelThresh) {
                currentPosIndex = lastPosIndex;
            }
            else {
                // An improvement here would be to binary search on only the RHS of lastPosIndex
                currentPosIndex = Arrays.binarySearch(positiveSet, modelThresh);

                // if currentPosIndex >= 0, then it was found and currentPosIndex=IndexWithValue
                // if currentPosIndex < 0, then abs(currentPosIndex) - 1 == FirstIndex[ > Value]

                // If found last element exactly then +1 so that end at 0
                if (currentPosIndex == positiveSet.length - 1) {
                    currentPosIndex = positiveSet.length;
                }
                else if (currentPosIndex < 0) {
                    currentPosIndex = Math.abs(currentPosIndex) - 1;
                }

            }

            // Only perform binary Search if the index must change!!
            if (lastNegIndex != -1 && lastNegIndex != 0 && negativeSet[lastNegIndex-1] < modelThresh) {
                currentNegIndex = lastNegIndex;
            }
            else {
                // An improvement here would be to binary search on only the RHS of lastNegIndex
                currentNegIndex = Arrays.binarySearch(negativeSet, modelThresh);

                // if currentNegIndex >= 0, then it was found and currentNegIndex=indexWithThisValue
                if (currentNegIndex < 0) {
                    currentNegIndex = Math.abs(currentNegIndex) - 1;
                }
            }


            truePosFract = ((double) (positiveSet.length - currentPosIndex)) / (double)positiveSet.length;
            falsePosFract = ((double) (negativeSet.length - currentNegIndex)) / (double)negativeSet.length;

            // Should have start values of (0,0) and end values of (1,1)
            // 	    if (binCounter == 0 || binCounter == 1 || binCounter == bins-1 || binCounter == bins) {
            // 		System.out.println("thresh = "+thresh);
            // 		System.out.println("modelThresh = "+modelThresh);
            // 		System.out.println("positiveSet.length = "+positiveSet.length);
            // 		System.out.println("posSearchIndex = "+currentPosIndex);
            // 		System.out.println("negativeSet.length = "+negativeSet.length);
            // 		System.out.println("negSearchIndex = "+currentNegIndex);
            // 		System.out.println("truePosFract = "+truePosFract);
            // 		System.out.println("falsePosFract = "+falsePosFract);
            // 		System.out.println();
            // 	    }

            // Increment AUC
            if (binCounter > 0) {
                auc += getTrapezoidalArea(rocArray[0][binCounter-1], falsePosFract, rocArray[1][binCounter-1], truePosFract);
            }

            rocArray[0][binCounter] = falsePosFract;
            rocArray[1][binCounter] = truePosFract;

            // If found last element exactly then -1
            if (currentPosIndex == positiveSet.length) {
                currentPosIndex = positiveSet.length - 1;
            }

            lastPosIndex = currentPosIndex;
            lastNegIndex = currentNegIndex;
            binCounter++;
        }
    }

    private double getTrapezoidalArea(double x0, double x1, double y0, double y1) {
        double area = Math.abs(x1 - x0) * Math.min(y0, y1);
        area += .5 * Math.abs(y1 - y0) * Math.abs(x1 - x0);
        return(area);
    }


    public String[] getStringsFromFile(String filePathName) {

        java.lang.StringBuffer text = new java.lang.StringBuffer();

        try {
            java.lang.String buffer;
            BufferedReader inBuffer = new BufferedReader(new InputStreamReader(new FileInputStream(filePathName)));

            while((buffer = inBuffer.readLine()) != null) {
                if (!buffer.trim().startsWith("#")) {
                    text.append(buffer);
                    text.append("\n");
                }
            }
            inBuffer.close();
            //System.out.println(text);
        }
        catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }
        return((new String(text)).split("\\n"));
    }

    public static double[] doubleObjsToPrimitives(Double[] values) {
        double[] doubles = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            doubles[i] = values[i].doubleValue();
	    }
        return(doubles);
    }

    public static double[][] doubleObjsToPrimitives(Double[][] values) {
        double[][] doubles = new double[values.length][values[0].length];
        for (int i = 0; i < values.length; i++) {
            for (int j = 0; j < values[0].length; i++) {
                doubles[i][j] = values[i][j].doubleValue();
            }
        }
        return(doubles);
    }

    public static double[][] convertLinesToDoubleMatrix(String[] linesStrings) {

        // Find if the first line is all labels
        int start = 0;
        if (linesStrings[0].trim().startsWith("#")) {
            start = 1;
        }

        // Find out how many columns there are
        String[] strings = linesStrings[start].split("\\t+");
        int numColumns = strings.length;

        // Create Double matrix
        double[][] doubles = new double[numColumns][linesStrings.length - start];

        // Each column represents a dataset of y's
        // x[]  = doubles[0]
        // y1[] = doubles[1]
        // y2[] = doubles[2]
        // ...
        for (int i = start; i < linesStrings.length; i++) {
            String[] doubleStrings = linesStrings[i].split("\\t+");
            for (int j = 0; j < numColumns; j++) {
                doubles[j][i - start] = Double.parseDouble(doubleStrings[j]);
            }
        }
        return(doubles);
    }

}
