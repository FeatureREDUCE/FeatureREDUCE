import java.util.*;
import java.util.StringTokenizer;
import java.io.File;
import java.text.NumberFormat;
import java.util.regex.Pattern;
import java.math.BigDecimal;

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
import org.biojava.bio.symbol.Symbol;

// coefficients are as returned in R
// First coefficient is the intercept, then the rest are in independent var order
// a motif is either a SymbolList, WeightMatrix, or FeaturedWeightMatrix
class RegressionFitData {
    public double[] coefficients;
    public double alpha;
    public double beta;
    public double nonSpecKa;
    public double residuals[];

    public int[] startPositions;
    public double rSquared;
    public double stddev;
    public double rss;
    public LinkedHashMap<Object, Double> motifsWithCoeffsMap;
    public LinkedHashSet motifsWithoutCoeffs;
    public LinkedHashMap<Object, Symbol[]> motifsToModsMap;
    public LinkedHashMap<Object, boolean[][]> motifsToMandatoriesMap;

    public RegressionFitData(
        LinkedHashMap aMotifsWithCoeffsMap,
        LinkedHashSet aMotifsWithoutCoeffs,
        double[] aCoefficients,
        double[] aResiduals,
        int[] aStartPositions,
        double anRSquared,
        double anStddev,
        double anRss) {
        this(aMotifsWithCoeffsMap,
            aMotifsWithoutCoeffs,
            null,
            null,
            aCoefficients,
            aResiduals,
            aStartPositions,
            anRSquared,
            anStddev,
            anRss);
    }

    public RegressionFitData(
        LinkedHashMap aMotifsWithCoeffsMap,
        LinkedHashSet aMotifsWithoutCoeffs,
        LinkedHashMap aMotifsToModsMap,
        LinkedHashMap aMotifsToMandatoriesMap,
        double[] aCoefficients,
        double[] aResiduals,
        int[] aStartPositions,
        double anRSquared,
        double anStddev,
        double anRss) {
        this(aMotifsWithCoeffsMap,
            aMotifsWithoutCoeffs,
            aMotifsToModsMap,
            aMotifsToMandatoriesMap,
            aCoefficients,
            aResiduals,
            Double.NaN,
            Double.NaN,
            Double.NaN,
            aStartPositions,
            anRSquared,
            anStddev,
            anRss);
    }

    public RegressionFitData(
        LinkedHashMap aMotifsWithCoeffsMap,
        LinkedHashSet aMotifsWithoutCoeffs,
        LinkedHashMap aMotifsToModsMap,
        LinkedHashMap aMotifsToMandatoriesMap,
        double[] aCoefficients,
        double[] aResiduals,
        double anAlpha,
        double aBeta,
        double aNonSpecKa,
        int[] aStartPositions,
        double anRSquared,
        double anStddev,
        double anRss) {
        this.motifsWithCoeffsMap = aMotifsWithCoeffsMap;
        this.motifsWithoutCoeffs = aMotifsWithoutCoeffs;
        this.motifsToModsMap = aMotifsToModsMap;
        this.motifsToMandatoriesMap = aMotifsToMandatoriesMap;
        this.coefficients = aCoefficients;
        this.residuals = aResiduals;
        this.alpha = anAlpha;
        this.beta = aBeta;
        this.nonSpecKa = aNonSpecKa;
        this.startPositions = aStartPositions;
        this.rSquared = anRSquared;
        this.stddev = anStddev;
        this.rss = anRss;
    }

    public double getAlpha() {
        return(alpha);
    }

    public double getBeta() {
        return(beta);
    }

    public double getNonSpecKa() {
        return(nonSpecKa);
    }

    public double[] getResiduals() {
        return(residuals);
    }

    public void setResiduals(double[] aResiduals) {
        this.residuals = aResiduals;
    }

    public double getIntercept() {
        return(coefficients[0]);
    }

//     public double getIntercept(Object seedMotif) {
//         double intercept = getIntercept();
//         double seedMotifCoeff = getCoefficient(seedMotif);
//         return(intercept / seedMotifCoeff);
//     }

    public double[] getCoefficients() {
        return(coefficients);
    }

    public double getCoefficient(Object motif) {
        double motifCoeff = motifsWithCoeffsMap.get(motif);
        return(motifCoeff);
    }

    public boolean isCoefficientZero(Object motif) {
        double motifCoefficient = getCoefficient(motif);
        if (motifCoefficient == 0) {
            return(true);
        }
        return(false);
    }

    public LinkedHashMap<Object, LinkedHashSet<SymbolList>> significantMotifs(Object seedMotif, double threshold) {
        LinkedHashMap<Object, LinkedHashSet<SymbolList>> motifToHammingSphereMap = null;

        try {
            motifToHammingSphereMap = new LinkedHashMap<Object, LinkedHashSet<SymbolList>>();
            motifToHammingSphereMap.put(seedMotif, null);

            // Go through the motifs that have a coefficient
            for (Object motif : motifsWithCoeffsMap.keySet()) {
                if (Math.abs(motifsWithCoeffsMap.get(motif)) > threshold) {
                    motifToHammingSphereMap.put(motif, null);
                }
            }

            // Go through the motifs that DON'T have a coefficient and find it
            if ((motifsWithoutCoeffs != null) && (!motifsWithoutCoeffs.isEmpty())) {
                Object revCompMotif = null;
                for (Object motif : motifsWithoutCoeffs) {
                    // get the revCompMotif
                    if (motif instanceof WeightMatrix) {
                        revCompMotif = WeightMatrixTools.reverseComplement((WeightMatrix)motif, DNATools.complementTable());
                    }
                    else if (motif instanceof FeaturedWeightMatrix) {
                        //revCompMotif = null;
                        continue;
                    }
                    else { // it's a SymbolList
                        revCompMotif = DNATools.reverseComplement((SymbolList)motif);
                    }

                    // lookup the coeff for the revCompMotif
                    if ((revCompMotif != null) && (Math.abs(motifsWithCoeffsMap.get(revCompMotif)) > threshold)) {
                        motifToHammingSphereMap.put(motif, null);
                    }
                }
            }
        }

        catch (Exception ex) {
            ex.printStackTrace();
        }

        return(motifToHammingSphereMap);
    }

    public ArrayList<WeightMatrixFeature> getFeatures(Object seedMotif) {
        return(getFeatures(seedMotif, 1.0));
    }

    public ArrayList<WeightMatrixFeature> getFeatures(Object seedMotif, double scaler) {
        return(getFeatures(seedMotif, scaler, false));
    }

    public ArrayList<WeightMatrixFeature> getFeatures(Object seedMotif, double scaler, boolean additive) {
        ArrayList<WeightMatrixFeature> featuresList = null;
        try {
            featuresList = new ArrayList<WeightMatrixFeature>();

            double seedMotifCoeff = getCoefficient(seedMotif);
            double seedMotifRevCompSim = 0;
            WeightMatrix seedPwm = null;
            WeightMatrix seedPwmRevComp = null;

            if (seedMotif instanceof WeightMatrix) {
                seedPwm = (WeightMatrix)seedMotif;
                seedPwmRevComp = WeightMatrixTools.reverseComplement(seedPwm, DNATools.complementTable());
                seedMotifRevCompSim = WeightMatrixTools.similarity(seedPwm, seedPwmRevComp);
            }

            // Go through the motifs that have a coefficient
            for (Object motif : motifsWithCoeffsMap.keySet()) {

                // the motif is a WeightMatrix
                if (motif instanceof WeightMatrix) {

                    // don't include the seedMotif as a feature
                    if (motif == seedMotif) {
                        continue;
                    }

                    WeightMatrix featurePwm = (WeightMatrix)motif;
                    WeightMatrix featurePwmRevComp = WeightMatrixTools.reverseComplement(featurePwm, DNATools.complementTable());
                    double featureRevCompSim = WeightMatrixTools.similarity(featurePwm, featurePwmRevComp);

                    // don't include the feature if the coefficient is 0
                    double featureCoeff = motifsWithCoeffsMap.get(featurePwm);
                    if (featureCoeff == 0) {
                        continue;
                    }
                    //                 else {
                    //                     //double[] mandatoriesArray = motifsToMandatoriesMap.get(featurePwm)[0];
                    //                     Symbol[] modsArray = motifsToModsMap.get(featurePwm);
                    //                     try {
                    //                         for (int i=0; i < seedPwm.columns(); i++) {
                    //                             if (modsArray[i*2] != null) {
                    //                                 double w = seedPwm.getColumn(i).getWeight(modsArray[i*2]);
                    //                                 System.out.println("i="+i+"; w="+w);
                    //                                 featureCoeff = featureCoeff / w;
                    //                                 System.out.println("modifiedFeatureCoeff="+featureCoeff);
                    //                             }
                    //                         }
                    //                     }
                    //                     catch (Exception ex) {
                    //                         ex.printStackTrace();
                    //                     }
                    //                 }

                    String name = featurePwm.getName();
                    String label;
                    if (name.contains("insert")) {
                        label = "insertion";
                    }
                    else if (name.contains("delete")) {
                        label = "deletion";
                    }
                    else {
                        label = "dependency";
                    }

                    boolean isGoodFeature;
                    double relAffinity;

                    if (seedMotifCoeff == 0) {
                        seedMotifCoeff = MathTools.max(coefficients);
                    }

                    if (additive) {
                        //relAffinity = Math.abs(featureCoeff / seedMotifCoeff) * scaler;
                        relAffinity = (featureCoeff / seedMotifCoeff) * scaler;

                        if (Math.signum(seedMotifCoeff) == Math.signum(featureCoeff)) {
                            isGoodFeature = true;
                        }
                        else {
                            isGoodFeature = false;
                        }
                    }
                    else {
                        relAffinity = ((featureCoeff + seedMotifCoeff) / seedMotifCoeff) * scaler;
                        //relAffinity = ( ((featureCoeff/(featureRevCompSim+1)) + (seedMotifCoeff/(seedMotifRevCompSim+1))) / (seedMotifCoeff/(seedMotifRevCompSim+1)) ) * scaler;
                        if (relAffinity < 0) {
                            relAffinity = 0;
                        }

                        if (relAffinity >= 1.0) {
                            isGoodFeature = true;
                        }
                        else {
                            isGoodFeature = false;
                        }
                    }

                    //                 if (relAffinity < 1.5) {
                    //                     continue;
                    //                 }

//                     System.out.println("featureCoeff="+featureCoeff);
//                     System.out.println("featureRevCompSim="+featureRevCompSim);
//                     System.out.println("seedMotifCoeff="+seedMotifCoeff);
//                     System.out.println("seedMotifRevCompSim="+seedMotifRevCompSim);
//                     System.out.println("relAffinity="+relAffinity);


                    //                 if (seedMotifCoeff != 0) {
                    //                     //relAffinity = Math.abs(featureCoeff / seedMotifCoeff) * scaler;
                    //                     //relAffinity = (featureCoeff / seedMotifCoeff) * scaler;
                    //                     relAffinity = ((featureCoeff+seedMotifCoeff) / seedMotifCoeff) * scaler;

                    //                     //if (Math.signum(seedMotifCoeff) == Math.signum(featureCoeff)) {
                    //                     if (relAffinity >= 1.0) {
                    //                         isGoodFeature = true;
                    //                     }
                    //                     else {
                    //                         isGoodFeature = false;
                    //                     }
                    //                 }
                    //                 else { // seedMotifCoeff==0
                    //                     double maxCoeff = MathTools.max(coefficients);
                    //                     //relAffinity = Math.abs(featureCoeff / maxCoeff) * scaler;
                    //                     //relAffinity = (featureCoeff / maxCoeff) * scaler;
                    //                     relAffinity = ((featureCoeff+maxCoeff) / maxCoeff) * scaler;

                    //                     //if (Math.signum(featureCoeff) == 1.0) {
                    //                     if (relAffinity >= 1.0) {
                    //                         isGoodFeature = true;
                    //                     }
                    //                     else {
                    //                         isGoodFeature = false;
                    //                     }
                    //                 }

                    //double intercept = coefficients[0];
                    //double relAffinity = Math.abs((featureCoeff + intercept)/ (seedMotifCoeff + intercept)) * scaler;

                    WeightMatrixFeature feature = new WeightMatrixFeature(
                        label,
                        featurePwm,
                        motifsToMandatoriesMap.get(featurePwm)[0],
                        isGoodFeature,
                        false,
                        relAffinity,
                        motifsToModsMap.get(featurePwm),
                        DNATools.complementTable());

                    //System.out.println("Hi 1");
                    featuresList.add(feature);

                    // also add feature revcomp if motif is symmetric
                    if((motifsWithoutCoeffs.isEmpty()) && (seedMotifRevCompSim == 1.0)) {

                        WeightMatrix revCompPWM = WeightMatrixTools.reverseComplement((WeightMatrix)featurePwm, DNATools.complementTable());
                        revCompPWM.setName(new String(featurePwm.getName()+" revComp"));

                        feature = new WeightMatrixFeature(
                            label,
                            revCompPWM,
                            featurePwm,
                            motifsToMandatoriesMap.get(featurePwm)[1],
                            motifsToMandatoriesMap.get(featurePwm)[0],
                            isGoodFeature,
                            true,
                            relAffinity,
                            SymbolListTools.reverseComplement(motifsToModsMap.get(featurePwm), DNATools.complementTable()),
                            DNATools.complementTable());

                        //System.out.println("Hi 1");
                        featuresList.add(feature);

                    }
                }
            }

            // Go through the motifs that DON'T have a coefficient and find it
            for (Object motif : motifsWithoutCoeffs) {

                // the motif is a WeightMatrix
                if (motif instanceof WeightMatrix) {
                    WeightMatrix featurePwm = (WeightMatrix)motif;

                    WeightMatrix revCompPWM = WeightMatrixTools.reverseComplement((WeightMatrix)featurePwm, DNATools.complementTable());
                    double featureCoeff = motifsWithCoeffsMap.get(revCompPWM);
                    double featureRevCompSim = WeightMatrixTools.similarity(featurePwm, revCompPWM);

                    // don't include the feature if the coefficient is 0
                    if (featureCoeff == 0) {
                        continue;
                    }
                    //                 else {
                    //                     //double[] mandatoriesArray = motifsToMandatoriesMap.get(featurePwm)[0];
                    //                     Symbol[] modsArray = motifsToModsMap.get(featurePwm);
                    //                     try {
                    //                         for (int i=0; i < seedPwm.columns(); i++) {
                    //                             if (modsArray[i*2] != null) {
                    //                                 double w = seedPwm.getColumn(i).getWeight(modsArray[i*2]);
                    //                                 System.out.println("i="+i+"; w="+w);
                    //                                 featureCoeff = featureCoeff / w;
                    //                                 System.out.println("modifiedFeatureCoeff="+featureCoeff);
                    //                             }
                    //                         }
                    //                     }
                    //                     catch (Exception ex) {
                    //                         ex.printStackTrace();
                    //                     }
                    //                 }

                    String name = featurePwm.getName();
                    String label;
                    if (name.contains("insert")) {
                        label = "insertion";
                    }
                    else if (name.contains("delete")) {
                        label = "deletion";
                    }
                    else {
                        label = "dependency";
                    }

                    //double relAffinity = Math.abs(featureCoeff / seedMotifCoeff) * scaler;
                    //double relAffinity = (featureCoeff / seedMotifCoeff) * scaler;

                    boolean isGoodFeature;
                    double relAffinity = 0;

                    if (additive) {
                        //relAffinity = Math.abs(featureCoeff / seedMotifCoeff) * scaler;
                        relAffinity = (featureCoeff / seedMotifCoeff) * scaler;

                        if (Math.signum(seedMotifCoeff) == Math.signum(featureCoeff)) {
                            isGoodFeature = true;
                        }
                        else {
                            isGoodFeature = false;
                        }
                    }
                    else {
                        relAffinity = ((featureCoeff + seedMotifCoeff) / seedMotifCoeff) * scaler;
                        //relAffinity = ( ((featureCoeff/(featureRevCompSim+1)) + (seedMotifCoeff/(seedMotifRevCompSim+1))) / (seedMotifCoeff/(seedMotifRevCompSim+1)) ) * scaler;
                        if (relAffinity < 0) {
                            relAffinity = 0;
                        }

                        if (relAffinity >= 1.0) {
                            isGoodFeature = true;
                        }
                        else {
                            isGoodFeature = false;
                        }
                    }

                    //double intercept = coefficients[0];
                    //double relAffinity = Math.abs((featureCoeff + intercept)/ (seedMotifCoeff + intercept)) * scaler;

                    //                 if (relAffinity < 1.5) {
                    //                     continue;
                    //                 }

//                     System.out.println("featureCoeff="+featureCoeff);
//                     System.out.println("featureRevCompSim="+featureRevCompSim);
//                     System.out.println("seedMotifCoeff="+seedMotifCoeff);
//                     System.out.println("seedMotifRevCompSim="+seedMotifRevCompSim);
//                     System.out.println("relAffinity="+relAffinity);

                    WeightMatrixFeature feature = new WeightMatrixFeature(
                        label,
                        featurePwm,
                        revCompPWM,
                        motifsToMandatoriesMap.get(featurePwm)[0],
                        motifsToMandatoriesMap.get(featurePwm)[1],
                        isGoodFeature,
                        true,
                        relAffinity,
                        motifsToModsMap.get(featurePwm),
                        DNATools.complementTable());

                    //System.out.println("Hi 2");
                    featuresList.add(feature);
                }
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

        return(featuresList);
    }


//     public static RegressionFitData average(ArrayList<RegressionFitData> fitDataArrayList) {
//         BigDecimal[] coefficientsAvg = new double[fitDataArrayList.get(0).coefficients.length];
//         double rSquaredSum = 0;
//         double stddevSum = 0;
//         double rssSum = 0;
//         LinkedHashMap<Object, BigDecimal> motifsToCoeffsAvgMap = new LinkedHashMap<Object, Double>();
//         int numRegressions = fitDataArrayList.size();

//         // Create all the sums over all the model fits in the ArrayList
//         for (int i=0; i < numRegressions; i++) {

//             RegressionFitData aRegressionFitData = fitDataArrayList.get(i);

//             if (aRegressionFitData.coefficients.length != coefficientsSum.length) {
//                 System.out.println("\nRegressionFitData.average(): Coefficients Array of wrong length!");
//                 System.exit(-1);
//             }

//             for (int j=0; j < coefficientsSum.length; j++) {
//                 coefficientsSum[j] += aRegressionFitData.coefficients[j];
//             }
//             rSquaredSum += aRegressionFitData.rSquared;
//             stddevSum += aRegressionFitData.stddev;
//             rssSum += aRegressionFitData.rss;

//             for(Object motif : aRegressionFitData.motifsWithCoeffsMap.keySet()) {
//                 Double currentSumObject = motifsToCoeffsSumMap.get(motif);
//                 double currentSum = 0;
//                 if (currentSumObject != null) {
//                     currentSum = currentSumObject.doubleValue();
//                 }
//                 motifsToCoeffsSumMap.put(motif, currentSum + aRegressionFitData.motifsWithCoeffsMap.get(motif));
//             }
//         }

//         // Normalize the sums by the size
//         for (int j=0; j < coefficientsSum.length; j++) {
//             coefficientsAvg[j] += coefficientsSum[j] / numRegressions;
//         }

//         for(Object motif : motifsToCoeffsSumMap.keySet()) {
//             double sum = motifsToCoeffsSumMap.get(motif);
//             motifsToCoeffsAvgMap.put(motif, sum / numRegressions);
//         }

//         // create the averaged RegressionFitData Object
//         RegressionFitData avgRegressionFitData = new RegressionFitData(
//             motifsToCoeffsAvgMap,
//             fitDataArrayList.get(0).motifsWithoutCoeffs,
//             fitDataArrayList.get(0).motifsToModsMap,
//             coefficientsAvg,
//             rSquaredSum / numRegressions,
//             stddevSum / numRegressions,
//             rssSum / numRegressions);

//         return(avgRegressionFitData);
//     }

    public void averageRevComps() {
        try {
            // doesn't touch the coefficients[] array!

            // average revComps found in the motifsWithCoeffsMap
            for(Object motif : motifsWithCoeffsMap.keySet()) {
                Object revCompMotif = null;

                // get the revCompMotif
                if (motif instanceof WeightMatrix) {
                    revCompMotif = WeightMatrixTools.reverseComplement((WeightMatrix)motif, DNATools.complementTable());
                }
                else if (motif instanceof FeaturedWeightMatrix) {
                    continue;
                }
                else { // it's a SymbolList
                    revCompMotif = DNATools.reverseComplement((SymbolList)motif);
                }

                // lookup the coeff for the revCompMotif
                Double revCompCoeff = motifsWithCoeffsMap.get(revCompMotif);
                if (revCompCoeff != null) {
                    double motifCoeff = motifsWithCoeffsMap.get(motif);
                    double average = (motifCoeff + revCompCoeff.doubleValue()) / 2;
                    motifsWithCoeffsMap.put(motif, average);
                    motifsWithCoeffsMap.put(revCompMotif, average);
                }
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public static RegressionFitData average(ArrayList<RegressionFitData> fitDataArrayList) {
        double[] coefficientsSum = new double[fitDataArrayList.get(0).coefficients.length];
        double[] coefficientsAvg = new double[fitDataArrayList.get(0).coefficients.length];

        double[] residualsSum = new double[fitDataArrayList.get(0).residuals.length];
        double[] residualsAvg = new double[fitDataArrayList.get(0).residuals.length];

        double alphaSum = 0;
        double betaSum = 0;
        double nonSpecKaSum = 0;

        double rSquaredSum = 0;
        double stddevSum = 0;
        double rssSum = 0;

        LinkedHashMap<Object, Double> motifsToCoeffsSumMap = new LinkedHashMap<Object, Double>();
        LinkedHashMap<Object, Double> motifsToCoeffsAvgMap = new LinkedHashMap<Object, Double>();
        int numRegressions = fitDataArrayList.size();

        // Create all the sums over all the model fits in the ArrayList
        for (int i=0; i < numRegressions; i++) {

            RegressionFitData aRegressionFitData = fitDataArrayList.get(i);

            if (aRegressionFitData.coefficients.length != coefficientsSum.length) {
                System.out.println("\nRegressionFitData.average(): Coefficients Array of wrong length!");
                System.exit(-1);
            }

            for (int j=0; j < coefficientsSum.length; j++) {
                coefficientsSum[j] += aRegressionFitData.coefficients[j];
            }

            for (int j=0; j < residualsSum.length; j++) {
                residualsSum[j] += aRegressionFitData.residuals[j];
            }

            alphaSum += aRegressionFitData.alpha;
            betaSum += aRegressionFitData.beta;
            nonSpecKaSum += aRegressionFitData.nonSpecKa;

            rSquaredSum += aRegressionFitData.rSquared;
            stddevSum += aRegressionFitData.stddev;
            rssSum += aRegressionFitData.rss;

            for(Object motif : aRegressionFitData.motifsWithCoeffsMap.keySet()) {
                Double currentSumObject = motifsToCoeffsSumMap.get(motif);
                double currentSum = 0;
                if (currentSumObject != null) {
                    currentSum = currentSumObject.doubleValue();
                }
                motifsToCoeffsSumMap.put(motif, currentSum + aRegressionFitData.motifsWithCoeffsMap.get(motif));
            }

        }

        // Normalize the sums by the size
        for (int j=0; j < coefficientsSum.length; j++) {
            coefficientsAvg[j] += coefficientsSum[j] / numRegressions;
        }

        // Normalize the sums by the size
        for (int j=0; j < residualsSum.length; j++) {
            residualsAvg[j] += residualsSum[j] / numRegressions;
        }

        // Normalize the sums by the size
        for(Object motif : motifsToCoeffsSumMap.keySet()) {
            double sum = motifsToCoeffsSumMap.get(motif);
            motifsToCoeffsAvgMap.put(motif, sum / numRegressions);
        }

        // create the averaged RegressionFitData Object
        RegressionFitData avgRegressionFitData = new RegressionFitData(
            motifsToCoeffsAvgMap,
            fitDataArrayList.get(0).motifsWithoutCoeffs,
            fitDataArrayList.get(0).motifsToModsMap,
            fitDataArrayList.get(0).motifsToMandatoriesMap,
            coefficientsAvg,
            residualsAvg,
            alphaSum / numRegressions,
            betaSum / numRegressions,
            nonSpecKaSum / numRegressions,
            fitDataArrayList.get(0).startPositions,
            rSquaredSum / numRegressions,
            stddevSum / numRegressions,
            rssSum / numRegressions);

        return(avgRegressionFitData);
    }

    public String toString() {
        StringBuffer outBuffer = new StringBuffer();

        try {
            outBuffer.append("\nThe coefficients are :");
            outBuffer.append("\n\t" + "Intercept"+ ": " + getIntercept() +"");

            if (startPositions == null) {
            //if (true) {
                // Go through the motifs that have a coefficient
                for (Object motif : motifsWithCoeffsMap.keySet()) {
                    String independentVar;
                    if (motif instanceof WeightMatrix) {
                        independentVar = ((WeightMatrix)motif).getName();
                    }
                    else if (motif instanceof FeaturedWeightMatrix) {
                        independentVar = ((FeaturedWeightMatrix)motif).getName();
                    }
                    else { // it's a SymbolList
                        independentVar = ((SymbolList)motif).seqString();
                    }
                    outBuffer.append("\n\t" + independentVar +": " + motifsWithCoeffsMap.get(motif) +"");
                }
            }
            else {
                WeightMatrixTools.BindingStrand[] strandsArray = WeightMatrixTools.BindingStrand.values();

                // look at both strands as default
                int maxStrandIndex = 2;

                // if (coefficients < 2*startPositions.len) then use only the positive strand
                if (coefficients.length < 2*startPositions.length) {
                    maxStrandIndex = 1;
                }

                // loop for the pos and neg strand
                for (int strandIndex = 0; strandIndex < maxStrandIndex; strandIndex++) {
                    WeightMatrixTools.BindingStrand currentStrand = strandsArray[strandIndex];
                    // loop for each startPosition
                    for (int posIndex=0; posIndex < startPositions.length; posIndex++) {
                        outBuffer.append("\n\tmotif(strand"+strandIndex+", "+startPositions[posIndex]+")"+": " + coefficients[posIndex+1 + (strandIndex*startPositions.length)]+"");
                    }
                }
            }

            // Go through the motifs that DON'T have a coefficient and find it
            if ((motifsWithoutCoeffs != null) && (!motifsWithoutCoeffs.isEmpty())) {
                outBuffer.append("\n\nThe additional coefficient mappings are :");
                for (Object motif : motifsWithoutCoeffs) {
                    String independentVar;
                    Double coeff;
                    if (motif instanceof WeightMatrix) {
                        independentVar = ((WeightMatrix)motif).getName();
                        coeff = motifsWithCoeffsMap.get(WeightMatrixTools.reverseComplement((WeightMatrix)motif, DNATools.complementTable()));
                    }
                    else if (motif instanceof FeaturedWeightMatrix) {
                        independentVar = ((FeaturedWeightMatrix)motif).getName();
                        coeff = 0.0;
                    }
                    else {
                        independentVar = ((SymbolList)motif).seqString();
                        coeff = motifsWithCoeffsMap.get(DNATools.reverseComplement((SymbolList)motif));
                    }
                    outBuffer.append("\n\t" + independentVar +": " + coeff +"");
                }
            }

            outBuffer.append("\n");
            if (!Double.isNaN(alpha)) {
                outBuffer.append("\nThe alpha is "+alpha+".");
            }
            if (!Double.isNaN(beta)) {
                outBuffer.append("\nThe beta is "+beta+".");
            }
            if (!Double.isNaN(nonSpecKa)) {
                outBuffer.append("\nThe nonSpecKa is "+nonSpecKa+".");
            }

            outBuffer.append("\n");
            if (!Double.isNaN(rSquared)) {
                outBuffer.append("\nThe R-squared is "+rSquared+".");
            }
            if (!Double.isNaN(stddev)) {
                outBuffer.append("\nThe stddev is "+stddev+".");
            }
            if (!Double.isNaN(rss)) {
                outBuffer.append("\nThe RSS is "+rss+".");
            }

            outBuffer.append("\n");
        }

        catch (Exception ex) {
            ex.printStackTrace();
        }

        return(outBuffer.toString());
    }
}
