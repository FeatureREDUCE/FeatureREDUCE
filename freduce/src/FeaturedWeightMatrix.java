
import java.io.Serializable;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SymbolListTools;
import org.biojava.bio.symbol.ReversibleTranslationTable;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.mydp.*;
import org.biojava.bio.dp.*;
import org.biojava.bio.dist.*;
import java.util.*;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.RNATools;

public class FeaturedWeightMatrix  implements Serializable {

    //private static final long serialVersionUID = 7923870404512774813L;
    //private static final long serialVersionUID = 7923870404512774814L; // 12/22/2011
    private static final long serialVersionUID = 7923870404512774815L;


    String name;
    WeightMatrix posStrandPWM;
    WeightMatrix negStrandPWM;

    WeightMatrix oldAMposStrandPWM = null; // additive model PSAM!!
    WeightMatrix oldAMnegStrandPWM = null; // additive model PSAM!!

    boolean[][] mandatoryColumns; // [strand][pos] for the posStrandPWM and negStrandPWM
    boolean isRevCompPalindrome;
    WeightMatrixTools.BindingStrand strand;
    WeightMatrixTools.BothStrandsCalc calc;

    Table featuresTable = null; //temporary sorted table above threshold
    LinkedHashMap<FeatureKey, WeightMatrixFeature> featureKeyToAddFeatureMap = null;
    LinkedHashMap<FeatureKey, WeightMatrixFeature> featureKeyToMultFeatureMap = null;

    boolean sorted = false;
    double[][] positionalWeights = null; // [strand][pos]

    double intercept = 0;
    double eToMu = 0;
    double nonSpecKa = 0;
    double revCompSimilarity;
    boolean includePwm = true;
    int columns = 0;
    Alphabet alphabet;
    boolean additive = false;

    double intensitiesScaler = 1;

    /////////////////////////////////////////////////
    // all-kmer models
    /////////////////////////////////////////////////
    char[] kmerLengths; // [kmerIndex]
    double[][][] kmerPositionalWeights = null; // [kmerIndex][strand][pos]
    int[][] kmerStartPositions; // [kmerIndex][pos]
    double[][] kmerToAffinityMatrix = null; // [kmerIndex][kmer]
    char[][] revCompMatrix = null; // [kmerIndex][kmer]
    boolean[] includeRevComps = null; // [kmerIndex]
    //TObjectCharHashMap symListToWordMap = null;


    //    public FeaturedWeightMatrix(String aLabel, boolean aGoodFlag, double aRelativeAffinity, double anR2, String[] aSymbolArray, boolean[] aDeleteArray) {
    public FeaturedWeightMatrix(
        String aName,
        WeightMatrix aPosStrandPWM)
        //        boolean additive)

    {
        WeightMatrixTools.BindingStrand aStrand = null;
        WeightMatrixTools.BothStrandsCalc aCalc = null;
        ReversibleTranslationTable complementTable = null;

        if (aPosStrandPWM.getAlphabet() == DNATools.getDNA()) {
            aStrand = WeightMatrixTools.BindingStrand.BOTH;
            aCalc = WeightMatrixTools.BothStrandsCalc.NORMED_SUM;
            complementTable = DNATools.complementTable();
        }
        else { //RNA
            aStrand = WeightMatrixTools.BindingStrand.BOTH;
            complementTable = RNATools.complementTable();
        }

        setName(aName);
        putPWMs(aPosStrandPWM, complementTable);
        setStrand(aStrand);
        setCalc(aCalc);
    }

    public FeaturedWeightMatrix(
        String aName,
        WeightMatrix aPosStrandPWM,
        WeightMatrixTools.BindingStrand aStrand,
        WeightMatrixTools.BothStrandsCalc aCalc,
        ReversibleTranslationTable complementTable)
        //        boolean additive)

    {
        setName(aName);
        putPWMs(aPosStrandPWM, complementTable);
        setStrand(aStrand);
        setCalc(aCalc);
        //setBaseModel(additive);
        //includePwm = true;
    }

    public FeaturedWeightMatrix(
        String aName,
        WeightMatrix aPosStrandPWM,
        WeightMatrixTools.BindingStrand aStrand,
        WeightMatrixTools.BothStrandsCalc aCalc,
        ReversibleTranslationTable complementTable,
        boolean anIncludePwm,
        ArrayList<WeightMatrixFeature> aFeaturesList,
        boolean additive)
    {
        //this(aName, aPosStrandPWM, aStrand, aCalc, complementTable, additive);
        this(aName, aPosStrandPWM, aStrand, aCalc, complementTable);
        if (aFeaturesList != null) {
            setFeatures(anIncludePwm, aFeaturesList, additive);
        }
    }

    public FeaturedWeightMatrix(
        String aName,
        WeightMatrix aPosStrandPWM,
        WeightMatrixTools.BindingStrand aStrand,
        WeightMatrixTools.BothStrandsCalc aCalc,
        ReversibleTranslationTable complementTable,
        boolean anIncludePwm,
        ArrayList<WeightMatrixFeature> additiveFeaturesList,
        ArrayList<WeightMatrixFeature> multiplicativeFeaturesList)
    {
        //this(aName, aPosStrandPWM, aStrand, aCalc, complementTable, additive);
        this(aName, aPosStrandPWM, aStrand, aCalc, complementTable);
        if (additiveFeaturesList != null) {
            setFeatures(anIncludePwm, additiveFeaturesList, true);
        }
        if (multiplicativeFeaturesList != null) {
            setFeatures(anIncludePwm, multiplicativeFeaturesList, false);
        }
    }

    // Objects are not copied!! So the clone will contain same objects as source!
    public FeaturedWeightMatrix clone() {
        return(clone(0.0, null));
    }

    public FeaturedWeightMatrix clone(double threshold, boolean[] excludePositions) {
        ArrayList<WeightMatrixFeature> addFeaturesList = getFeatures(threshold, true, excludePositions);
        ArrayList<WeightMatrixFeature> multFeaturesList = getFeatures(threshold, false, excludePositions);

        return(new FeaturedWeightMatrix(
                this.name,
                this.posStrandPWM,
                this.strand,
                this.calc,
                DNATools.complementTable(),
                this.includePwm,
                addFeaturesList,
                multFeaturesList));
    }

    // This is a convenience method that will generate a WeightMatrix -> cloned(FeaturedWeightMatrix)) map
    // where all the cloned objects have the same features removed
    // and their PSAMs set according to the keys on the aMotifToHammingSphereMap
    public LinkedHashMap<Object, FeaturedWeightMatrix> clones(
        LinkedHashMap<Object, LinkedHashSet<SymbolList>> aMotifToHammingSphereMap,
        LinkedHashMap<Object, boolean[][]> aMotifsToMandatoriesMap,
        ArrayList<FeatureKey> aFeatureRemovalList)
    {
        System.out.println("\nRemoving "+aFeatureRemovalList.size()+" features for each FSAM clone.");

        LinkedHashMap<Object, FeaturedWeightMatrix> motifToFwmMap = new LinkedHashMap<Object, FeaturedWeightMatrix>();

        for (Object motif : aMotifToHammingSphereMap.keySet()) {
            WeightMatrix aPwm = (WeightMatrix)motif;

            FeaturedWeightMatrix clone = this.clone();
            clone.putPWMs(aPwm, DNATools.complementTable());
            clone.setMandatoryColumns(aMotifsToMandatoriesMap.get(motif));
            clone.removeFeatures(aFeatureRemovalList);

            motifToFwmMap.put(motif, clone);

        }
        return(motifToFwmMap);
    }

    public void reHashFeatures() {

        if (featureKeyToMultFeatureMap == null) {
            return;
        }

        LinkedHashMap<FeatureKey, WeightMatrixFeature> newFeatureKeyToFeatureMap = new LinkedHashMap<FeatureKey, WeightMatrixFeature>();
        for (WeightMatrixFeature feature : featureKeyToMultFeatureMap.values()) {
            FeatureKey newFeatureKey = new FeatureKey(feature.getModsArray());
            feature.setKey(newFeatureKey);
            newFeatureKeyToFeatureMap.put(newFeatureKey, feature);
        }
        featureKeyToMultFeatureMap = newFeatureKeyToFeatureMap;


        newFeatureKeyToFeatureMap = new LinkedHashMap<FeatureKey, WeightMatrixFeature>();
        for (WeightMatrixFeature feature : featureKeyToAddFeatureMap.values()) {
            FeatureKey newFeatureKey = new FeatureKey(feature.getModsArray());
            feature.setKey(newFeatureKey);
            newFeatureKeyToFeatureMap.put(newFeatureKey, feature);
        }
        featureKeyToAddFeatureMap = newFeatureKeyToFeatureMap;

    }

    public String toString() {
        return(toString(false));
    }

    public String toString(boolean additive) {
        StringBuffer stringBuffer = new StringBuffer();
        stringBuffer.append("\nname = "+name);
        stringBuffer.append("\ncolumns = "+columns);
        stringBuffer.append("\nisRevCompPalindrome = "+isRevCompPalindrome);
        stringBuffer.append("\nrevCompSimilarity = "+revCompSimilarity);
        stringBuffer.append("\nstrand = "+strand);
        stringBuffer.append("\ncalc = "+calc);
        stringBuffer.append("\nincludePwm = "+includePwm);
        stringBuffer.append("\neToMu = "+eToMu);
        stringBuffer.append("\nnonSpecKa = "+nonSpecKa);

        ///////////////////////////////////////////////////
        // always display the multiplicative features
        // additive = false;
        ///////////////////////////////////////////////////


        if (posStrandPWM != null) {
            if (additive && (oldAMposStrandPWM != null)) {
                stringBuffer.append("\n\nAMposStrandPWM:\n"+WeightMatrixTools.toString(oldAMposStrandPWM));
            }
            else if (posStrandPWM != null) {
                stringBuffer.append("\n\nposStrandPWM:\n"+WeightMatrixTools.toString(posStrandPWM));
            }
        }

        if (positionalWeights != null) {
            String[] rowLabels = {"posStrandRelBias", "negStrandRelBias"};
            int trainingProbeSeqLengths = this.positionalWeights[0].length - this.columns + 1;
            int[] startPositions = getStartPositions(trainingProbeSeqLengths);
            String[] columnLabels = StringTools.toStrings(startPositions);
            String outputTable = StringTools.toString(positionalWeights, rowLabels, columnLabels);
            stringBuffer.append("\n\n"+outputTable+"\n");
        }

        // get the model features by additive flag
        getFeatures(0.0, additive);

        if (featuresTable != null) {
            // sort the table by the relative affinities if the table is in an unsorted state
            //sortIfNeeded();
            if (additive) {
                stringBuffer.append("\n\nAdditive Features:");
            }
            else {
                stringBuffer.append("\n\nMultiplicative Features:");
            }
            int count=0;
            for (Iterator iter = featuresTable.iterator(); iter.hasNext(); ) {
                count++;
                List rowList = (List)iter.next();
                WeightMatrixFeature feature = (WeightMatrixFeature)rowList.get(0);
                stringBuffer.append("\n\n("+count+") "+feature.toString());
            }
        }
        return(stringBuffer.toString());
    }

    public int[] getStartPositions(int probeSeqLengths) {
        int numKmerWindows = probeSeqLengths + this.columns - 1;
        int[] startPositions = new int[numKmerWindows];

        int firstStartPos = -1 * (this.columns - 1);
        for (int posIndex= 0; posIndex < numKmerWindows; posIndex++) {
            startPositions[posIndex] = firstStartPos + posIndex;
        }
        return(startPositions);
    }


    public WeightMatrixFeature getFeature(FeatureKey aFeatureKey) {
        return(getFeature(aFeatureKey, false));
    }

    public WeightMatrixFeature getFeature(FeatureKey aFeatureKey, boolean additive) {
        if (additive) {
            if (featureKeyToAddFeatureMap != null) {
                return(featureKeyToAddFeatureMap.get(aFeatureKey));
            }
        }
        else {
            if (featureKeyToMultFeatureMap != null) {
                return(featureKeyToMultFeatureMap.get(aFeatureKey));
            }
        }
        return(null);
    }


    public void removeFeatures() {
        featuresTable = null;
        featureKeyToMultFeatureMap = null;
        featureKeyToAddFeatureMap = null;

        sorted = true;
        includePwm = true;
    }

    public void removeFeature(WeightMatrixFeature aPwmFeature, boolean additive) {
        if (additive) {
            if (featureKeyToAddFeatureMap != null) {
                featureKeyToAddFeatureMap.remove(aPwmFeature.getKey());
                sorted = false;
            }
        }
        else {
            if (featureKeyToMultFeatureMap != null) {
                featureKeyToMultFeatureMap.remove(aPwmFeature.getKey());
                sorted = false;
            }
        }
    }

//     public void removeFeatures(ArrayList<WeightMatrixFeature> aFeaturesList) {
//         if (featureKeyToFeatureMap != null) {
//             // populate the featuresMap
//             for (WeightMatrixFeature feature : aFeaturesList) {
//                 featureKeyToFeatureMap.remove(feature.getKey());
//             }
//             sorted = false;
//         }
//     }

    public void removeFeature(FeatureKey aFeatureKey, boolean additive) {
        if (additive) {
            if (featureKeyToAddFeatureMap != null) {
                featureKeyToAddFeatureMap.remove(aFeatureKey);
                sorted = false;
            }
        }
        else {
            if (featureKeyToMultFeatureMap != null) {
                featureKeyToMultFeatureMap.remove(aFeatureKey);
                sorted = false;
            }
        }
    }

    public void removeFeatures(ArrayList<FeatureKey> aFeatureKeysList, boolean additive) {
        if (additive) {
            if (featureKeyToAddFeatureMap != null) {
                // traverse the featuresMap
                for (FeatureKey featureKey : aFeatureKeysList) {
                    featureKeyToAddFeatureMap.remove(featureKey);
                }
                sorted = false;
            }
        }
        else {
            if (featureKeyToMultFeatureMap != null) {
                // traverse the featuresMap
                for (FeatureKey featureKey : aFeatureKeysList) {
                    featureKeyToMultFeatureMap.remove(featureKey);
                }
                sorted = false;
            }
        }
    }

    public void removeFeatures(ArrayList<FeatureKey> aFeatureKeysList) {
        if (featureKeyToAddFeatureMap != null) {
            // traverse the featuresMap
            for (FeatureKey featureKey : aFeatureKeysList) {
                featureKeyToAddFeatureMap.remove(featureKey);
            }
            sorted = false;
        }

        if (featureKeyToMultFeatureMap != null) {
            // traverse the featuresMap
            for (FeatureKey featureKey : aFeatureKeysList) {
                featureKeyToMultFeatureMap.remove(featureKey);
            }
            sorted = false;
        }
    }

    public boolean hasHighOrderFeatures() {
        if (featureKeyToMultFeatureMap != null) {
            return(true);
        }
        return(false);
    }

    public Alphabet getAlphabet() {
        return(alphabet);
    }

    public char[] getKmerLengths() {
        return(kmerLengths);
    }

    public boolean[] getRevCompEquivFlags() {
        return(includeRevComps);
    }

    public double [] getKmerToAffinityMatrix(int index) {
        return(kmerToAffinityMatrix[index]);
    }

    public double [][] getKmerToAffinityMatrix() {
        return(kmerToAffinityMatrix);
    }

    public double[][] getKmerPositionalWeights(int index) {
        return(kmerPositionalWeights[index]);
    }

    public WeightMatrix getPosStrandPWM() {
        return(posStrandPWM);
    }

    public WeightMatrix getNegStrandPWM() {
        return(negStrandPWM);
    }

    public WeightMatrix getPWM(WeightMatrixTools.BindingStrand strand) {
        if (strand == WeightMatrixTools.BindingStrand.POS) {
            return(posStrandPWM);
        }
        else {
            return(negStrandPWM);
        }
    }

    public double getRevCompSimilarity() {
        return(revCompSimilarity);
    }

    public void setRevCompSimilarity(double aRevCompSimilarity) {
        revCompSimilarity = aRevCompSimilarity;
    }

    public void setBaseModel(boolean anAdditive) {
        additive = anAdditive;
    }

    public boolean[][] getMandatoryColumns() {
        return(mandatoryColumns);
    }

    public void setMandatoryColumns(boolean[][] aMandatoryColumns) {
        mandatoryColumns = aMandatoryColumns;
    }

    public double getIntercept() {
        return(intercept);
    }

    public void setIntercept(double anIntercept) {
        intercept = anIntercept;
    }

    public double getScaler() {
        return(intensitiesScaler);
    }

    public void setScaler(double aScaler) {
        intensitiesScaler = aScaler;
    }

    public double getEToMu() {
        return(eToMu);
    }

    public void setEToMu(double anEToMu) {
        eToMu = anEToMu;
    }

    public double getNonSpecKa() {
        return(nonSpecKa);
    }

    public void setNonSpecKa(double aNonSpecKa) {
        nonSpecKa = aNonSpecKa;
    }

    public double[][] getPositionalWeights() {
        return(positionalWeights);
    }

    public void setPositionalWeights(double[][] aPositionalWeights) {
        positionalWeights = aPositionalWeights;
    }

    public String getName() {
        return(this.name);
    }

    public int columns() {
        return(this.columns);
    }

    public void setName(String aName) {
        this.name = aName;
        return;
    }

    public WeightMatrixTools.BindingStrand getStrand() {
        return(strand);
    }

    public WeightMatrixTools.BothStrandsCalc getCalc() {
        return(calc);
    }

    public void putPWMs(WeightMatrix aPosStrandPWM, ReversibleTranslationTable complementTable) {
        try {
            WeightMatrix aNegStrandPWM = WeightMatrixTools.reverseComplement(aPosStrandPWM, complementTable);
            putPWMs(aPosStrandPWM, aNegStrandPWM);
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public void putPWMs(WeightMatrix aPosStrandPWM, WeightMatrix aNegStrandPWM) {
        try {
            posStrandPWM = aPosStrandPWM;
            negStrandPWM = aNegStrandPWM;
            isRevCompPalindrome = WeightMatrixTools.areEmissionSpectraEqual(posStrandPWM, negStrandPWM);
            revCompSimilarity = WeightMatrixTools.similarity(posStrandPWM, negStrandPWM);
            columns = aPosStrandPWM.columns();
            alphabet = posStrandPWM.getAlphabet();
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public void setStrand(WeightMatrixTools.BindingStrand aStrand) {
        strand = aStrand;
    }

    public void setCalc(WeightMatrixTools.BothStrandsCalc aCalc) {
        calc = aCalc;
    }

    public void setFeatures(boolean anIncludePwm, ArrayList<WeightMatrixFeature> aFeaturesList, boolean additive) {
        //featuresTable = new Table();
        if (additive) {
            featureKeyToAddFeatureMap = new LinkedHashMap<FeatureKey, WeightMatrixFeature>();
        }
        else {
            featureKeyToMultFeatureMap = new LinkedHashMap<FeatureKey, WeightMatrixFeature>();
        }
        addFeatures(anIncludePwm, aFeaturesList, additive);
        sorted = false;
    }

    public void addFeatures(boolean anIncludePwm, ArrayList<WeightMatrixFeature> aFeaturesList) {
        addFeatures(anIncludePwm, aFeaturesList, false);
    }

    public void addFeatures(boolean anIncludePwm, ArrayList<WeightMatrixFeature> aFeaturesList, boolean anAdditive) {
        includePwm = anIncludePwm;
        this.additive = anAdditive;

        if (anAdditive) {
            //if (featuresTable == null) {
            if (featureKeyToAddFeatureMap == null) {
                //featuresTable = new Table();
                featureKeyToAddFeatureMap = new LinkedHashMap<FeatureKey, WeightMatrixFeature>();
            }

            // populate the featuresMap
            for (WeightMatrixFeature feature : aFeaturesList) {
                featureKeyToAddFeatureMap.put(feature.getKey(), feature);
            }
            sorted = false;
        }
        else {
            //if (featuresTable == null) {
            if (featureKeyToMultFeatureMap == null) {
                //featuresTable = new Table();
                featureKeyToMultFeatureMap = new LinkedHashMap<FeatureKey, WeightMatrixFeature>();
            }

            // populate the featuresMap
            for (WeightMatrixFeature feature : aFeaturesList) {
                featureKeyToMultFeatureMap.put(feature.getKey(), feature);
            }
            sorted = false;
        }
    }

    public void avgFeatures(
        boolean anIncludePwm,
        ArrayList<WeightMatrixFeature> aFeaturesList,
        boolean anAdditive,
        double newFeatureWeight)
    {

        LinkedHashMap<FeatureKey, WeightMatrixFeature> featureKeyToFeatureMap = null;
        if (anAdditive) {
            featureKeyToFeatureMap = featureKeyToAddFeatureMap;
        }
        else {
            featureKeyToFeatureMap = featureKeyToMultFeatureMap;
        }

        if (featureKeyToFeatureMap == null) {
            return;
        }

        for (WeightMatrixFeature newFeature : aFeaturesList) {
            WeightMatrixFeature oldFeature = featureKeyToFeatureMap.get(newFeature.getKey());

            if (oldFeature == null) {

                System.out.println("\nMultiplicative NNDD feature "+oldFeature.getPosStrandPWM().getName()+" not found! Creating a new one.");

                // scale down the multiplicative estimate so that we don't over-estimate
                newFeature.setRelativeAffinity((1 + newFeature.getRelativeAffinity()) *.5); // go half-way

                // put the new feature into the map
                featureKeyToFeatureMap.put(newFeature.getKey(), newFeature);
                sorted = false;
            }
            else {
                // weighted average of the relative affinities
                double newRelAffinity =
                    (newFeatureWeight * newFeature.getRelativeAffinity())
                    + ((1-newFeatureWeight) * oldFeature.getRelativeAffinity());

                oldFeature.setRelativeAffinity(newRelAffinity);

                if (newRelAffinity >= 1.0) {
                    oldFeature.setIsGood(true);
                }
                else {
                    oldFeature.setIsGood(false);
                }
            }

        }
    }

    public void calcAvgRelativeAffinities() {
        if (featureKeyToMultFeatureMap != null) {
            for(WeightMatrixFeature feature : featureKeyToMultFeatureMap.values()) {
                ArrayList<Double> modsRelAffinities = new ArrayList<Double>();
                int modsStartPos = feature.getModsStartPos();
                Symbol[] tempModsArray = null;

                if (modsStartPos == 0) {
                    for (Iterator dnaIter1 = ((FiniteAlphabet)alphabet).iterator(); dnaIter1.hasNext(); ) {
                        Symbol mut1Symbol = (Symbol)dnaIter1.next();
                        double relAffinitySum = 0;

                        // Left
                        relAffinitySum += feature.getRelativeAffinity();

                        // Right
                        tempModsArray = new Symbol[(this.columns*2)-1];
                        tempModsArray[modsStartPos+2] = feature.getModsArray()[modsStartPos+2];
                        tempModsArray[modsStartPos+4] = mut1Symbol;
                        FeatureKey featureKey =  new FeatureKey(tempModsArray);
                        WeightMatrixFeature neighborFeature = (WeightMatrixFeature)featureKeyToMultFeatureMap.get(featureKey);
                        if (neighborFeature != null) {
                            relAffinitySum += neighborFeature.getRelativeAffinity();
                        }

                        modsRelAffinities.add(relAffinitySum);
                    }

                }

                else if (modsStartPos == (feature.getModsArray().length-3)) {
                    for (Iterator dnaIter1 = ((FiniteAlphabet)alphabet).iterator(); dnaIter1.hasNext(); ) {
                        Symbol mut1Symbol = (Symbol)dnaIter1.next();
                        double relAffinitySum = 0;

                        // Left
                        tempModsArray = new Symbol[(this.columns*2)-1];
                        tempModsArray[modsStartPos-2] = mut1Symbol;
                        tempModsArray[modsStartPos] = feature.getModsArray()[modsStartPos];
                        FeatureKey featureKey = new FeatureKey(tempModsArray);
                        WeightMatrixFeature neighborFeature = (WeightMatrixFeature)featureKeyToMultFeatureMap.get(featureKey);
                        if (neighborFeature != null) {
                            relAffinitySum += neighborFeature.getRelativeAffinity();
                        }

                        // Right
                        relAffinitySum += feature.getRelativeAffinity();

                        modsRelAffinities.add(relAffinitySum);
                    }

                }
                else {
                    for (Iterator dnaIter1 = ((FiniteAlphabet)alphabet).iterator(); dnaIter1.hasNext(); ) {
                        Symbol mut1Symbol = (Symbol)dnaIter1.next();

                        for (Iterator dnaIter2 = ((FiniteAlphabet)alphabet).iterator(); dnaIter2.hasNext(); ) {
                            Symbol mut2Symbol = (Symbol)dnaIter2.next();
                            double relAffinitySum = 0;

                            // Left
                            tempModsArray = new Symbol[(this.columns*2)-1];
                            tempModsArray[modsStartPos-2] = mut1Symbol;
                            tempModsArray[modsStartPos] = feature.getModsArray()[modsStartPos];
                            FeatureKey featureKey = new FeatureKey(tempModsArray);
                            WeightMatrixFeature neighborFeature = (WeightMatrixFeature)featureKeyToMultFeatureMap.get(featureKey);
                            if (neighborFeature != null) {
                                relAffinitySum += neighborFeature.getRelativeAffinity();
                            }

                            // Center
                            relAffinitySum += feature.getRelativeAffinity();

                            // Right
                            tempModsArray = new Symbol[(this.columns*2)-1];
                            tempModsArray[modsStartPos+2] = feature.getModsArray()[modsStartPos+2];
                            tempModsArray[modsStartPos+4] = mut2Symbol;
                            featureKey = new FeatureKey(tempModsArray);
                            neighborFeature = (WeightMatrixFeature)featureKeyToMultFeatureMap.get(featureKey);
                            if (neighborFeature != null) {
                                relAffinitySum += neighborFeature.getRelativeAffinity();
                            }

                            modsRelAffinities.add(relAffinitySum);
                        }
                    }
                }

                double modsRelAffinitiesAvg = MathTools.mean(ArrayTools.toDoubleArray(modsRelAffinities));
                feature.setAvgRelativeAffinity(modsRelAffinitiesAvg);

            }
        }
    }

    // sort the table by the relative affinities if the table is in an unsorted state
    private void sortIfNeeded() {
        if (!sorted && featuresTable!=null) {
            featuresTable.sort(1); //sort by ascending relativeAffinity values
            //featuresTable.sort(2); //sort by ascending avgRelativeAffinity values
            featuresTable.reverse(); // get descending order
            sorted = true;
        }
    }

    // default is get Mult Features
    public ArrayList<WeightMatrixFeature> getFeatures() {
        return(getFeatures(0.0, false, null));
    }

    // default is get Mult Features
    public ArrayList<WeightMatrixFeature> getFeatures(double threshold) {
        return(getFeatures(threshold, false, null));
    }

    public ArrayList<WeightMatrixFeature> getFeatures(double threshold, boolean additive) {
        return(getFeatures(threshold, additive, null));
    }

    public ArrayList<WeightMatrixFeature> getFeatures(double threshold, boolean additive, boolean[] excludePositions) {

        LinkedHashMap<FeatureKey, WeightMatrixFeature> featureKeyToFeatureMap = null;
        if (additive) {
            featureKeyToFeatureMap = featureKeyToAddFeatureMap;
        }
        else {
            featureKeyToFeatureMap = featureKeyToMultFeatureMap;
        }

        if (featureKeyToFeatureMap == null) {
            return(null);
        }

        // Calculate the average relative affinities for each feature
        calcAvgRelativeAffinities();

        featuresTable = new Table();

        // populate the featureTable
        for (WeightMatrixFeature feature : featureKeyToFeatureMap.values()) {

            Symbol[] posStrandModsArray = feature.getPosStrandModsArray();

            boolean exclude = false;

            if (excludePositions != null) {

                // if the modsArray contains a mod that should not be excluded then keep this feature!!
                exclude = true;
                for (int i=0; i < posStrandModsArray.length; i += 2) { // skip insertions
                    // modsArray includes 1bp insertions
                    if ((posStrandModsArray[i] != null) && (excludePositions[(int)i/2] == false)) {
                        exclude = false;
                        break;
                    }
                }
            }

            if (!exclude) {
                featuresTable.add(Arrays.asList(
                        feature,
                        //                     new Double(Math.abs(feature.getRelativeAffinity())),
                        //                     new Double(Math.abs(feature.getAvgRelativeAffinity()))
                        new Double(Math.abs(1 - feature.getRelativeAffinity())),
                        new Double(Math.abs(1 - feature.getAvgRelativeAffinity()))
                                                ));
            }
        }
        sorted = false;

        // sort the table by the relative affinities if the table is in an unsorted state
        sortIfNeeded();

        ArrayList<WeightMatrixFeature> featuresList = new ArrayList<WeightMatrixFeature>();

        if (featuresTable != null) {
            for (Iterator iter = featuresTable.iterator(); iter.hasNext(); ) {
                List rowList = (List)iter.next();
                WeightMatrixFeature feature = (WeightMatrixFeature)rowList.get(0);
                if (Math.abs(feature.getRelativeAffinity()) < threshold) {
                    break;
                }
                featuresList.add(feature);
            }
        }
        return(featuresList);
    }



    public double getHighestScore() {
        double highestScore = -1;

        try {
            Table allScoresTable = new Table(getAllScores());
            allScoresTable.sort(1, false); // descending
            highestScore = ((Double)allScoresTable.getElement(0,1)).doubleValue();
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

        return(highestScore);
    }

    public LinkedHashMap<SymbolList, Double>  getAllScores() {
        LinkedHashMap<SymbolList, Double> relAffinityLHM = new LinkedHashMap<SymbolList, Double>();

        try {

//             double bestScore = score(DNATools.createDNA("ATGATTAATGAC"), null, featureThreshold);
//             double goodScore = score(DNATools.createDNA("ATGATTAATGGC"), null, featureThreshold);
//             double badScore  = score(DNATools.createDNA("AGCATTAATGGC"), null, featureThreshold);

//             System.out.println("\nbestScore = "+bestScore+"\tgoodScore = "+goodScore+"\tbadScore = "+badScore);

//             highestScore = bestScore;

            //create the cross-product Alphabet
            Alphabet orderNalphabet = AlphabetManager.getCrossProductAlphabet(Collections.nCopies(this.columns, alphabet));

            //use the OrderNDistributionFactory to create the Distribution
            //OrderNDistribution orderNDist = (OrderNDistribution)OrderNDistributionFactory.DEFAULT.createDistribution(orderNalphabet);

            double highestScore = -1;
            for (Iterator dnaIter = ((FiniteAlphabet)orderNalphabet).iterator(); dnaIter.hasNext(); ) {
                Symbol symbol = (Symbol)dnaIter.next();
                java.util.List symbols = ((AtomicSymbol)symbol).getSymbols();
                SymbolList symList = new SimpleSymbolList(alphabet, symbols);

                //double score = score(symList, (ScoreType)null, 0);
                //double score = score(symList, null);
                double score = score(symList, (ScoreType)null);

                relAffinityLHM.put(symList, score);

                if (score > highestScore) {
                    highestScore = score;
                }
            }

            // normalize all the scores by the highest score
            for (SymbolList symList : relAffinityLHM.keySet()) {
                double quotient = relAffinityLHM.get(symList) / highestScore;
                relAffinityLHM.put(symList, quotient);
            }

        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

        return(relAffinityLHM);
    }

    public Table getRelAffinitySphere(
        ScoreType scoreType,
        double featureThreshold,
        double relAffinityThreshold)
    {
        return(getRelAffinitySphere(scoreType, featureThreshold, relAffinityThreshold, null));
    }

    public Table getRelAffinitySphere(
        ScoreType scoreType,
        double featureThreshold,
        double relAffinityThreshold,
        ReversibleTranslationTable aComplementTable)
    {
        Table allScoresTable = null;
        try {
            LinkedHashMap<SymbolList, Double> relAffinityLHM = getAllScores();
            allScoresTable = new Table(relAffinityLHM);
            allScoresTable.sort(1, false); // descending

            // remove all entires < threshold
            int removeIndex = 0;
            for (int i=0; i < allScoresTable.rows(); i++) {
                if (((Double)allScoresTable.getElement(i, 1)).doubleValue() < relAffinityThreshold) {
                    removeIndex = i;
                    break;
                }
            }
            allScoresTable.removeRows(removeIndex, allScoresTable.rows());
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

        return(allScoresTable);
    }

//     public LinkedHashMap<SymbolList, Double> getRelAffinitySphere(
//         ScoreType scoreType,
//         double featureThreshold,
//         double relAffinityThreshold,
//         ReversibleTranslationTable aComplementTable)
//     {
//         LinkedHashMap<SymbolList, Double> relAffinityLHM = null;

//         try {
//             Table allScoresTable = new Table(getAllScores());
//             allScoresTable.sort(1, false); // descending

//             // remove all entires < threshold
//             int removeIndex = 0;
//             for (int i=0; i < allScoresTable.rows(); i++) {
//                 if (((Double)allScoresTable.getElement(i, 1)).doubleValue() < relAffinityThreshold) {
//                     removeIndex = i;
//                     break;
//                 }
//             }
//             allScoresTable.removeRows(removeIndex, allScoresTable.rows());

//             relAffinityLHM = allScoresTable.toMap(0, 1);
//         }
//         catch (Exception ex) {
//             ex.printStackTrace();
//         }

//         return(relAffinityLHM);
//     }

//     public LinkedHashMap<SymbolList, Double> getRelAffinitySphere(
//         ScoreType scoreType,
//         double featureThreshold,
//         double relAffinityThreshold,
//         ReversibleTranslationTable aComplementTable)
//     {

//         LinkedHashMap<SymbolList, Double> relAffinityLHM = new LinkedHashMap<SymbolList, Double>();
//         LinkedHashMap<SymbolList, Double> relAffinityRevCompsLHM = new LinkedHashMap<SymbolList, Double>();

//         try {

//             double highestScore = getHighestScore();

//             //SymbolList seedSymList = DNATools.createDNA("ATGATTAATGGC");
//             //SymbolList seedSymList = DNATools.createDNA("ATGATTAATGAC");
//             SymbolList seedSymList = WeightMatrixTools.getMostLikelySymList(posStrandPWM);
//             SymbolList revCompSeedSymList = SymbolListTools.reverseComplement(seedSymList, aComplementTable);

//             relAffinityLHM.put(seedSymList, 1.0);
//             relAffinityLHM.put(revCompSeedSymList, 1.0);

//             relAffinityRevCompsLHM.put(seedSymList, 1.0);
//             relAffinityRevCompsLHM.put(revCompSeedSymList, 1.0);

//             // iterate length(weightMatrix) number of times
//             for (int j = 0; j < this.columns; j++) {
//                 ArrayList<SymbolList> seedList = new ArrayList<SymbolList>(relAffinityLHM.keySet());

//                 // iterate over each symbolList currently in the set
//                 for (SymbolList seed : seedList) {

//                     // iterate over each position
//                     // 1-based indexing
//                     //for (int i = 1; i <= seed.length(); i++) {
//                     for (int i = j; i <= j; i++) {
//                         //for (int i = 1; i <= 5; i++) {

//                         // iterate over every symbol at this position
//                         for (Iterator dnaIter = ((FiniteAlphabet)alphabet).iterator(); dnaIter.hasNext(); ) {
//                             Symbol symbol = (Symbol)dnaIter.next();

//                             // if symbol is different from the seed symbol then create new sequence and test and add to the set
//                             if (symbol != seed.symbolAt(i)) {
//                                 //Symbol[] symbolArray = {symbol};
//                                 //make a copy of the seed
//                                 SimpleSymbolList editedSeed = new SimpleSymbolList(seed);
//                                 //Edit edit = new Edit(i, 1, new SimpleSymbolList(symbolArray, 1, alphabet));
//                                 Edit edit = new Edit(i, alphabet, symbol);
//                                 editedSeed.edit(edit);
//                                 //String seqName = new String("mutation-"+i+symbol.getName().substring(0,1));

//                                 double revCompSimilarity = WeightMatrixTools.similarity(posStrandPWM, negStrandPWM);

//                                 //double score = score(posStrandPWM, editedSeed, 0);

//                                 double score = score(
//                                     editedSeed,
//                                     null);
// //                                 double score = score(
// //                                     editedSeed,
// //                                     scoreType,
// //                                     featureThreshold);

//                                 score = score / highestScore;

//                                 if (score >= relAffinityThreshold) {
//                                     if (!relAffinityLHM.containsKey(editedSeed)) {
//                                         relAffinityLHM.put(editedSeed, score);
//                                     }

//                                     if (negStrandPWM != null) {
//                                         SymbolList revComp = SymbolListTools.reverseComplement(editedSeed, aComplementTable);
//                                         if (!relAffinityRevCompsLHM.containsKey(revComp)) {
//                                             relAffinityRevCompsLHM.put(revComp, score);
//                                         }
//                                     }

//                                 }

//                                 //                                 if ((negStrandPWM != null) && ((score = score(negStrandPWM, editedSeed, 0)) >= relAffinityThreshold)) {
//                                 //                                     SymbolList revComp = SymbolListTools.reverseComplement(editedSeed, aComplementTable);
//                                 //                                     if (!relAffinityLHM.containsKey(revComp)) {
//                                 //                                         relAffinityLHM.put(revComp, score);
//                                 //                                     }

//                                 //                                 }

//                             }
//                         }
//                     }
//                 }
//             }

//         }
//         catch (Exception ex) {
//             ex.printStackTrace();
//         }

//         relAffinityLHM.putAll(relAffinityRevCompsLHM);
//         return(relAffinityLHM);
//     }

    // strand should be POS or NEG
    public ArrayList<WeightMatrixFeature> getNNDDFeatures(
        SymbolList aSymbolList,
        WeightMatrixTools.BindingStrand aStrand,
        int startPos,
        int startColumn,
        int endColumn) {

        ArrayList<WeightMatrixFeature> features = new ArrayList<WeightMatrixFeature>();

        try {

//             if (featureKeyToMultFeatureMap == null) {
//                 return(features);
//             }

            //SymbolTokenization parser = DNATools.getDNA().getTokenization("default");
            //System.out.println("startPos="+startPos+"; startColumn="+startColumn+"; endColumn="+endColumn+"; seq="+aSymbolList.seqString()+"\n");

            for (int pos1 = startColumn; pos1 < endColumn-1; pos1++) { // from 0 to (L-1)
                for (int pos2 = pos1 + 1; pos2 <= pos1 + 1; pos2++) { // Nearest-Neighbor
                    Symbol[] featureArray = new Symbol[(columns*2)-1];
                    featureArray[(pos1*2)] = aSymbolList.symbolAt(startPos+pos1+1);
                    featureArray[(pos2*2)] = aSymbolList.symbolAt(startPos+pos2+1);

                    if (aStrand == WeightMatrixTools.BindingStrand.NEG) {
                        featureArray = SymbolListTools.reverseComplement(featureArray, DNATools.complementTable());
                        //featureArray = ArrayTools.reverse(featureArray);
                    }

                    //System.out.println("featureArray="+StringTools.toString(featureArray, parser, "", "-")+"\n");

                    FeatureKey featureKey = new FeatureKey(featureArray);

//                     if (featureKeyToMultFeatureMap.isEmpty()) {
//                         System.out.println("featureKeyToMultFeatureMap is empty!");
//                     }

//                     int keyNum=1;
//                     for (FeatureKey aFeatureKey: featureKeyToMultFeatureMap.keySet()) {
//                         System.out.println("FeatureKey "+keyNum+" ="+StringTools.toString(aFeatureKey.getModsArray(), parser, "", "-")+"  ; hashcode="+aFeatureKey.hashCode());
//                         keyNum++;
//                     }

                    WeightMatrixFeature feature = (WeightMatrixFeature)featureKeyToMultFeatureMap.get(featureKey);
                    if (feature != null) {
                        //System.out.println("feature found.");
                        features.add(feature);
                    }
//                     else {
//                         System.out.print("feature NOT found!");
//                         System.out.println("  featureArray="+StringTools.toString(featureArray, parser, "", "-")+"  ; hashcode="+featureKey.hashCode()+"\n");
//                     }
                }
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

        return(features);
    }

    // strand should be POS or NEG
    public ArrayList<WeightMatrixFeature> getNNNDDFeatures(
        SymbolList aSymbolList,
        WeightMatrixTools.BindingStrand aStrand,
        int startPos,
        int startColumn,
        int endColumn) {
        ArrayList<WeightMatrixFeature> features = new ArrayList<WeightMatrixFeature>();

        try {

//             if (featureKeyToMultFeatureMap == null) {
//                 return(features);
//             }

            //SymbolTokenization parser = DNATools.getDNA().getTokenization("default");
            //System.out.println("startPos="+startPos+"; startColumn="+startColumn+"; endColumn="+endColumn+"; seq="+aSymbolList.seqString()+"\n");

            for (int pos1 = startColumn; pos1 < endColumn-2; pos1++) { // from 0 to (L-1)
                for (int pos2 = pos1 + 2; pos2 < endColumn; pos2++) { // Nearest-Neighbor
                    Symbol[] featureArray = new Symbol[(columns*2)-1];
                    featureArray[(pos1*2)] = aSymbolList.symbolAt(startPos+pos1+1);
                    featureArray[(pos2*2)] = aSymbolList.symbolAt(startPos+pos2+1);

                    if (aStrand == WeightMatrixTools.BindingStrand.NEG) {
                        featureArray = SymbolListTools.reverseComplement(featureArray, DNATools.complementTable());
                        //featureArray = ArrayTools.reverse(featureArray);
                    }

                    //System.out.println("featureArray="+StringTools.toString(featureArray, parser, "", "-")+"\n");

                    FeatureKey featureKey = new FeatureKey(featureArray);

//                     if (featureKeyToMultFeatureMap.isEmpty()) {
//                         System.out.println("featureKeyToMultFeatureMap is empty!");
//                     }

//                     int keyNum=1;
//                     for (FeatureKey aFeatureKey: featureKeyToMultFeatureMap.keySet()) {
//                         System.out.println("FeatureKey "+keyNum+" ="+StringTools.toString(aFeatureKey.getModsArray(), parser, "", "-")+"  ; hashcode="+aFeatureKey.hashCode());
//                         keyNum++;
//                     }

                    WeightMatrixFeature feature = (WeightMatrixFeature)featureKeyToMultFeatureMap.get(featureKey);
                    if (feature != null) {
                        //System.out.println("feature found.");
                        features.add(feature);
                    }
//                     else {
//                         System.out.print("feature NOT found!");
//                         System.out.println("  featureArray="+StringTools.toString(featureArray, parser, "", "-")+"  ; hashcode="+featureKey.hashCode()+"\n");
//                     }
                }
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

        return(features);
    }


    // strand should be POS or NEG
    public ArrayList<WeightMatrixFeature> getNNTDFeatures(
        SymbolList aSymbolList,
        WeightMatrixTools.BindingStrand aStrand,
        int startPos,
        int startColumn,
        int endColumn) {
        ArrayList<WeightMatrixFeature> features = new ArrayList<WeightMatrixFeature>();

        try {

//             if (featureKeyToMultFeatureMap == null) {
//                 return(features);
//             }

            //SymbolTokenization parser = DNATools.getDNA().getTokenization("default");
            //System.out.println("startPos="+startPos+"; startColumn="+startColumn+"; endColumn="+endColumn+"; seq="+aSymbolList.seqString()+"\n");

            for (int pos1 = startColumn; pos1 < endColumn-2; pos1++) { // from 0 to (L-1)
                for (int pos2 = pos1 + 1; pos2 <= pos1 + 1; pos2++) { // Nearest-Neighbor
                    for (int pos3 = pos2 + 1; pos3 <= pos2 + 1; pos3++) { // Nearest-Neighbor
                        Symbol[] featureArray = new Symbol[(columns*2)-1];
                        featureArray[(pos1*2)] = aSymbolList.symbolAt(startPos+pos1+1);
                        featureArray[(pos2*2)] = aSymbolList.symbolAt(startPos+pos2+1);
                        featureArray[(pos3*2)] = aSymbolList.symbolAt(startPos+pos3+1);

                        if (aStrand == WeightMatrixTools.BindingStrand.NEG) {
                            featureArray = SymbolListTools.reverseComplement(featureArray, DNATools.complementTable());
                            //featureArray = ArrayTools.reverse(featureArray);
                        }

                        //System.out.println("featureArray="+StringTools.toString(featureArray, parser, "", "-")+"\n");

                        FeatureKey featureKey = new FeatureKey(featureArray);

                        //                     if (featureKeyToMultFeatureMap.isEmpty()) {
                        //                         System.out.println("featureKeyToMultFeatureMap is empty!");
                        //                     }

                        //                     int keyNum=1;
                        //                     for (FeatureKey aFeatureKey: featureKeyToMultFeatureMap.keySet()) {
                        //                         System.out.println("FeatureKey "+keyNum+" ="+StringTools.toString(aFeatureKey.getModsArray(), parser, "", "-")+"  ; hashcode="+aFeatureKey.hashCode());
                        //                         keyNum++;
                        //                     }

                        WeightMatrixFeature feature = (WeightMatrixFeature)featureKeyToMultFeatureMap.get(featureKey);
                        if (feature != null) {
                            //System.out.println("feature found.");
                            features.add(feature);
                        }
                        //                     else {
                        //                         System.out.print("feature NOT found!");
                        //                         System.out.println("  featureArray="+StringTools.toString(featureArray, parser, "", "-")+"  ; hashcode="+featureKey.hashCode()+"\n");
                        //                     }
                    }
                }
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

        return(features);
    }

    // score WITH the weights
    public double psamScore(SymbolList posStrandSeq, double[][] weights) {
        if (mandatoryColumns == null) {
            mandatoryColumns = new boolean[2][];
        }

        double score = WeightMatrixTools.score(
            posStrandPWM,
            negStrandPWM,
            posStrandSeq,
            weights,
            mandatoryColumns,
            strand,
            calc,
            eToMu,
            nonSpecKa,
            revCompSimilarity);

        return(score);
    }

    // score WITHOUT the weights
    public double psamScore(SymbolList posStrandSeq) {
        return(psamScore(posStrandSeq, this.strand));
    }

    // score WITHOUT the weights
    public double psamScore(SymbolList posStrandSeq, WeightMatrixTools.BindingStrand aStrand) {
        double score = WeightMatrixTools.score(
            posStrandPWM,
            negStrandPWM,
            posStrandSeq,
            //positionalWeights,
            //mandatoryColumns,
            aStrand,
            calc,
            eToMu,
            nonSpecKa,
            revCompSimilarity);

        return(score);
    }

    // score WITH the weights
    public double oldPsamScore(SymbolList posStrandSeq, double[][] weights) {
        if (mandatoryColumns == null) {
            mandatoryColumns = new boolean[2][];
        }

        double score = WeightMatrixTools.score(
            oldAMposStrandPWM,
            oldAMnegStrandPWM,
            posStrandSeq,
            weights,
            mandatoryColumns,
            strand,
            calc,
            eToMu,
            nonSpecKa,
            revCompSimilarity);

        return(score);
    }

    // score WITHOUT the weights
    public double oldPsamScore(SymbolList posStrandSeq) {
        double score = WeightMatrixTools.score(
            oldAMposStrandPWM,
            oldAMnegStrandPWM,
            posStrandSeq,
            //positionalWeights,
            //mandatoryColumns,
            strand,
            calc,
            eToMu,
            nonSpecKa,
            revCompSimilarity);

        return(score);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // score() = Multiplicative Model
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

//     // score WITH the weights
//     // use the mandatory columns of each feature
//     public double score(SymbolList posStrandSeq, double[][] weights, FeatureKey withholdFeatureKey) {

//         this.positionalWeights = weights;

//         //double score = WeightMatrixTools.score(
//         double score = score(
//             this,
//             posStrandSeq,
//             positionalWeights,
//             strand,
//             calc,
//             eToMu,
//             nonSpecKa,
//             //this.getRevCompSimilarity(),
//             revCompSimilarity,
//             withholdFeatureKey);


//         if (score < 0) {
//             if (nonSpecKa != 0) {
//                 score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
//             }
//             else {
//                 if (mandatoryColumns == null) {
//                     //mandatoryColumns = new boolean[2][this.columns];
//                     mandatoryColumns = new boolean[2][];
//                 }
//                 // score with PSAM only
//                 score = WeightMatrixTools.score(
//                     posStrandPWM,
//                     negStrandPWM,
//                     posStrandSeq,
//                     positionalWeights,
//                     mandatoryColumns,
//                     strand,
//                     calc,
//                     eToMu,
//                     nonSpecKa,
//                     revCompSimilarity);

//             }
//         }

//         return(score);
//     }

//     // score WITHOUT weights
//     // does not use the mandatory columns of each feature
//     public double score(SymbolList posStrandSeq, FeatureKey withholdFeatureKey, int startPos) {


//         //double score = WeightMatrixTools.score(
//         double score = score(
//             this,
//             posStrandSeq,
//             strand,
//             calc,
//             eToMu,
//             nonSpecKa,
//             //this.getRevCompSimilarity(),
//             revCompSimilarity,
//             startPos,
//             withholdFeatureKey);


//         if (score < 0) {
//             if (nonSpecKa != 0) {
//                 score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
//             }
//             else {
//                 if (mandatoryColumns == null) {
//                     //mandatoryColumns = new boolean[2][this.columns];
//                     mandatoryColumns = new boolean[2][];
//                 }
//                 // score with PSAM only
//                 score = WeightMatrixTools.score(
//                     posStrandPWM,
//                     negStrandPWM,
//                     posStrandSeq,
//                     strand,
//                     calc,
//                     eToMu,
//                     nonSpecKa,
//                     revCompSimilarity,
//                     startPos);

//             }
//         }

//         return(score);
//     }

//     // score WITHOUT weights
//     // does not use the mandatory columns of each feature
//     public double score(SymbolList posStrandSeq, FeatureKey withholdFeatureKey) {


//         //double score = WeightMatrixTools.score(
//         double score = score(
//             this,
//             posStrandSeq,
//             strand,
//             calc,
//             eToMu,
//             nonSpecKa,
//             //this.getRevCompSimilarity(),
//             revCompSimilarity,
//             withholdFeatureKey);


//         if (score < 0) {
//             if (nonSpecKa != 0) {
//                 score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
//             }
//             else {
//                 if (mandatoryColumns == null) {
//                     //mandatoryColumns = new boolean[2][this.columns];
//                     mandatoryColumns = new boolean[2][];
//                 }
//                 // score with PSAM only
//                 score = WeightMatrixTools.score(
//                     posStrandPWM,
//                     negStrandPWM,
//                     posStrandSeq,
//                     strand,
//                     calc,
//                     eToMu,
//                     nonSpecKa,
//                     revCompSimilarity);

//             }
//         }

//         return(score);
//     }


//     public double score(
//         FeaturedWeightMatrix anFWM,
//         SymbolList aSeq,
//         WeightMatrixTools.BindingStrand strand,
//         int startPos,
//         boolean[] mandatoryCols,
//         FeatureKey withholdFeatureKey) {

//         double logScore = 0;
//         try {
//             WeightMatrix pwm = anFWM.getPWM(strand);

//             // set the startColumn
//             int startColumn = 0;
//             if (startPos < 0) {
//                 startColumn = -1 * startPos;
//             }

//             // set the endColumn
//             int seqLength = aSeq.length();
//             int endColumn = anFWM.columns();
//             if (startPos + endColumn > seqLength) {
//                 endColumn = seqLength - startPos;
//             }

//             if (mandatoryCols != null) {
//                 // if mandatory column skipped, then return 0
//                 for (int c = 0; c < startColumn; c++) {
//                     if (mandatoryCols[c]) {
//                         return(0.0);
//                     }
//                 }
//                 for (int c = endColumn; c < pwm.columns(); c++) {
//                     if (mandatoryCols[c]) {
//                         return(0.0);
//                     }
//                 }
//             }

//             // calculate score
//             for (int c = startColumn; c < endColumn; c++) {
//                 double prob = pwm.getColumn(c).getWeight(aSeq.symbolAt(startPos+c+1));
//                 if (prob == 0.0) {
//                     return(0.0);
//                 }
//                 logScore += Math.log(prob);
//             }

//             //System.out.println("Hi 1");

//             if (featureKeyToMultFeatureMap != null) {
//                 // iterate through all the features found in this window
//                 //ArrayList<WeightMatrixFeature> nNDDFeatures = anFWM.getNNDDFeatures(aSeq, strand, startColumn, endColumn);
//                 //ArrayList<WeightMatrixFeature> nNDDFeatures = anFWM.getNNDDFeatures(aSeq, strand, startPos+startColumn, endColumn);

//                 ArrayList<WeightMatrixFeature> features = new ArrayList<WeightMatrixFeature>();

//                 features.addAll(anFWM.getNNDDFeatures(aSeq, strand, startPos, startColumn, endColumn));

//                 //features.addAll(anFWM.getNNNDDFeatures(aSeq, strand, startPos, startColumn, endColumn));

// //                 features.addAll(anFWM.getNNTDFeatures(aSeq, strand, startPos, startColumn, endColumn));

//                 //System.out.println("Hi 2");

//                 for (WeightMatrixFeature feature : features) {
//                     //boolean[] mandatoryColumns = feature.getMandatoryColumns();

//                     //System.out.println("Hi 3");

//                     // if should withhold this feature then continue to next feature
//                     if ((withholdFeatureKey!=null) && (withholdFeatureKey==feature.getKey())) {
//                         continue;
//                     }

//                     double featureRelAffinity = feature.getRelativeAffinity();

//                     //System.out.println(feature.toString()+"\n");
//                     //System.out.println("featureRelAffinity="+featureRelAffinity+"\n");

// //                     SymbolTokenization parser = pwm.getAlphabet().getTokenization("default");
// //                     System.out.println("window="+aSeq.subStr(startPos+startColumn+1, startPos+endColumn)+"; feature="+StringTools.toString(feature.getModsArray(), parser, "", "-")+"; strand="+strand+"; seq="+aSeq.seqString()+"\n");

//                     if (featureRelAffinity == 0.0) {
//                         return(0.0);
//                     }
//                     logScore += Math.log(featureRelAffinity);
//                 }
//             }

//         }
//         catch (Exception e) {
//             e.printStackTrace();
//         }
//         return(Math.exp(logScore));
//     }

//     // weights are applied after occupancy calculation
//     public double score(
//         FeaturedWeightMatrix anFWM,
//         SymbolList posStrandSeq,
//         double[][] weights,
//         WeightMatrixTools.BindingStrand strand,
//         WeightMatrixTools.BothStrandsCalc calc,
//         double eToMu,
//         double nonSpecKa,
//         double revCompSimilarity,
//         int startPos,
//         FeatureKey withholdFeatureKey) {

//         double posScore = 0;
//         double negScore = 0;
//         int weightsIndex = startPos + anFWM.columns() - 1;
//         double assocKa;

//         if (mandatoryColumns == null) {
//             mandatoryColumns = new boolean[2][];
//         }

//         // return the right score
//         switch (strand) {

//         case POS:
//             // score the window with the anFWM and posStrand
//             assocKa = score(anFWM, posStrandSeq, strand, startPos, mandatoryColumns[0], withholdFeatureKey);
//             if (eToMu == 0) {
//                 posScore = assocKa + nonSpecKa;
//             }
//             else {
//                 double prodKa = eToMu * (assocKa + nonSpecKa);
//                 posScore = prodKa / (prodKa + 1);
//             }
//             posScore = weights[0][weightsIndex] * posScore;
//             return(posScore);

//         case NEG:
//             // score the window with the negStrandPWM
//             assocKa = score(anFWM, posStrandSeq, strand, startPos, mandatoryColumns[1], withholdFeatureKey);
//             if (eToMu == 0) {
//                 negScore = assocKa + nonSpecKa;
//             }
//             else {
//                 double prodKa = eToMu * (assocKa + nonSpecKa);
//                 negScore = prodKa / (prodKa + 1);
//             }
//             negScore = weights[1][weightsIndex] * negScore;
//             return(negScore);

//         case BOTH:
//             double posAssocKa = score(anFWM, posStrandSeq, WeightMatrixTools.BindingStrand.POS, startPos, mandatoryColumns[0], withholdFeatureKey);
//             double negAssocKa = score(anFWM, posStrandSeq, WeightMatrixTools.BindingStrand.NEG, startPos, mandatoryColumns[1], withholdFeatureKey);
//             if (eToMu == 0) {
//                 posScore = posAssocKa + nonSpecKa;
//                 negScore = negAssocKa + nonSpecKa;
//             }
//             else {
//                 double posProdKa = eToMu * (posAssocKa + nonSpecKa);
//                 double negProdKa = eToMu * (negAssocKa + nonSpecKa);
//                 posScore = posProdKa / (posProdKa + 1);
//                 negScore = negProdKa / (negProdKa + 1);
//             }
//             posScore = weights[0][weightsIndex] * posScore;
//             negScore = weights[1][weightsIndex] * negScore;

//             switch (calc) {
//             case MAX:
//                 return(Math.max(posScore, negScore));
//             case SUM:
//                 return(posScore + negScore);
//             case UNION:
//                 return(posScore + negScore - (posScore*negScore));
//             case NORMED_SUM:
//                 return((posScore + negScore) / (revCompSimilarity + 1));
//             }
//             break;

//         }
//         return(0);
//     }

//     // weights are applied after occupancy calculation
//     public double score(
//         FeaturedWeightMatrix anFWM,
//         SymbolList posStrandSeq,
//         WeightMatrixTools.BindingStrand strand,
//         WeightMatrixTools.BothStrandsCalc calc,
//         double eToMu,
//         double nonSpecKa,
//         double revCompSimilarity,
//         int startPos,
//         FeatureKey withholdFeatureKey) {

//         double posScore = 0;
//         double negScore = 0;
//         double assocKa;

//         boolean[] emptyMandatoryColumns = null;

//         // return the right score
//         switch (strand) {

//         case POS:
//             // score the window with the anFWM and posStrand
//             assocKa = score(anFWM, posStrandSeq, strand, startPos, emptyMandatoryColumns, withholdFeatureKey);
//             if (eToMu == 0) {
//                 posScore = assocKa + nonSpecKa;
//             }
//             else {
//                 double prodKa = eToMu * (assocKa + nonSpecKa);
//                 posScore = prodKa / (prodKa + 1);
//             }
//             return(posScore);

//         case NEG:
//             // score the window with the negStrandPWM
//             assocKa = score(anFWM, posStrandSeq, strand, startPos, emptyMandatoryColumns, withholdFeatureKey);
//             if (eToMu == 0) {
//                 negScore = assocKa + nonSpecKa;
//             }
//             else {
//                 double prodKa = eToMu * (assocKa + nonSpecKa);
//                 negScore = prodKa / (prodKa + 1);
//             }
//             return(negScore);

//         case BOTH:
//             double posAssocKa = score(anFWM, posStrandSeq, WeightMatrixTools.BindingStrand.POS, startPos, emptyMandatoryColumns, withholdFeatureKey);
//             double negAssocKa = score(anFWM, posStrandSeq, WeightMatrixTools.BindingStrand.NEG, startPos, emptyMandatoryColumns, withholdFeatureKey);
//             if (eToMu == 0) {
//                 posScore = posAssocKa + nonSpecKa;
//                 negScore = negAssocKa + nonSpecKa;
//             }
//             else {
//                 double posProdKa = eToMu * (posAssocKa + nonSpecKa);
//                 double negProdKa = eToMu * (negAssocKa + nonSpecKa);
//                 posScore = posProdKa / (posProdKa + 1);
//                 negScore = negProdKa / (negProdKa + 1);
//             }

//             switch (calc) {
//             case MAX:
//                 return(Math.max(posScore, negScore));
//             case SUM:
//                 return(posScore + negScore);
//             case UNION:
//                 return(posScore + negScore - (posScore*negScore));
//             case NORMED_SUM:
//                 return((posScore + negScore) / (revCompSimilarity + 1));
//             }
//             break;

//         }
//         return(0);
//     }

//     public double score(FeaturedWeightMatrix anFWM, SymbolList posStrandSeq, WeightMatrixTools.BindingStrand strand, WeightMatrixTools.BothStrandsCalc calc, double eToMu, double nonSpecKa, double revCompSimilarity, FeatureKey withholdFeatureKey) {
//         double scoreSum = 0;
//         // Sweep down the posStrandSeq scoring each window
//         // 0-based indexing
//         for (int startPos = 0; startPos <= posStrandSeq.length() - anFWM.columns(); startPos++) {
//             scoreSum += score(anFWM, posStrandSeq, strand, calc, eToMu, nonSpecKa, revCompSimilarity, startPos, withholdFeatureKey);
//         }
//         return(scoreSum);
//     }

//     public double score(FeaturedWeightMatrix anFWM, SymbolList posStrandSeq, double[][] weights, WeightMatrixTools.BindingStrand strand, WeightMatrixTools.BothStrandsCalc calc, double eToMu, double nonSpecKa, double revCompSimilarity, FeatureKey withholdFeatureKey) {
//         double scoreSum = 0;
//         // Sweep down the posStrandSeq scoring each window
//         // 0-based indexing
//         for (int startPos = -1 * (anFWM.columns() - 1); startPos < posStrandSeq.length(); startPos++) {
//             scoreSum += score(anFWM, posStrandSeq, weights, strand, calc, eToMu, nonSpecKa, revCompSimilarity, startPos, withholdFeatureKey);
//         }
//         return(scoreSum);
//     }

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// score() = New AddModelScores
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

    // // score WITH the weights
    // // use the mandatory columns of each feature
    // public double score(SymbolList posStrandSeq, double[][] weights, FeatureKey withholdFeatureKey) {

    //     this.positionalWeights = weights;

    //     //double score = WeightMatrixTools.score(
    //     double score = score(
    //         posStrandSeq,
    //         positionalWeights,
    //         strand,
    //         calc,
    //         eToMu,
    //         nonSpecKa,
    //         //this.getRevCompSimilarity(),
    //         revCompSimilarity,
    //         withholdFeatureKey);


    //     if (score < 0) {
    //         if (nonSpecKa != 0) {
    //             score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
    //         }
    //         else {
    //             if (mandatoryColumns == null) {
    //                 //mandatoryColumns = new boolean[2][this.columns];
    //                 mandatoryColumns = new boolean[2][];
    //             }
    //             // score with PSAM only
    //             score = WeightMatrixTools.score(
    //                 posStrandPWM,
    //                 negStrandPWM,
    //                 posStrandSeq,
    //                 positionalWeights,
    //                 mandatoryColumns,
    //                 strand,
    //                 calc,
    //                 eToMu,
    //                 nonSpecKa,
    //                 revCompSimilarity);

    //         }
    //     }

    //     return(score);
    // }

    // // score WITHOUT weights
    // // does not use the mandatory columns of each feature
    // public double score(SymbolList posStrandSeq, FeatureKey withholdFeatureKey, int startPos) {


    //     //double score = WeightMatrixTools.score(
    //     double score = score(
    //         posStrandSeq,
    //         strand,
    //         calc,
    //         eToMu,
    //         nonSpecKa,
    //         //this.getRevCompSimilarity(),
    //         revCompSimilarity,
    //         startPos,
    //         withholdFeatureKey);


    //     if (score < 0) {
    //         if (nonSpecKa != 0) {
    //             score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
    //         }
    //         else {
    //             if (mandatoryColumns == null) {
    //                 //mandatoryColumns = new boolean[2][this.columns];
    //                 mandatoryColumns = new boolean[2][];
    //             }
    //             // score with PSAM only
    //             score = WeightMatrixTools.score(
    //                 posStrandPWM,
    //                 negStrandPWM,
    //                 posStrandSeq,
    //                 strand,
    //                 calc,
    //                 eToMu,
    //                 nonSpecKa,
    //                 revCompSimilarity,
    //                 startPos);

    //         }
    //     }

    //     return(score);
    // }

    // // score WITHOUT weights
    // // does not use the mandatory columns of each feature
    // public double score(SymbolList posStrandSeq, ScoreType scoreType) {
    //     return(score(posStrandSeq, (FeatureKey)null));
    // }

    // // score WITHOUT weights
    // // does not use the mandatory columns of each feature
    // public double score(SymbolList posStrandSeq, FeatureKey withholdFeatureKey) {


    //     //double score = WeightMatrixTools.score(
    //     double score = score(
    //         posStrandSeq,
    //         strand,
    //         calc,
    //         eToMu,
    //         nonSpecKa,
    //         //this.getRevCompSimilarity(),
    //         revCompSimilarity,
    //         withholdFeatureKey);


    //     if (score < 0) {
    //         if (nonSpecKa != 0) {
    //             score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
    //         }
    //         else {
    //             if (mandatoryColumns == null) {
    //                 //mandatoryColumns = new boolean[2][this.columns];
    //                 mandatoryColumns = new boolean[2][];
    //             }
    //             // score with PSAM only
    //             score = WeightMatrixTools.score(
    //                 posStrandPWM,
    //                 negStrandPWM,
    //                 posStrandSeq,
    //                 strand,
    //                 calc,
    //                 eToMu,
    //                 nonSpecKa,
    //                 revCompSimilarity);

    //         }
    //     }

    //     return(score);
    // }

    // public double score(
    //     //FeaturedWeightMatrix anFWM,
    //     SymbolList posStrandSeq,
    //     WeightMatrixTools.BindingStrand strand,
    //     WeightMatrixTools.BothStrandsCalc calc,
    //     double eToMu,
    //     double nonSpecKa,
    //     double revCompSimilarity,
    //     int startPos,
    //     FeatureKey withholdFeatureKey)
    // {
    //     double scoreSum = 0;
    //     double psamScore = 0;

    //     WeightMatrix addPosStrandPWM;
    //     WeightMatrix addNegStrandPWM;
    //     if (oldAMposStrandPWM == null) {
    //         addPosStrandPWM = posStrandPWM;
    //         addNegStrandPWM = negStrandPWM;
    //     }
    //     else {
    //         addPosStrandPWM = oldAMposStrandPWM;
    //         addNegStrandPWM = oldAMnegStrandPWM;
    //     }

    //     if (includePwm) {
    //         scoreSum += WeightMatrixTools.score(
    //             addPosStrandPWM,
    //             addNegStrandPWM,
    //             posStrandSeq,
    //             strand,
    //             calc,
    //             eToMu,
    //             nonSpecKa,
    //             revCompSimilarity,
    //             startPos);
    //         psamScore = scoreSum;
    //     }

    //     // add any features
    //     if (featureKeyToMultFeatureMap != null) {

    //         // set the startColumn
    //         int startColumn = 0;
    //         if (startPos < 0) {
    //             startColumn = -1 * startPos;
    //         }

    //         // set the endColumn
    //         int seqLength = posStrandSeq.length();
    //         int endColumn = this.columns;
    //         if (startPos + endColumn > seqLength) {
    //             endColumn = seqLength - startPos;
    //         }

    //         ArrayList<WeightMatrixFeature> features = new ArrayList<WeightMatrixFeature>();
    //         features.addAll(getNNDDFeatures(posStrandSeq, strand, startPos, startColumn, endColumn));

    //         for (WeightMatrixFeature feature : features) {

    //             if ((withholdFeatureKey!=null) && (withholdFeatureKey==feature.getKey())) {
    //                 continue;
    //             }

    //             double featureRelAffinity = feature.getRelativeAffinity();

    //             double featureScore = WeightMatrixTools.score(
    //                 feature.getPosStrandPWM(),
    //                 feature.getNegStrandPWM(),
    //                 posStrandSeq,
    //                 strand,
    //                 calc,
    //                 eToMu,
    //                 0, //nonSpecKa,
    //                 feature.getRevCompSimilarity(),
    //                 startPos);

    //             scoreSum += featureRelAffinity * featureScore;
    //         }
    //     }

    //     if (scoreSum < 0) {
    //         if (nonSpecKa != 0) {
    //             scoreSum = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
    //         }
    //         else {
    //             scoreSum = psamScore;
    //         }
    //     }

    //     return(scoreSum);
    // }

    // public double score(
    //     //FeaturedWeightMatrix anFWM,
    //     SymbolList posStrandSeq,
    //     WeightMatrixTools.BindingStrand strand,
    //     WeightMatrixTools.BothStrandsCalc calc,
    //     double eToMu,
    //     double nonSpecKa,
    //     double revCompSimilarity,
    //     FeatureKey withholdFeatureKey)
    // {
    //     double scoreSum = 0;
    //     double psamScore = 0;

    //     WeightMatrix addPosStrandPWM;
    //     WeightMatrix addNegStrandPWM;
    //     if (oldAMposStrandPWM == null) {
    //         addPosStrandPWM = posStrandPWM;
    //         addNegStrandPWM = negStrandPWM;
    //     }
    //     else {
    //         addPosStrandPWM = oldAMposStrandPWM;
    //         addNegStrandPWM = oldAMnegStrandPWM;
    //     }

    //     // Sweep down the posStrandSeq scoring each window
    //     // 0-based indexing
    //     for (int startPos = 0; startPos <= posStrandSeq.length() - this.columns; startPos++) {

    //         if (includePwm) {
    //             scoreSum += WeightMatrixTools.score(
    //                 addPosStrandPWM,
    //                 addNegStrandPWM,
    //                 posStrandSeq,
    //                 strand,
    //                 calc,
    //                 eToMu,
    //                 nonSpecKa,
    //                 revCompSimilarity,
    //                 startPos);
    //             psamScore = scoreSum;
    //         }

    //         // add any features
    //         if (featureKeyToMultFeatureMap != null) {

    //             // set the startColumn
    //             int startColumn = 0;
    //             if (startPos < 0) {
    //                 startColumn = -1 * startPos;
    //             }

    //             // set the endColumn
    //             int seqLength = posStrandSeq.length();
    //             int endColumn = this.columns;
    //             if (startPos + endColumn > seqLength) {
    //                 endColumn = seqLength - startPos;
    //             }

    //             ArrayList<WeightMatrixFeature> features = new ArrayList<WeightMatrixFeature>();
    //             features.addAll(getNNDDFeatures(posStrandSeq, strand, startPos, startColumn, endColumn));

    //             for (WeightMatrixFeature feature : features) {

    //                 if ((withholdFeatureKey!=null) && (withholdFeatureKey==feature.getKey())) {
    //                     continue;
    //                 }

    //                 double featureRelAffinity = feature.getRelativeAffinity();

    //                 double featureScore = WeightMatrixTools.score(
    //                     feature.getPosStrandPWM(),
    //                     feature.getNegStrandPWM(),
    //                     posStrandSeq,
    //                     strand,
    //                     calc,
    //                     eToMu,
    //                     0, //nonSpecKa,
    //                     feature.getRevCompSimilarity(),
    //                     startPos);

    //                 scoreSum += featureRelAffinity * featureScore;
    //             }
    //         }
    //     }

    //     if (scoreSum < 0) {
    //         if (nonSpecKa != 0) {
    //             scoreSum = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
    //         }
    //         else {
    //             scoreSum = psamScore;
    //         }
    //     }

    //     return(scoreSum);
    // }

    // public double score(
    //     //FeaturedWeightMatrix anFWM,
    //     SymbolList posStrandSeq,
    //     double[][] weights,
    //     WeightMatrixTools.BindingStrand strand,
    //     WeightMatrixTools.BothStrandsCalc calc,
    //     double eToMu,
    //     double nonSpecKa,
    //     double revCompSimilarity,
    //     FeatureKey withholdFeatureKey)
    // {
    //     double scoreSum = 0;
    //     double psamScore = 0;

    //     WeightMatrix addPosStrandPWM;
    //     WeightMatrix addNegStrandPWM;
    //     if (oldAMposStrandPWM == null) {
    //         addPosStrandPWM = posStrandPWM;
    //         addNegStrandPWM = negStrandPWM;
    //     }
    //     else {
    //         addPosStrandPWM = oldAMposStrandPWM;
    //         addNegStrandPWM = oldAMnegStrandPWM;
    //     }

    //     // create an empty mandatoryColumns
    //     boolean[][] emptyMandatoryColumns = new boolean[2][this.columns];

    //     // Sweep down the posStrandSeq scoring each window
    //     // 0-based indexing
    //     for (int startPos = -1 * (this.columns - 1); startPos < posStrandSeq.length(); startPos++) {

    //         if (includePwm) {
    //             scoreSum += WeightMatrixTools.score(
    //                 addPosStrandPWM,
    //                 addNegStrandPWM,
    //                 posStrandSeq,
    //                 weights,
    //                 emptyMandatoryColumns,
    //                 strand,
    //                 calc,
    //                 eToMu,
    //                 nonSpecKa,
    //                 revCompSimilarity,
    //                 startPos);
    //             psamScore = scoreSum;
    //         }

    //         // add any features
    //         if (featureKeyToMultFeatureMap != null) {

    //             // set the startColumn
    //             int startColumn = 0;
    //             if (startPos < 0) {
    //                 startColumn = -1 * startPos;
    //             }

    //             // set the endColumn
    //             int seqLength = posStrandSeq.length();
    //             int endColumn = this.columns;
    //             if (startPos + endColumn > seqLength) {
    //                 endColumn = seqLength - startPos;
    //             }

    //             ArrayList<WeightMatrixFeature> features = new ArrayList<WeightMatrixFeature>();
    //             features.addAll(getNNDDFeatures(posStrandSeq, strand, startPos, startColumn, endColumn));

    //             for (WeightMatrixFeature feature : features) {

    //                 if ((withholdFeatureKey!=null) && (withholdFeatureKey==feature.getKey())) {
    //                     continue;
    //                 }

    //                 double featureRelAffinity = feature.getRelativeAffinity();

    //                 double featureScore = WeightMatrixTools.score(
    //                     feature.getPosStrandPWM(),
    //                     feature.getNegStrandPWM(),
    //                     posStrandSeq,
    //                     weights,
    //                     feature.getMandatoryColumns(),
    //                     strand,
    //                     calc,
    //                     eToMu,
    //                     0, //nonSpecKa
    //                     feature.getRevCompSimilarity(),
    //                     startPos);

    //                 scoreSum += featureRelAffinity * featureScore;
    //             }
    //         }
    //     }

    //     if (scoreSum < 0) {
    //         if (nonSpecKa != 0) {
    //             scoreSum = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
    //         }
    //         else {
    //             scoreSum = psamScore;
    //         }
    //     }

    //     return(scoreSum);
    // }


// /////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////////////////
// // score() = 0ld AddModelScores
// /////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////////////////


//     /////////////////////////////////////////////
//     // score WITH the weights
//     // use the mandatory columns of each feature
//     /////////////////////////////////////////////
//     public double score(SymbolList posStrandSeq, double[][] weights, FeatureKey blah) {
//         return(score(posStrandSeq, weights, null, 0));
//     }

//     public double score(SymbolList posStrandSeq, double[][] weights, ScoreType scoreType) {
//         return(score(posStrandSeq, weights, scoreType, 0));
//     }

//     public double score(SymbolList posStrandSeq, double[][] weights, ScoreType scoreType, double threshold) {

//         this.positionalWeights = weights;

//         //System.out.println("\nthreshold = "+threshold);
//         //System.out.println("\nweights = "+StringTools.toString(weights, "\n", "\t"));

//         // If the PSAM is an reverse-complement palindrome, then we should look at only the positive
//         // or negative strand (BUT NOT BOTH)
// //         if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.BOTH)) {
// //             //out.println("\n\nSince the motif is self reverse-complement, we should look only at the positive strand.");
// //             strand = WeightMatrixTools.BindingStrand.POS;
// //         }

//         // create an empty mandatoryColumns
//         boolean[][] mandatoryColumns = new boolean[2][this.columns];

//         // score the posStrandSeq with the posStrandPWM
//         double score = 0;
//         double psamScore = 0;

//         WeightMatrix addPosStrandPWM;
//         WeightMatrix addNegStrandPWM;
//         if (oldAMposStrandPWM == null) {
//             addPosStrandPWM = posStrandPWM;
//             addNegStrandPWM = negStrandPWM;
//         }
//         else {
//             addPosStrandPWM = oldAMposStrandPWM;
//             addNegStrandPWM = oldAMnegStrandPWM;
//         }

//         if (includePwm) {
//             if (scoreType != null) {
//                 score += WeightMatrixTools.score(
//                     addPosStrandPWM,
//                     addNegStrandPWM,
//                     posStrandSeq,
//                     positionalWeights,
//                     mandatoryColumns,
//                     strand,
//                     calc,
//                     scoreType,
//                     eToMu,
//                     nonSpecKa,
//                     revCompSimilarity);
//             }
//             else {
//                 score += WeightMatrixTools.score(
//                     addPosStrandPWM,
//                     addNegStrandPWM,
//                     posStrandSeq,
//                     positionalWeights,
//                     mandatoryColumns,
//                     strand,
//                     calc,
//                     eToMu,
//                     nonSpecKa,
//                     revCompSimilarity);
//             }
//             psamScore = score;
//             // System.out.println("\nPSAM score = "+score);
//         }

//         // Iterate over each feature and add their contribution to the affinity score
//         getFeatures(0, true);
//         if (featuresTable != null) {

//             // sort the table by the relative affinities if the table is in an unsorted state
//             sortIfNeeded();

//             double lastScore = score;

// //             ArrayList<WeightMatrixFeature> nNDDFeatures = getNNDDFeatures(posStrandSeq);
// //             for (WeightMatrixFeature feature : nNDDFeatures) {

// //             ArrayList<WeightMatrixFeature> features = new ArrayList<WeightMatrixFeature>();
// //             features.addAll(anFWM.getNNDDFeatures(aSeq, strand, startPos, startColumn, endColumn));

//             for (Iterator iter = featuresTable.iterator(); iter.hasNext(); ) {
//                 List rowList = (List)iter.next();
//                 WeightMatrixFeature feature = (WeightMatrixFeature)rowList.get(0);

//                 // skip the feature if strand==BOTH and isRevComp
//                 //             if (strand==WeightMatrixTools.BindingStrand.BOTH && feature.isRevComp()) {
//                 //                 continue;
//                 //             }

//                 double relativeAffinity = feature.getRelativeAffinity();

//                 // if the relativeAffinity is below the threshold then stop adding feature corrections
//                 if (Math.abs(relativeAffinity) < threshold) {
//                     System.out.println("\nRelative affinity "+relativeAffinity+" is less than the threshold "+threshold);
//                     break;
//                 }

//                 // score the posStrandSeq with the featurePosStrandPWM
//                 double featureScore;
//                 if (scoreType != null) {
//                     featureScore = WeightMatrixTools.score(
//                         feature.getPosStrandPWM(),
//                         feature.getNegStrandPWM(),
//                         posStrandSeq,
//                         positionalWeights,
//                         feature.getMandatoryColumns(),
//                         strand,
//                         calc,
//                         scoreType,
//                         eToMu,
//                         //nonSpecKa,
//                         0,
//                         feature.getRevCompSimilarity());
//                 }
//                 else {
//                     featureScore = WeightMatrixTools.score(
//                         feature.getPosStrandPWM(),
//                         feature.getNegStrandPWM(),
//                         posStrandSeq,
//                         positionalWeights,
//                         feature.getMandatoryColumns(),
//                         strand,
//                         calc,
//                         eToMu,
//                         //nonSpecKa,
//                         0,
//                         feature.getRevCompSimilarity());
//                 }

// //                 if (!feature.isGood()) {
// //                     relativeAffinity *= -1.0;
// //                 }
//                 score += relativeAffinity * featureScore;

//                 // This feature applies in this sequence context iff the score has changed
//                 if (score != lastScore) {
//                     // System.out.println("feature = "+feature.modsToString()+"\tupdated score = "+score);
//                 }
//                 lastScore = score;
//             }
//         }

//         // If our motif isRevCompPalindrome and we looked at only the POS strand, then double our score
// //         if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.POS)) {
// //             score *= 2.0;
// //         }

//         if (score < 0) {
//             if (nonSpecKa != 0) {
//                 score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
//             }
//             else {
//                 score = psamScore;
//             }
//         }

//         return(score);
//     }

//     ///////////////////////////////////////////
//     // score with a startPos
//     // use the mandatory columns of each feature
//     ///////////////////////////////////////////
//     public double score(SymbolList posStrandSeq, FeatureKey blah, int startPos) {
//         return(score(posStrandSeq, null, 0, startPos));
//     }

//     public double score(SymbolList posStrandSeq, ScoreType scoreType, int startPos) {
//         return(score(posStrandSeq, scoreType, 0, startPos));
//     }

//     public double score(SymbolList posStrandSeq, ScoreType scoreType, double threshold, int startPos) {

//         // If the PSAM is an reverse-complement palindrome, then we should look at only the positive
//         // or negative strand (BUT NOT BOTH)
// //         if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.BOTH)) {
// //             //out.println("\n\nSince the motif is self reverse-complement, we should look only at the positive strand.");
// //             strand = WeightMatrixTools.BindingStrand.POS;
// //         }

//         // score the posStrandSeq with the addPosStrandPWM
//         double score = 0;
//         double psamScore = 0;

//         WeightMatrix addPosStrandPWM;
//         WeightMatrix addNegStrandPWM;
//         if (oldAMposStrandPWM == null) {
//             addPosStrandPWM = posStrandPWM;
//             addNegStrandPWM = negStrandPWM;
//         }
//         else {
//             addPosStrandPWM = oldAMposStrandPWM;
//             addNegStrandPWM = oldAMnegStrandPWM;
//         }

//         if (includePwm) {
//             if (scoreType != null) {
//                 score += WeightMatrixTools.score(
//                     addPosStrandPWM,
//                     addNegStrandPWM,
//                     posStrandSeq,
//                     strand,
//                     calc,
//                     scoreType,
//                     eToMu,
//                     nonSpecKa,
//                     revCompSimilarity,
//                     startPos);
//             }
//             else {
//                 score += WeightMatrixTools.score(
//                     addPosStrandPWM,
//                     addNegStrandPWM,
//                     posStrandSeq,
//                     strand,
//                     calc,
//                     eToMu,
//                     nonSpecKa,
//                     revCompSimilarity,
//                     startPos);
//             }
//             psamScore = score;
//             // System.out.println("\nPSAM score = "+score);
//         }

//         // make an array of unity weights
//         double[][] weights = new double[2][posStrandSeq.length() + this.columns - 1];
//         Arrays.fill(weights[0], 1.0);
//         Arrays.fill(weights[1], 1.0);

//         // Iterate over each feature and add their contribution to the affinity score
//         getFeatures(0, true);
//         if (featuresTable != null) {

//             // sort the table by the relative affinities if the table is in an unsorted state
//             sortIfNeeded();

//             double lastScore = score;

//             for (Iterator iter = featuresTable.iterator(); iter.hasNext(); ) {
//                 List rowList = (List)iter.next();
//                 WeightMatrixFeature feature = (WeightMatrixFeature)rowList.get(0);

//                 // skip the feature if strand==BOTH and isRevComp
//                 //             if (strand==WeightMatrixTools.BindingStrand.BOTH && feature.isRevComp()) {
//                 //                 continue;
//                 //             }

//                 double relativeAffinity = feature.getRelativeAffinity();

//                 // if the relativeAffinity is below the threshold then stop adding feature corrections
//                 if (Math.abs(relativeAffinity) < threshold) {
//                     System.out.println("\nRelative affinity "+relativeAffinity+" is less than the threshold "+threshold);
//                     break;
//                 }

//                 // score the posStrandSeq with the featurePosStrandPWM
//                 double featureScore;
//                 if (scoreType != null) {
//                     featureScore = WeightMatrixTools.score(
//                         feature.getPosStrandPWM(),
//                         feature.getNegStrandPWM(),
//                         posStrandSeq,
//                         weights,
//                         feature.getMandatoryColumns(),
//                         strand,
//                         calc,
//                         scoreType,
//                         eToMu,
//                         //nonSpecKa,
//                         0,
//                         feature.getRevCompSimilarity(),
//                         startPos);
//                 }
//                 else {
//                     featureScore = WeightMatrixTools.score(
//                         feature.getPosStrandPWM(),
//                         feature.getNegStrandPWM(),
//                         posStrandSeq,
//                         weights,
//                         feature.getMandatoryColumns(),
//                         strand,
//                         calc,
//                         eToMu,
//                         //nonSpecKa,
//                         0,
//                         feature.getRevCompSimilarity(),
//                         startPos);
//                 }

// //                 if (!feature.isGood()) {
// //                     relativeAffinity *= -1.0;
// //                 }
//                 score += relativeAffinity * featureScore;

//                 // This feature applies in this sequence context iff the score has changed
//                 if (score != lastScore) {
//                     // System.out.println("feature = "+feature.modsToString()+"\tupdated score = "+score);
//                 }
//                 lastScore = score;
//             }
//         }

//         // If our motif isRevCompPalindrome and we looked at only the POS strand, then double our score
// //         if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.POS)) {
// //             score *= 2.0;
// //         }

//         if (score < 0) {
//             if (nonSpecKa != 0) {
//                 score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
//             }
//             else {
//                 score = psamScore;
//             }
//         }

//         return(score);
//     }

//     ////////////////////////////////////////////////////
//     // score WITHOUT weights
//     // does not use the mandatory columns of each feature
//     ////////////////////////////////////////////////////
//     public double score(SymbolList posStrandSeq, FeatureKey blah) {
//         return(score(posStrandSeq, (ScoreType)null, 0));
//     }

//     public double score(SymbolList posStrandSeq, ScoreType scoreType) {
//         return(score(posStrandSeq, scoreType, 0));
//     }

//     public double score(SymbolList posStrandSeq, ScoreType scoreType, double threshold) {

//         // If the PSAM is an reverse-complement palindrome, then we should look at only the positive
//         // or negative strand (BUT NOT BOTH)
// //         if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.BOTH)) {
// //             //out.println("\n\nSince the motif is self reverse-complement, we should look only at the positive strand.");
// //             strand = WeightMatrixTools.BindingStrand.POS;
// //         }

//         // score the posStrandSeq with the posStrandPWM
//         double score = 0;
//         double psamScore = 0;

//         WeightMatrix addPosStrandPWM;
//         WeightMatrix addNegStrandPWM;
//         if (oldAMposStrandPWM == null) {
//             addPosStrandPWM = posStrandPWM;
//             addNegStrandPWM = negStrandPWM;
//         }
//         else {
//             addPosStrandPWM = oldAMposStrandPWM;
//             addNegStrandPWM = oldAMnegStrandPWM;
//         }

//         if (includePwm) {
//             if (scoreType != null) {
//                 score += WeightMatrixTools.score(
//                     addPosStrandPWM,
//                     addNegStrandPWM,
//                     posStrandSeq,
//                     strand,
//                     calc,
//                     scoreType,
//                     eToMu,
//                     nonSpecKa,
//                     revCompSimilarity);
//             }
//             else {
//                 score += WeightMatrixTools.score(
//                     addPosStrandPWM,
//                     addNegStrandPWM,
//                     posStrandSeq,
//                     strand,
//                     calc,
//                     eToMu,
//                     nonSpecKa,
//                     revCompSimilarity);
//             }
//             psamScore = score;
//             //System.out.println("\nPSAM score = "+score);
//         }

//         // Iterate over each feature and add their contribution to the affinity score
//         getFeatures(0, true);
//         if (featuresTable != null) {

//             // sort the table by the relative affinities if the table is in an unsorted state
//             sortIfNeeded();

//             double lastScore = score;

//             for (Iterator iter = featuresTable.iterator(); iter.hasNext(); ) {
//                 List rowList = (List)iter.next();
//                 WeightMatrixFeature feature = (WeightMatrixFeature)rowList.get(0);

//                 // skip the feature if strand==BOTH and isRevComp
//                 //             if (strand==WeightMatrixTools.BindingStrand.BOTH && feature.isRevComp()) {
//                 //                 continue;
//                 //             }

//                 double relativeAffinity = feature.getRelativeAffinity();

//                 // if the relativeAffinity is below the threshold then stop adding feature corrections
//                 if (Math.abs(relativeAffinity) < threshold) {
//                     System.out.println("\nRelative affinity "+relativeAffinity+" is less than the threshold "+threshold);
//                     break;
//                 }

//                 // score the posStrandSeq with the featurePosStrandPWM
//                 double featureScore;
//                 if (scoreType != null) {
//                     featureScore = WeightMatrixTools.score(
//                         feature.getPosStrandPWM(),
//                         feature.getNegStrandPWM(),
//                         posStrandSeq,
//                         strand,
//                         calc,
//                         scoreType,
//                         eToMu,
//                         //nonSpecKa,
//                         0,
//                         feature.getRevCompSimilarity());
//                 }
//                 else {
//                     featureScore = WeightMatrixTools.score(
//                         feature.getPosStrandPWM(),
//                         feature.getNegStrandPWM(),
//                         posStrandSeq,
//                         strand,
//                         calc,
//                         eToMu,
//                         //nonSpecKa,
//                         0,
//                         feature.getRevCompSimilarity());
//                 }

// //                 if (!feature.isGood()) {
// //                     relativeAffinity *= -1.0;
// //                 }
//                 score += relativeAffinity * featureScore;

//                 // This feature applies in this sequence context iff the score has changed
//                 if (score != lastScore) {
//                     //System.out.println("feature = "+feature.modsToString()+"\tupdated score = "+score);
//                 }
//                 lastScore = score;
//             }
//         }

//         // If our motif isRevCompPalindrome and we looked at only the POS strand, then double our score
// //         if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.POS)) {
// //             score *= 2.0;
// //         }

//         if (score < 0) {
//             if (nonSpecKa != 0) {
//                 score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
//             }
//             else {
//                 score = psamScore;
//             }
//         }

//         return(score);
//     }

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// score() = call multModelScore()
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////////
    // a threshold is no longer necessary with the multModel!!!!!
    /////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////
    // score WITH the weights
    // use the mandatory columns of each feature
    /////////////////////////////////////////////
    public double score(SymbolList posStrandSeq, double[][] weights, FeatureKey withholdFeatureKey) {
        return(multModelScore(posStrandSeq, weights, withholdFeatureKey)); // done
    }

    public double score(SymbolList posStrandSeq, double[][] weights, ScoreType scoreType) {
        return(multModelScore(posStrandSeq, weights, null)); // ignores scoreType!!!
        //return(multModelScore(posStrandSeq, weights, scoreType, 0));
    }

    public double score(SymbolList posStrandSeq, double[][] weights, ScoreType scoreType, double threshold) {
        return(multModelScore(posStrandSeq, weights, null)); // ignores scoreType and threshold!!!
        //return(multModelScore(posStrandSeq, weights, scoreType, threshold));
    }

    ///////////////////////////////////////////
    // score with a startPos
    // use the mandatory columns of each feature
    ///////////////////////////////////////////
    public double score(SymbolList posStrandSeq, FeatureKey withholdFeatureKey, int startPos) {
        return(multModelScore(posStrandSeq, withholdFeatureKey, startPos)); // done
    }

    public double score(SymbolList posStrandSeq, ScoreType scoreType, int startPos) {
        return(multModelScore(posStrandSeq, null, startPos)); // ignores scoreType!!!!
        //return(multModelScore(posStrandSeq, scoreType, 0, startPos));
    }

    public double score(SymbolList posStrandSeq, ScoreType scoreType, double threshold, int startPos) {
        return(multModelScore(posStrandSeq, null, startPos)); // ignores scoreType and threshold!!!!
        //return(multModelScore(posStrandSeq, scoreType, threshold, startPos));
    }

    ////////////////////////////////////////////////////
    // score WITHOUT weights
    // does not use the mandatory columns of each feature
    ////////////////////////////////////////////////////
    public double score(SymbolList posStrandSeq, FeatureKey withholdFeatureKey) {
        return(multModelScore(posStrandSeq, withholdFeatureKey)); // done
    }

    public double score(SymbolList posStrandSeq, ScoreType scoreType) {
        return(multModelScore(posStrandSeq, scoreType)); // done
    }

    public double score(SymbolList posStrandSeq, ScoreType scoreType, double threshold) {
        return(multModelScore(posStrandSeq, scoreType)); // ignores threshold!!!!
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Multiplicative Model
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    // score WITH the weights
    // use the mandatory columns of each feature
    public double multModelScore(SymbolList posStrandSeq, double[][] weights, FeatureKey withholdFeatureKey) {

        this.positionalWeights = weights;

        //double score = WeightMatrixTools.score(
        double score = multModelScore(
            this,
            posStrandSeq,
            positionalWeights,
            strand,
            calc,
            eToMu,
            nonSpecKa,
            //this.getRevCompSimilarity(),
            revCompSimilarity,
            withholdFeatureKey);


        // if the multModelScore is negative!!
        if (score < 0) {
            System.out.println("Error: MultModelScore is negative!!!");
            if (nonSpecKa != 0) {
                score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
            }
            else {
                if (mandatoryColumns == null) {
                    //mandatoryColumns = new boolean[2][this.columns];
                    mandatoryColumns = new boolean[2][];
                }
                // score with PSAM only
                score = WeightMatrixTools.score(
                    posStrandPWM,
                    negStrandPWM,
                    posStrandSeq,
                    positionalWeights,
                    mandatoryColumns,
                    strand,
                    calc,
                    eToMu,
                    nonSpecKa,
                    revCompSimilarity);

            }
        }

        return(score);
    }

    // score WITHOUT weights
    // does not use the mandatory columns of each feature
    public double multModelScore(SymbolList posStrandSeq, FeatureKey withholdFeatureKey, int startPos) {


        //double score = WeightMatrixTools.score(
        double score = multModelScore(
            this,
            posStrandSeq,
            strand,
            calc,
            eToMu,
            nonSpecKa,
            //this.getRevCompSimilarity(),
            revCompSimilarity,
            startPos,
            withholdFeatureKey);


        // if the multModelScore is negative!!
        if (score < 0) {
            System.out.println("Error: MultModelScore is negative!!!");
            if (nonSpecKa != 0) {
                score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
            }
            else {
                // if (mandatoryColumns == null) {
                //     //mandatoryColumns = new boolean[2][this.columns];
                //     mandatoryColumns = new boolean[2][];
                // }

                // score with PSAM only
                score = WeightMatrixTools.score(
                    posStrandPWM,
                    negStrandPWM,
                    posStrandSeq,
                    strand,
                    calc,
                    eToMu,
                    nonSpecKa,
                    revCompSimilarity,
                    startPos);

            }
        }

        return(score);
    }

    // score WITHOUT weights
    // does not use the mandatory columns of each feature
    public double multModelScore(SymbolList posStrandSeq, ScoreType scoreType) {
        return(multModelScore(posStrandSeq, (FeatureKey)null));
    }

    public double multModelScore(SymbolList posStrandSeq) {
        return(multModelScore(posStrandSeq, (FeatureKey)null));
    }

    public double multModelScore(SymbolList posStrandSeq, FeatureKey withholdFeatureKey) {

        //double score = WeightMatrixTools.score(
        double score = multModelScore(
            this,
            posStrandSeq,
            strand,
            calc,
            eToMu,
            nonSpecKa,
            //this.getRevCompSimilarity(),
            revCompSimilarity,
            withholdFeatureKey);


        // if the multModelScore is negative!!
        if (score < 0) {
            System.out.println("Error: MultModelScore is negative!!!");
            if (nonSpecKa != 0) {
                score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
            }
            else {
                // if (mandatoryColumns == null) {
                //     //mandatoryColumns = new boolean[2][this.columns];
                //     mandatoryColumns = new boolean[2][];
                // }

                // score with PSAM only
                score = WeightMatrixTools.score(
                    posStrandPWM,
                    negStrandPWM,
                    posStrandSeq,
                    strand,
                    calc,
                    eToMu,
                    nonSpecKa,
                    revCompSimilarity);

            }
        }

        return(score);
    }


    // No weights, with startPos, one strand
    public double multModelScore(
        FeaturedWeightMatrix anFWM,
        SymbolList aSeq,
        WeightMatrixTools.BindingStrand strand,
        int startPos,
        boolean[] mandatoryCols,
        FeatureKey withholdFeatureKey) {

        double logScore = 0;
        try {
            WeightMatrix pwm = anFWM.getPWM(strand);

            // set the startColumn
            int startColumn = 0;
            if (startPos < 0) {
                startColumn = -1 * startPos;
            }

            // set the endColumn
            int seqLength = aSeq.length();
            int endColumn = anFWM.columns();
            if (startPos + endColumn > seqLength) {
                endColumn = seqLength - startPos;
            }

            if (mandatoryCols != null) {
                // if mandatory column skipped, then return 0
                for (int c = 0; c < startColumn; c++) {
                    if (mandatoryCols[c]) {
                        return(0.0);
                    }
                }
                for (int c = endColumn; c < pwm.columns(); c++) {
                    if (mandatoryCols[c]) {
                        return(0.0);
                    }
                }
            }

            // calculate score
            for (int c = startColumn; c < endColumn; c++) {
                double prob = pwm.getColumn(c).getWeight(aSeq.symbolAt(startPos+c+1));
                if (prob == 0.0) {
                    return(0.0);
                }
                logScore += Math.log(prob);
            }

            //System.out.println("Hi 1");

            if (featureKeyToMultFeatureMap != null) {
                // iterate through all the features found in this window
                //ArrayList<WeightMatrixFeature> nNDDFeatures = anFWM.getNNDDFeatures(aSeq, strand, startColumn, endColumn);
                //ArrayList<WeightMatrixFeature> nNDDFeatures = anFWM.getNNDDFeatures(aSeq, strand, startPos+startColumn, endColumn);

                ArrayList<WeightMatrixFeature> features = new ArrayList<WeightMatrixFeature>();

                features.addAll(anFWM.getNNDDFeatures(aSeq, strand, startPos, startColumn, endColumn));

                //features.addAll(anFWM.getNNNDDFeatures(aSeq, strand, startPos, startColumn, endColumn));

//                 features.addAll(anFWM.getNNTDFeatures(aSeq, strand, startPos, startColumn, endColumn));

                //System.out.println("Hi 2");

                for (WeightMatrixFeature feature : features) {
                    //boolean[] mandatoryColumns = feature.getMandatoryColumns();

                    //System.out.println("Hi 3");

                    // if should withhold this feature then continue to next feature
                    if ((withholdFeatureKey!=null) && (withholdFeatureKey==feature.getKey())) {
                        continue;
                    }

                    double featureRelAffinity = feature.getRelativeAffinity();

                    //System.out.println(feature.toString()+"\n");
                    //System.out.println("featureRelAffinity="+featureRelAffinity+"\n");

//                     SymbolTokenization parser = pwm.getAlphabet().getTokenization("default");
//                     System.out.println("window="+aSeq.subStr(startPos+startColumn+1, startPos+endColumn)+"; feature="+StringTools.toString(feature.getModsArray(), parser, "", "-")+"; strand="+strand+"; seq="+aSeq.seqString()+"\n");

                    if (featureRelAffinity == 0.0) {
                        return(0.0);
                    }
                    logScore += Math.log(featureRelAffinity);
                }
            }

        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return(Math.exp(logScore));
    }

    // With weights, with startPos, both strands
    // weights are applied after occupancy calculation
    public double multModelScore(
        FeaturedWeightMatrix anFWM,
        SymbolList posStrandSeq,
        double[][] weights,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double eToMu,
        double nonSpecKa,
        double revCompSimilarity,
        int startPos,
        FeatureKey withholdFeatureKey) {

        double posScore = 0;
        double negScore = 0;
        int weightsIndex = startPos + anFWM.columns() - 1;
        double assocKa;

        if (mandatoryColumns == null) {
            mandatoryColumns = new boolean[2][];
        }

        // return the right score
        switch (strand) {

        case POS:
            // score the window with the anFWM and posStrand
            assocKa = multModelScore(anFWM, posStrandSeq, strand, startPos, mandatoryColumns[0], withholdFeatureKey);
            if (eToMu == 0) {
                posScore = assocKa + nonSpecKa;
            }
            else {
                double prodKa = eToMu * (assocKa + nonSpecKa);
                posScore = prodKa / (prodKa + 1);
            }
            posScore = weights[0][weightsIndex] * posScore;
            return(posScore);

        case NEG:
            // score the window with the negStrandPWM
            assocKa = multModelScore(anFWM, posStrandSeq, strand, startPos, mandatoryColumns[1], withholdFeatureKey);
            if (eToMu == 0) {
                negScore = assocKa + nonSpecKa;
            }
            else {
                double prodKa = eToMu * (assocKa + nonSpecKa);
                negScore = prodKa / (prodKa + 1);
            }
            negScore = weights[1][weightsIndex] * negScore;
            return(negScore);

        case BOTH:
            double posAssocKa = multModelScore(anFWM, posStrandSeq, WeightMatrixTools.BindingStrand.POS, startPos, mandatoryColumns[0], withholdFeatureKey);
            double negAssocKa = multModelScore(anFWM, posStrandSeq, WeightMatrixTools.BindingStrand.NEG, startPos, mandatoryColumns[1], withholdFeatureKey);
            if (eToMu == 0) {
                posScore = posAssocKa + nonSpecKa;
                negScore = negAssocKa + nonSpecKa;
            }
            else {
                double posProdKa = eToMu * (posAssocKa + nonSpecKa);
                double negProdKa = eToMu * (negAssocKa + nonSpecKa);
                posScore = posProdKa / (posProdKa + 1);
                negScore = negProdKa / (negProdKa + 1);
            }
            posScore = weights[0][weightsIndex] * posScore;
            negScore = weights[1][weightsIndex] * negScore;

            switch (calc) {
            case MAX:
                return(Math.max(posScore, negScore));
            case SUM:
                return(posScore + negScore);
            case UNION:
                return(posScore + negScore - (posScore*negScore));
            case NORMED_SUM:
                return((posScore + negScore) / (revCompSimilarity + 1));
            }
            break;

        }
        return(0);
    }

    // No weights, with startPos, both strands
    public double multModelScore(
        FeaturedWeightMatrix anFWM,
        SymbolList posStrandSeq,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double eToMu,
        double nonSpecKa,
        double revCompSimilarity,
        int startPos,
        FeatureKey withholdFeatureKey) {

        double posScore = 0;
        double negScore = 0;
        double assocKa;

        boolean[] emptyMandatoryColumns = null;

        // return the right score
        switch (strand) {

        case POS:
            // score the window with the anFWM and posStrand
            assocKa = multModelScore(anFWM, posStrandSeq, strand, startPos, emptyMandatoryColumns, withholdFeatureKey);
            if (eToMu == 0) {
                posScore = assocKa + nonSpecKa;
            }
            else {
                double prodKa = eToMu * (assocKa + nonSpecKa);
                posScore = prodKa / (prodKa + 1);
            }
            return(posScore);

        case NEG:
            // score the window with the negStrandPWM
            assocKa = multModelScore(anFWM, posStrandSeq, strand, startPos, emptyMandatoryColumns, withholdFeatureKey);
            if (eToMu == 0) {
                negScore = assocKa + nonSpecKa;
            }
            else {
                double prodKa = eToMu * (assocKa + nonSpecKa);
                negScore = prodKa / (prodKa + 1);
            }
            return(negScore);

        case BOTH:
            double posAssocKa = multModelScore(anFWM, posStrandSeq, WeightMatrixTools.BindingStrand.POS, startPos, emptyMandatoryColumns, withholdFeatureKey);
            double negAssocKa = multModelScore(anFWM, posStrandSeq, WeightMatrixTools.BindingStrand.NEG, startPos, emptyMandatoryColumns, withholdFeatureKey);
            if (eToMu == 0) {
                posScore = posAssocKa + nonSpecKa;
                negScore = negAssocKa + nonSpecKa;
            }
            else {
                double posProdKa = eToMu * (posAssocKa + nonSpecKa);
                double negProdKa = eToMu * (negAssocKa + nonSpecKa);
                posScore = posProdKa / (posProdKa + 1);
                negScore = negProdKa / (negProdKa + 1);
            }

            switch (calc) {
            case MAX:
                return(Math.max(posScore, negScore));
            case SUM:
                return(posScore + negScore);
            case UNION:
                return(posScore + negScore - (posScore*negScore));
            case NORMED_SUM:
                return((posScore + negScore) / (revCompSimilarity + 1));
            }
            break;

        }
        return(0);
    }

    // sliding window
    public double multModelScore(
        FeaturedWeightMatrix anFWM,
        SymbolList posStrandSeq,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double eToMu,
        double nonSpecKa,
        double revCompSimilarity,
        FeatureKey withholdFeatureKey)
    {
        double scoreSum = 0;
        // Sweep down the posStrandSeq scoring each window
        // 0-based indexing
        for (int startPos = 0; startPos <= posStrandSeq.length() - anFWM.columns(); startPos++) {
            scoreSum += multModelScore(anFWM, posStrandSeq, strand, calc, eToMu, nonSpecKa, revCompSimilarity, startPos, withholdFeatureKey);
        }
        return(scoreSum);
    }

    // sliding window
    public double[] multModelScoreWindows(
        FeaturedWeightMatrix anFWM,
        SymbolList posStrandSeq,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double eToMu,
        double nonSpecKa,
        double revCompSimilarity,
        FeatureKey withholdFeatureKey
        // double minThreshold
                                        )
    {
        double[] scores = new double[posStrandSeq.length() - anFWM.columns() + 1];

        // Sweep down the posStrandSeq scoring each window
        // 0-based indexing
        for (int startPos = 0; startPos <= posStrandSeq.length() - anFWM.columns(); startPos++) {
            // scores[startPos] = multModelScore(anFWM, posStrandSeq, strand, calc, eToMu, nonSpecKa, revCompSimilarity, startPos, withholdFeatureKey, minThreshold);
            scores[startPos] = multModelScore(anFWM, posStrandSeq, strand, calc, eToMu, nonSpecKa, revCompSimilarity, startPos, withholdFeatureKey);
        }
        return(scores);
    }

    // sliding window
    public double multModelScore(
        FeaturedWeightMatrix anFWM,
        SymbolList posStrandSeq,
        double[][] weights,
        WeightMatrixTools.BindingStrand strand,
        WeightMatrixTools.BothStrandsCalc calc,
        double eToMu,
        double nonSpecKa,
        double revCompSimilarity,
        FeatureKey withholdFeatureKey)
    {
        double scoreSum = 0;
        // Sweep down the posStrandSeq scoring each window
        // 0-based indexing
        for (int startPos = -1 * (anFWM.columns() - 1); startPos < posStrandSeq.length(); startPos++) {
            scoreSum += multModelScore(anFWM, posStrandSeq, weights, strand, calc, eToMu, nonSpecKa, revCompSimilarity, startPos, withholdFeatureKey);
        }
        return(scoreSum);
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Additive Model
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    public SymbolList getHighestAffinitySymList(LinkedHashSet<SymbolList> hammingSphere, boolean additive) {

        SymbolList highestAffinitySymList = null;
        double highestAffinity = 0;

        for (SymbolList aSymList : hammingSphere) {

            double aScore;
            if (additive) {
                // get Additive Model Score(featureSeq)
                aScore = addModelScore(
                    aSymList,
                    //this.positionalWeights,
                    (ScoreType)null);
            }
            else {
                // get Multiplicative Model Score(featureSeq)
                aScore = multModelScore(
                    aSymList,
                    //this.positionalWeights,
                    (ScoreType)null);
            }

            if (aScore > highestAffinity) {
                highestAffinitySymList = aSymList;
                highestAffinity = aScore;
            }
        }

        return(highestAffinitySymList);
    }

    public double getR2(LinkedHashSet<SymbolList> hammingSphere) {
        int length = hammingSphere.size();
        double[] addScores = new double[length];
        double[] multScores = new double[length];

        int i = 0;
        for (SymbolList aSymList : hammingSphere) {

            // get Additive Model Score(featureSeq)
            addScores[i] = addModelScore(
                aSymList,
                //this.positionalWeights,
                (ScoreType)null);

            // get current Mult Model Score(featureSeq)
            multScores[i] = multModelScore(
                aSymList,
                //this.positionalWeights,
                (ScoreType)null);

            i++;

        }
        return(MathTools.rSquared(addScores, multScores));
    }

    public double getR2psam(LinkedHashSet<SymbolList> hammingSphere) {
        return(getR2psam(hammingSphere, this.posStrandPWM, this.negStrandPWM));
    }

    public double getR2psam(
        LinkedHashSet<SymbolList> hammingSphere,
        WeightMatrix aPosStrandPWM,
        WeightMatrix aNegStrandPWM)
    {
        int length = hammingSphere.size();
        double[] addScores = new double[length];
        double[] psamScores = new double[length];

        int i = 0;
        for (SymbolList aSymList : hammingSphere) {

            // get Additive Model Score(featureSeq)
            addScores[i] = addModelScore(
                aSymList,
                //this.positionalWeights,
                (ScoreType)null);

            // get new PSAM Model Score(featureSeq)
            psamScores[i] = WeightMatrixTools.score(
                aPosStrandPWM,
                aNegStrandPWM,
                aSymList,
                this.strand,
                this.calc,
                this.eToMu, //eToMu
                this.nonSpecKa, //nonSpecKa
                this.revCompSimilarity);

            i++;

        }
        return(MathTools.rSquared(addScores, psamScores));
    }

    public void equilibrate() {

        try {

            sorted = false;

            // set the oldAMPSAMs
            // if (this.oldAMposStrandPWM == null) {
                this.oldAMposStrandPWM = WeightMatrixTools.copy(this.posStrandPWM);
                this.oldAMnegStrandPWM = WeightMatrixTools.copy(this.negStrandPWM);
            // }

            WeightMatrixTools.setMinAffinity(this.posStrandPWM, .001);
            WeightMatrixTools.setMinAffinity(this.negStrandPWM, .001);


            // remove mult features
            featureKeyToMultFeatureMap = new LinkedHashMap<FeatureKey, WeightMatrixFeature>();


            // get all the symbol lists, and the hamming Sphere around the highest affinity seq
            // System.out.print("\nGenerating all sequences of length "+this.columns+"...");
            // LinkedHashSet<SymbolList> allSymLists = SymbolListTools.getAllSymbolLists(this.alphabet, this.columns);
            // System.out.println("Done.");

            // System.out.print("\nScoring all the "+this.columns+"-bp sequences with the FSAM (additive model)...");
            // SymbolList seedSymList = getHighestAffinitySymList(allSymLists, true);
            // System.out.println("Done.");

            // SymbolList seedSymList = WeightMatrixTools.getMostLikelySymList(this.posStrandPWM);
            SymbolList seedSymList = DNATools.createDNA("ATGATTAATTAC");

            SymbolList revComp = DNATools.reverseComplement(seedSymList);
            double seedScore = psamScore(seedSymList, WeightMatrixTools.BindingStrand.POS);
            double revCompScore = psamScore(revComp, WeightMatrixTools.BindingStrand.POS);
            if (revCompScore > seedScore) {
                seedSymList = revComp;
            }

            //SymbolList seedSymList = WeightMatrixTools.getMostLikelySymList(this.posStrandPWM);
            //SymbolList seedSymList = WeightMatrixTools.getMostLikelySymList(this.oldAMposStrandPWM);

            LinkedHashSet<SymbolList> hammingSphere = null;
            hammingSphere = SymbolListTools.getHammingSphere(seedSymList, 2, true, DNATools.complementTable());

            // get new PSAM based on additive feature model
            int iters = 2;
            double[] r2s = new double[iters];
            double maxR2 = 0;

            System.out.println("\nCalculating rSquareds between the Additive and PSAMs (over Hamming-2 Sphere):");
            for (int i = 0; i < iters; i++) {

                /////////////////////////////////////////////////////////////////////////////////
                //WeightMatrix newPosStrandPWM = equilibratePSAM(seedSymList, true);
                WeightMatrix newPosStrandPWM = this.posStrandPWM;
                /////////////////////////////////////////////////////////////////////////////////

                WeightMatrix newNegStrandPWM = WeightMatrixTools.reverseComplement(newPosStrandPWM, DNATools.complementTable());

                r2s[i] = getR2psam(hammingSphere, newPosStrandPWM, newNegStrandPWM);

                if (r2s[i] > maxR2) {
                    // set the fsam.psam to the new PSAM
                    this.posStrandPWM = newPosStrandPWM;
                    this.negStrandPWM = newNegStrandPWM;
                    maxR2 = r2s[i];
                }

                System.out.println("\t"+  (i+1) +": "+r2s[i]);
            }
            System.out.println("\nrSquareds = "+StringTools.toString(r2s, ", ")+"\n");

            if (this.featureKeyToMultFeatureMap == null) {
                this.featureKeyToMultFeatureMap = new LinkedHashMap<FeatureKey, WeightMatrixFeature>();
            }

            /////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////

            // calc multiplicative nnDD features repeatably in order to converge to additive model.
            // seedSymList = WeightMatrixTools.getMostLikelySymList(this.posStrandPWM);
            seedSymList = DNATools.createDNA("ATGATTAATTAC");

            hammingSphere = SymbolListTools.getHammingSphere(seedSymList, 2, true, DNATools.complementTable());

            //iters = 1;
            iters = 20;
            r2s = new double[iters];
            System.out.println("\nCalculating rSquareds between the Additive and Multiplicative Models (over Hamming-2 Sphere):");
            for (int i = 0; i < iters; i++) {

                ////////////////////////////////////////
                // nnDD
                ////////////////////////////////////////
                // get new multiplicative (nnDD) features BEFORE resetting the fsam.psam to the new PSAM
                // LinkedHashMap<FeatureKey, WeightMatrixFeature> newFeatureKeyToMultFeatureMap = equilibrateMultNNDD();
                LinkedHashMap<FeatureKey, WeightMatrixFeature> newFeatureKeyToMultFeatureMap = equilibrateMultNNDD(seedSymList);

                // set the multiplicative fsam.NNDD features to the new nnDDs
                // this.featureKeyToMultFeatureMap = newFeatureKeyToMultFeatureMap;


                ////////////////////////////////////////
                // nnTD
                ////////////////////////////////////////
                // get new multiplicative (nnDD) features BEFORE resetting the fsam.psam to the new PSAM
                //             newFeatureKeyToMultFeatureMap = equilibrateMultNNTD();

                // set the multiplicative fsam.NNDD features to the new nnDDs
                // this.featureKeyToMultFeatureMap = newFeatureKeyToMultFeatureMap;

                ///////////////////////////////////////////////////////////
                // get r2 between Additive and Multiplicative models
                ///////////////////////////////////////////////////////////
                r2s[i] = getR2(hammingSphere);
                System.out.println("\t"+  (i+1) +": "+r2s[i]);
            }
            System.out.println("\nrSquareds = "+StringTools.toString(r2s, ", ")+"\n");

        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

    }

    public void equilibrateMultModel(LinkedHashSet<SymbolList> symLists) {

        try {

            System.out.println("\nRe-equilibrating the PSAM based on the multiplicative nnDDs...");

            System.out.print("\n\tScoring all the "+this.columns+"-bp sequences with the FSAM (multiplicative model)...");
            SymbolList seedSymList = getHighestAffinitySymList(symLists, false); //multModel
            System.out.println("Done.");

            SymbolList revComp = DNATools.reverseComplement(seedSymList);
            double seedScore = psamScore(seedSymList, WeightMatrixTools.BindingStrand.POS);
            double revCompScore = psamScore(revComp, WeightMatrixTools.BindingStrand.POS);
            if (revCompScore > seedScore) {
                seedSymList = revComp;
            }

            // Re-equilibrate with latest nnDDs to Create a new PSAM
            /////////////////////////////////////////////////////////////////////////////////
            // WeightMatrix newPSAM = equilibratePSAM(seedSymList, false); // multModel
            // newPSAM = WeightMatrixTools.weightedAverage(.7, .3, posStrandPWM, newPSAM);

            WeightMatrix newPSAM = this.posStrandPWM;
            /////////////////////////////////////////////////////////////////////////////////

            // Re-equilibrate new nnDDs
            LinkedHashMap<FeatureKey, WeightMatrixFeature> newMultFeatures = equilibrateMultNNDD(newPSAM);

            // Re-set this FSAM with the new PSAM and MultFeatures
            putPWMs(newPSAM, DNATools.complementTable());
            this.featureKeyToMultFeatureMap = newMultFeatures;

            System.out.println("Done.");
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

    }

    // create a new PSAM by performing HammingReduce on the Additive Model
    public WeightMatrix equilibratePSAM(SymbolList seedSymList, boolean additive) {
        WeightMatrix newPosStrandPWM = null;

        try {

            this.nonSpecKa = 0;

            newPosStrandPWM = new SimpleWeightMatrix(this.alphabet, this.columns, DistributionFactory.DEFAULT);
            newPosStrandPWM.setName("PSAM");

            //SymbolList seedSymList = WeightMatrixTools.getMostLikelySymList(posStrandPWM);
            //SymbolList seedSymList = WeightMatrixTools.getMostLikelySymList(this.oldAMposStrandPWM);

            int endColumn;
            if (this.isRevCompPalindrome) {
                endColumn = this.columns / 2;
            }
            else {
                endColumn = this.columns;
            }

            // 1-based indexing
            for (int i = 1; i <= endColumn; i++) {
            // for (int i = 1; i <= posStrandPWM.columns(); i++) {

                boolean[] freezeArray = new boolean[this.columns];
                freezeArray[i-1] = true;

                int hammingDistance = 0;
                if (additive) {
                    hammingDistance = 1;
                }

                LinkedHashSet<SymbolList> hammingSphere = SymbolListTools.getHammingSphere(
                    seedSymList,
                    hammingDistance,
                    true, // include reverse complements
                    freezeArray,
                    DNATools.complementTable());

                LinkedHashMap<Symbol, ArrayList<Double>> ratiosHM = new LinkedHashMap<Symbol, ArrayList<Double>>();
                for (SymbolList hammingSphereSymList : hammingSphere) {

                    ArrayList<Symbol> symbols = new ArrayList<Symbol>(4);
                    //ArrayList<Double> weights = new ArrayList<Double>(4);
                    List<Number> weights = new ArrayList<Number>(4);

                    // iterate over every symbol at this position
                    for (Iterator dnaIter = ((FiniteAlphabet)alphabet).iterator(); dnaIter.hasNext(); ) {
                        Symbol symbol = (Symbol)dnaIter.next();
                        SimpleSymbolList editedSeed = new SimpleSymbolList(hammingSphereSymList);
                        Edit edit = new Edit(i, alphabet, symbol);
                        editedSeed.edit(edit);

                        double score;
                        if (additive) {
                            score = addModelScore(
                                editedSeed,
                                //this.positionalWeights,
                                (ScoreType)null);
                        }
                        else {
                            score = multModelScore(
                                editedSeed,
                                //this.positionalWeights,
                                (ScoreType)null);
                        }

                        symbols.add(symbol);
                        weights.add(score);
                    }

                    double maxScore = MathTools.max(weights);
                    for (int j=0; j < weights.size(); j++) {
                        double newRatio = weights.get(j).doubleValue() / maxScore;
                        Symbol symbol = symbols.get(j);
                        ArrayList<Double> ratiosArray = ratiosHM.get(symbol);
                        if (ratiosArray == null) {
                            ratiosArray = new ArrayList<Double>();
                            ratiosHM.put(symbol, ratiosArray);
                        }
                        ratiosArray.add(newRatio);
//                         Double oldSum = ratiosHM.get(symbol);
//                         if (oldSum == null) {
//                             oldSum = new Double(0);
//                         }
//                         ratiosHM.put(symbol, oldSum + newRatio);
                    }
                }

                // calculate trimmed mean for all the relative affinities across the hamming sphere
                LinkedHashMap<Symbol, Double> weightsHM = new LinkedHashMap<Symbol, Double>();
                for (Symbol symbol : ratiosHM.keySet()) {
                    double estimator;
                    if (additive) {
                        estimator = MathTools.trimmedMean(ArrayTools.toDoubleArray(ratiosHM.get(symbol)), .20);
                    }
                    else {
                        //System.out.println("coefficients length = "+ArrayTools.toDoubleArray(ratiosHM.get(symbol)).length);
                        estimator = MathTools.mean(ArrayTools.toDoubleArray(ratiosHM.get(symbol)));
                    }
                    weightsHM.put(symbol, estimator);
//                     double avgRatio = ratiosHM.get(symbol) / hammingSphere.size();
//                     newPosStrandPWM.getColumn(i-1).setWeight(symbol, avgRatio);
                }

                // re-normalize affinities
                //double maxWeight = MathTools.max(weightsHM.values().toArray(new Double[]));
                double maxWeight = MathTools.max(weightsHM.values());
                for (Symbol symbol : ratiosHM.keySet()) {
                    double newWeight = weightsHM.get(symbol) / maxWeight;
                    newPosStrandPWM.getColumn(i-1).setWeight(symbol, newWeight);
                }


            }

            // if revCompPalidrome then need to copy the reverseComplement of the
            // second half of the newPSAM
            if (isRevCompPalindrome) {
                for (int destPos = endColumn; destPos < this.columns; destPos++) {
                    int sourcePos = (this.columns - 1) - destPos;
                    System.out.println("Copying complement of column "+(sourcePos+1)+" distribution to column "+(destPos+1)+".");
                    //dists[destPos] = DistributionTools.complement(dists[sourcePos], complementTable);
                    newPosStrandPWM.setColumn(destPos, DistributionTools.complement(newPosStrandPWM.getColumn(sourcePos), DNATools.complementTable()));
                }
            }


        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(newPosStrandPWM);
    }

//     // creates a new multiplicative model (on top of the currentMult) based on the AdditiveModel, and currentMultiplicativeModel
//     public LinkedHashMap<FeatureKey, WeightMatrixFeature> equilibrateMultNNDD() {
//         LinkedHashMap<FeatureKey, WeightMatrixFeature> newFeatureMap = new LinkedHashMap<FeatureKey, WeightMatrixFeature>();

//         try {
//             this.nonSpecKa = 0;

//             SymbolList seedSymList = WeightMatrixTools.getMostLikelySymList(this.posStrandPWM);

//             // get Additive Model Score(HighestAffinitySeq)
//             double seedAddScore = addModelScore(
//                 seedSymList,
//                 //this.positionalWeights,
//                 (ScoreType) null);


//             // 0-based indexing
//             //for (int mut1 = 0; mut1 < 5; mut1++) {
//             for (int mut1 = 0; mut1 < this.columns-1; mut1++) {
//                 for (int mut2 = mut1 + 1; mut2 <= mut1 + 1; mut2++) { // Nearest-Neighbor

//                     // iterate over every symbol at the mut1 position
//                     for (Iterator dnaIter1 = ((FiniteAlphabet)alphabet).iterator(); dnaIter1.hasNext(); ) {
//                         Symbol mut1Symbol = (Symbol)dnaIter1.next();

//                         // iterate over every symbol at the mut2 position
//                         for (Iterator dnaIter2 = ((FiniteAlphabet)alphabet).iterator(); dnaIter2.hasNext(); ) {
//                             Symbol mut2Symbol = (Symbol)dnaIter2.next();

// //                             Symbol avoidMut1Symbol = DistributionTools.getMostLikelySymbol(this.posStrandPWM.getColumn(mut1));
// //                             Symbol avoidMut2Symbol = DistributionTools.getMostLikelySymbol(this.posStrandPWM.getColumn(mut2));

// //                             if ( (mut1Symbol == avoidMut1Symbol) && (mut2Symbol == avoidMut2Symbol) ) {
// //                                 //out.println((byte)100, "\nNot including high affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression.");
// //                                 continue;
// //                             }

//                             SimpleSymbolList editedSeed = new SimpleSymbolList(seedSymList);
//                             // Edit uses 1-based indexing
//                             Edit edit1 = new Edit(mut1+1, alphabet, mut1Symbol);
//                             editedSeed.edit(edit1);
//                             Edit edit2 = new Edit(mut2+1, alphabet, mut2Symbol);
//                             editedSeed.edit(edit2);

//                             Symbol[] modsArray = new Symbol[(this.columns*2)-1];
//                             modsArray[(mut1)*2] = mut1Symbol;
//                             modsArray[(mut2)*2] = mut2Symbol;

//                             FeatureKey featureKey = new FeatureKey(modsArray);
//                             FeatureKey featureKeyRevComp = featureKey.reverseComplement(DNATools.complementTable());

//                             boolean[] mandatoryColumns = new boolean[this.columns];
//                             mandatoryColumns[mut1] = true;
//                             mandatoryColumns[mut2] = true;

//                             double newFeatureScore = 1;
//                             if ((strand == WeightMatrixTools.BindingStrand.BOTH)
//                                 && (this.isRevCompPalindrome)
//                                 && (newFeatureMap.containsKey(featureKeyRevComp))) {

//                                 newFeatureScore = newFeatureMap.get(featureKeyRevComp).getRelativeAffinity();
//                             }
//                             // else calculate the multiplicative Feature Score
//                             else {
//                                 LinkedHashSet<SymbolList> hammingSphere = SymbolListTools.getHammingSphere(
//                                     editedSeed,
//                                     0, //hamming distance
//                                     true, // include reverse complements
//                                     mandatoryColumns,
//                                     DNATools.complementTable());

//                                 ArrayList<Double> ratiosArray = new ArrayList<Double>();
//                                 for (SymbolList hammingSphereSymList : hammingSphere) {

//                                     /////////////////////////////////////////////////////////////
//                                     // get Additive Model Score(featureSeq)
//                                     /////////////////////////////////////////////////////////////
//                                     double addScore = addModelScore(
//                                         hammingSphereSymList,
//                                         //this.positionalWeights,
//                                         (ScoreType)null);

//                                     /////////////////////////////////////////////////////////////
//                                     // get current Mult Model Score(featureSeq) WHILE WITHHOLDING THIS FEATURE
//                                     /////////////////////////////////////////////////////////////
//                                     //double currentMultScore = score(
//                                     double currentMultScore = multModelScore(
//                                         hammingSphereSymList,
//                                         //this.positionalWeights,
//                                         featureKey);

//                                     double newRatio = (addScore/seedAddScore) / currentMultScore;
//                                     ratiosArray.add(newRatio);
//                                     //ratioSum += (addScore/seedAddScore) / currentMultScore;

//                                 }
//                                 newFeatureScore = MathTools.trimmedMean(ArrayTools.toDoubleArray(ratiosArray), .20);
//                                 //double newFeatureScore = ratioSum / hammingSphere.size();
//                                 //double newFeatureScore = (addScore/currentMultScore) / seedAddScore;
//                                 //System.out.println("\naddScore="+addScore+"\tcurrentMultScore="+currentMultScore+"\tseedAddScore="+seedAddScore+"\n");
//                             }

//                             if (newFeatureScore == 1.0) {
//                                 continue;
//                             }

//                             WeightMatrix featurePSAM = WeightMatrixTools.mutate(this.posStrandPWM, mut1, mut1Symbol, true);
//                             featurePSAM = WeightMatrixTools.mutate(featurePSAM, mut2, mut2Symbol, false);

//                             // name has 1-based indexing
//                             featurePSAM.setName(
//                                 "mutate_"+ (mut1+1) +"_"+ mut1Symbol.getName().substring(0,1).toUpperCase()
//                                 +"__"+ "mutate_"+ (mut2+1) +"_"+ mut2Symbol.getName().substring(0,1).toUpperCase() +"");

//                             String label = "dependency";

//                             boolean isGoodFeature;
//                             if (newFeatureScore >= 1.0) {
//                                 isGoodFeature = true;
//                             }
//                             else {
//                                 isGoodFeature = false;
//                             }

//                             WeightMatrixFeature feature = new WeightMatrixFeature(
//                                 label,
//                                 featurePSAM,
//                                 mandatoryColumns,
//                                 isGoodFeature,
//                                 false,
//                                 newFeatureScore,
//                                 modsArray,
//                                 DNATools.complementTable());

//                             newFeatureMap.put(feature.getKey(), feature);

//                             ////////////////////////////////////////////////////////////////////////////////////
//                             // put it into the current feature map if greedy
//                             ////////////////////////////////////////////////////////////////////////////////////
//                             featureKeyToMultFeatureMap.put(feature.getKey(), feature);

//                         }
//                     }
//                 }
//             }
//         }
//         catch (Exception ex) {
//             ex.printStackTrace();
//         }
//         return(newFeatureMap);
//     }


    // creates a new multiplicative model (on top of the currentMult) based on the new PSAM, and currentMultiplicativeModel
    public LinkedHashMap<FeatureKey, WeightMatrixFeature> equilibrateMultNNDD(WeightMatrix aNewPSAM) {
        LinkedHashMap<FeatureKey, WeightMatrixFeature> newFeatureMap = new LinkedHashMap<FeatureKey, WeightMatrixFeature>();

        try {

            for (WeightMatrixFeature feature : featureKeyToMultFeatureMap.values()) {
                Symbol[] modsArray = feature.getModsArray();
                FeatureKey newFeatureKey = new FeatureKey(modsArray);
                feature.setKey(newFeatureKey);

                double oldPosIndependScore = WeightMatrixTools.scoreModsArray(this.posStrandPWM, modsArray);
                double newPosIndependScore = WeightMatrixTools.scoreModsArray(aNewPSAM, modsArray);

                double newFeatureAffinity = (feature.getRelativeAffinity() * oldPosIndependScore) / newPosIndependScore;

                boolean isGoodFeature;
                if (newFeatureAffinity >= 1.0) {
                    isGoodFeature = true;
                }
                else {
                    isGoodFeature = false;
                }

                WeightMatrixFeature newFeature = new WeightMatrixFeature(
                    feature.getLabel(),
                    feature.getPosStrandPWM(),
                    feature.getPosStrandMandatoryColumns(),
                    isGoodFeature,
                    false,
                    newFeatureAffinity,
                    modsArray,
                    DNATools.complementTable());


                newFeatureMap.put(newFeatureKey, newFeature);
            }

        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(newFeatureMap);
    }

    // creates a new multiplicative model (on top of the currentMult) based on the AdditiveModel, and current MultiplicativeModel
    public LinkedHashMap<FeatureKey, WeightMatrixFeature> equilibrateMultNNDD(SymbolList seedSymList) {
        LinkedHashMap<FeatureKey, WeightMatrixFeature> newFeatureMap = new LinkedHashMap<FeatureKey, WeightMatrixFeature>();

        try {

            this.nonSpecKa = 0;

            //SymbolList seedSymList = WeightMatrixTools.getMostLikelySymList(this.posStrandPWM);

            // get Additive Model Score(HighestAffinitySeq)
            double seedAddScore = addModelScore(
                seedSymList,
                //this.positionalWeights,
                (ScoreType) null);


            // 0-based indexing
            //for (int mut1 = 0; mut1 < 5; mut1++) {
            for (int mut1 = 0; mut1 < this.columns-1; mut1++) {
                for (int mut2 = mut1 + 1; mut2 <= mut1 + 1; mut2++) { // Nearest-Neighbor

                    // iterate over every symbol at the mut1 position
                    for (Iterator dnaIter1 = ((FiniteAlphabet)alphabet).iterator(); dnaIter1.hasNext(); ) {
                        Symbol mut1Symbol = (Symbol)dnaIter1.next();

                        // iterate over every symbol at the mut2 position
                        for (Iterator dnaIter2 = ((FiniteAlphabet)alphabet).iterator(); dnaIter2.hasNext(); ) {
                            Symbol mut2Symbol = (Symbol)dnaIter2.next();

//                             Symbol avoidMut1Symbol = DistributionTools.getMostLikelySymbol(this.posStrandPWM.getColumn(mut1));
//                             Symbol avoidMut2Symbol = DistributionTools.getMostLikelySymbol(this.posStrandPWM.getColumn(mut2));

//                             if ( (mut1Symbol == avoidMut1Symbol) && (mut2Symbol == avoidMut2Symbol) ) {
//                                 //out.println((byte)100, "\nNot including high affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression.");
//                                 continue;
//                             }

                            SimpleSymbolList editedSeed = new SimpleSymbolList(seedSymList);
                            // Edit uses 1-based indexing
                            Edit edit1 = new Edit(mut1+1, alphabet, mut1Symbol);
                            editedSeed.edit(edit1);
                            Edit edit2 = new Edit(mut2+1, alphabet, mut2Symbol);
                            editedSeed.edit(edit2);

                            Symbol[] modsArray = new Symbol[(this.columns*2)-1];
                            modsArray[(mut1)*2] = mut1Symbol;
                            modsArray[(mut2)*2] = mut2Symbol;

                            FeatureKey featureKey = new FeatureKey(modsArray);
                            FeatureKey featureKeyRevComp = featureKey.reverseComplement(DNATools.complementTable());

                            boolean[] mandatoryColumns = new boolean[this.columns];
                            mandatoryColumns[mut1] = true;
                            mandatoryColumns[mut2] = true;

                            double newFeatureScore = 1;
                            if ((strand == WeightMatrixTools.BindingStrand.BOTH)
                                && (this.isRevCompPalindrome)
                                && (newFeatureMap.containsKey(featureKeyRevComp))) {

                                newFeatureScore = newFeatureMap.get(featureKeyRevComp).getRelativeAffinity();
                            }
                            // else calculate the multiplicative Feature Score
                            else {
                                LinkedHashSet<SymbolList> hammingSphere = SymbolListTools.getHammingSphere(
                                    editedSeed,
                                    1, //hamming distance
                                    true, // include reverse complements
                                    mandatoryColumns, // freeze where the mutations are made
                                    DNATools.complementTable());

                                ArrayList<Double> ratiosArray = new ArrayList<Double>();
                                for (SymbolList hammingSphereSymList : hammingSphere) {

                                    /////////////////////////////////////////////////////////////
                                    // get Additive Model Score(featureSeq)
                                    /////////////////////////////////////////////////////////////
                                    double addScore = addModelScore(
                                        hammingSphereSymList,
                                        //this.positionalWeights,
                                        (ScoreType)null);

                                    /////////////////////////////////////////////////////////////
                                    // get current Mult Model Score(featureSeq) WHILE WITHHOLDING THIS FEATURE
                                    /////////////////////////////////////////////////////////////
                                    double currentMultScore = multModelScore(
                                        hammingSphereSymList,
                                        //this.positionalWeights,
                                        featureKey);

                                    /////////////////////////////////////////////////////////////
                                    // get current seed Mult Model Score(featureSeq) WHILE WITHHOLDING THIS FEATURE
                                    /////////////////////////////////////////////////////////////
                                    double currentSeedMultScore = multModelScore(
                                        seedSymList,
                                        //this.positionalWeights,
                                        featureKey);

                                    double newRatio = (addScore/seedAddScore) / (currentMultScore/currentSeedMultScore);
                                    //double newRatio = (addScore/seedAddScore) / currentMultScore;
                                    //double newRatio = addScore/seedAddScore;

                                    ratiosArray.add(newRatio);
                                    //ratioSum += (addScore/seedAddScore) / currentMultScore;

                                    // if ((mut1==7) && (mut2==8)) {
                                    //     System.out.print("\n"+hammingSphereSymList.seqString().toUpperCase()+": addScore="+addScore+"\tcurrentMultScore="+currentMultScore+"\tseedAddScore="+seedAddScore+"\tcurrentSeedMultScore="+currentSeedMultScore+"\tnewRatio="+newRatio);
                                    // }
                                }

                                newFeatureScore = MathTools.trimmedMean(ArrayTools.toDoubleArray(ratiosArray), .20);
                                //newFeatureScore = MathTools.mean(ArrayTools.toDoubleArray(ratiosArray));

                                //double newFeatureScore = ratioSum / hammingSphere.size();
                                //double newFeatureScore = (addScore/currentMultScore) / seedAddScore;
                                //System.out.println("\nAfter trimmed mean: addScore="+addScore+"\tcurrentMultScore="+currentMultScore+"\tseedAddScore="+seedAddScore+"\n");


                                ////////////////////////////////////////////////////////////////////////////////
                                // put a cap on the max featureScore
                                ////////////////////////////////////////////////////////////////////////////////
                                // if (newFeatureScore > 10000)
                                // {
                                //     newFeatureScore = 10000;
                                // }
                                // if (newFeatureScore > 1000)
                                // {
                                //     newFeatureScore = 1000;
                                // }
                                // if (newFeatureScore > 20)
                                // {
                                //     newFeatureScore = 20;
                                // }

                                ////////////////////////////////////////////////////////////////////////////////
                                // if there already are mult features then perform a step in this direction
                                ////////////////////////////////////////////////////////////////////////////////
                                if (this.featureKeyToMultFeatureMap != null)
                                {
                                    WeightMatrixFeature oldFeature = featureKeyToMultFeatureMap.get(featureKey);

                                    double oldFeatureWeight = 1;
                                    if (oldFeature != null) {
                                        oldFeatureWeight = oldFeature.getRelativeAffinity();
                                    }
                                    else {
                                        System.out.println("\nMultiplicative NNDD feature not found!");
                                    }

                                    //double weightedAvg = ((9*oldFeatureWeight) + (1*newFeatureScore)) / 10;
                                    double weightedAvg = ((4*oldFeatureWeight) + (1*newFeatureScore)) / 5;
                                    //double weightedAvg = ((1*oldFeatureWeight) + (1*newFeatureScore)) / 2;
                                    //double weightedAvg = ((1*oldFeatureWeight) + (2*newFeatureScore)) / 3;
                                    //double weightedAvg = ((1*oldFeatureWeight) + (4*newFeatureScore)) / 5;
                                    newFeatureScore = weightedAvg;

                                }
                            }

                            if (newFeatureScore == 1.0) {
                                continue;
                            }

                            WeightMatrix featurePSAM = WeightMatrixTools.mutate(this.posStrandPWM, mut1, mut1Symbol, true);
                            featurePSAM = WeightMatrixTools.mutate(featurePSAM, mut2, mut2Symbol, false);

                            // name has 1-based indexing
                            featurePSAM.setName(
                                "mutate_"+ (mut1+1) +"_"+ mut1Symbol.getName().substring(0,1).toUpperCase()
                                +"__"+ "mutate_"+ (mut2+1) +"_"+ mut2Symbol.getName().substring(0,1).toUpperCase() +"");

                            String label = "dependency";

                            boolean isGoodFeature;
                            if (newFeatureScore >= 1.0) {
                                isGoodFeature = true;
                            }
                            else {
                                isGoodFeature = false;
                            }

                            WeightMatrixFeature feature = new WeightMatrixFeature(
                                label,
                                featurePSAM,
                                mandatoryColumns,
                                isGoodFeature,
                                false,
                                newFeatureScore,
                                modsArray,
                                DNATools.complementTable());

                            newFeatureMap.put(feature.getKey(), feature);

                            ////////////////////////////////////////////////////////////////////////////////////
                            // put it into the current feature map if greedy
                            ////////////////////////////////////////////////////////////////////////////////////
                            featureKeyToMultFeatureMap.put(feature.getKey(), feature);

                        }
                    }
                }
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(newFeatureMap);
    }

    // creates a new multiplicative model (on top of the currentMult) based on the AdditiveModel, and currentMultiplicativeModel
    public LinkedHashMap<FeatureKey, WeightMatrixFeature> equilibrateMultNNTD() {
        LinkedHashMap<FeatureKey, WeightMatrixFeature> newFeatureMap = new LinkedHashMap<FeatureKey, WeightMatrixFeature>();

        try {
            SymbolList seedSymList = WeightMatrixTools.getMostLikelySymList(this.posStrandPWM);

            // get Additive Model Score(HighestAffinitySeq)
            double seedAddScore = addModelScore(
                seedSymList,
                //this.positionalWeights,
                (ScoreType) null);

            // 0-based indexing
            //for (int mut1 = 0; mut1 < 5; mut1++) {
            for (int mut1 = 0; mut1 < this.columns-2; mut1++) {
                for (int mut2 = mut1 + 1; mut2 <= mut1 + 1; mut2++) { // Nearest-Neighbor
                    for (int mut3 = mut2 + 1; mut3 <= mut2 + 1; mut3++) { // Nearest-Neighbor

                        // iterate over every symbol at the mut1 position
                        for (Iterator dnaIter1 = ((FiniteAlphabet)alphabet).iterator(); dnaIter1.hasNext(); ) {
                            Symbol mut1Symbol = (Symbol)dnaIter1.next();

                            // iterate over every symbol at the mut2 position
                            for (Iterator dnaIter2 = ((FiniteAlphabet)alphabet).iterator(); dnaIter2.hasNext(); ) {
                                Symbol mut2Symbol = (Symbol)dnaIter2.next();

                                // iterate over every symbol at the mut3 position
                                for (Iterator dnaIter3 = ((FiniteAlphabet)alphabet).iterator(); dnaIter3.hasNext(); ) {
                                    Symbol mut3Symbol = (Symbol)dnaIter3.next();

                                    //                             Symbol avoidMut1Symbol = DistributionTools.getMostLikelySymbol(this.posStrandPWM.getColumn(mut1));
                                    //                             Symbol avoidMut2Symbol = DistributionTools.getMostLikelySymbol(this.posStrandPWM.getColumn(mut2));

                                    //                             if ( (mut1Symbol == avoidMut1Symbol) && (mut2Symbol == avoidMut2Symbol) ) {
                                    //                                 //out.println((byte)100, "\nNot including high affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression.");
                                    //                                 continue;
                                    //                             }

                                    SimpleSymbolList editedSeed = new SimpleSymbolList(seedSymList);
                                    // Edit uses 1-based indexing
                                    Edit edit1 = new Edit(mut1+1, alphabet, mut1Symbol);
                                    editedSeed.edit(edit1);
                                    Edit edit2 = new Edit(mut2+1, alphabet, mut2Symbol);
                                    editedSeed.edit(edit2);
                                    Edit edit3 = new Edit(mut3+1, alphabet, mut3Symbol);
                                    editedSeed.edit(edit3);

                                    Symbol[] modsArray = new Symbol[(this.columns*2)-1];
                                    modsArray[(mut1)*2] = mut1Symbol;
                                    modsArray[(mut2)*2] = mut2Symbol;
                                    modsArray[(mut3)*2] = mut3Symbol;

                                    FeatureKey featureKey = new FeatureKey(modsArray);
                                    FeatureKey featureKeyRevComp = featureKey.reverseComplement(DNATools.complementTable());

                                    boolean[] mandatoryColumns = new boolean[this.columns];
                                    mandatoryColumns[mut1] = true;
                                    mandatoryColumns[mut2] = true;
                                    mandatoryColumns[mut3] = true;

                                    double newFeatureScore = 1;
                                    if ((strand == WeightMatrixTools.BindingStrand.BOTH)
                                        && (this.isRevCompPalindrome)
                                        && (newFeatureMap.containsKey(featureKeyRevComp))) {

                                        newFeatureScore = newFeatureMap.get(featureKeyRevComp).getRelativeAffinity();
                                    }
                                    // else calculate the multiplicative Feature Score
                                    else {
                                        LinkedHashSet<SymbolList> hammingSphere = SymbolListTools.getHammingSphere(
                                            editedSeed,
                                            1, //hamming distance
                                            true, // include reverse complements
                                            mandatoryColumns,
                                            DNATools.complementTable());

                                        ArrayList<Double> ratiosArray = new ArrayList<Double>();
                                        for (SymbolList hammingSphereSymList : hammingSphere) {

                                            /////////////////////////////////////////////////////////////
                                            // get Additive Model Score(featureSeq)
                                            /////////////////////////////////////////////////////////////
                                            double addScore = addModelScore(
                                                hammingSphereSymList,
                                                //this.positionalWeights,
                                                (ScoreType)null);

                                            /////////////////////////////////////////////////////////////
                                            // get current Mult Model Score(featureSeq) WHILE WITHHOLDING THIS FEATURE
                                            /////////////////////////////////////////////////////////////
                                            double currentMultScore = multModelScore(
                                                hammingSphereSymList,
                                                //this.positionalWeights,
                                                featureKey);

                                            double newRatio = (addScore/seedAddScore) / currentMultScore;
                                            ratiosArray.add(newRatio);
                                            //ratioSum += (addScore/seedAddScore) / currentMultScore;

                                        }

                                        newFeatureScore = MathTools.trimmedMean(ArrayTools.toDoubleArray(ratiosArray), .20);
                                        //double newFeatureScore = ratioSum / hammingSphere.size();
                                        //double newFeatureScore = (addScore/currentMultScore) / seedAddScore;
                                        //System.out.println("\naddScore="+addScore+"\tcurrentMultScore="+currentMultScore+"\tseedAddScore="+seedAddScore+"\n");
                                    }

                                    if (newFeatureScore == 1.0) {
                                        continue;
                                    }

                                    WeightMatrix featurePSAM = WeightMatrixTools.mutate(this.posStrandPWM, mut1, mut1Symbol, true);
                                    featurePSAM = WeightMatrixTools.mutate(featurePSAM, mut2, mut2Symbol, false);

                                    // name has 1-based indexing
                                    featurePSAM.setName(
                                        "mutate_"+ (mut1+1) +"_"+ mut1Symbol.getName().substring(0,1).toUpperCase()
                                        +"__"+ "mutate_"+ (mut2+1) +"_"+ mut2Symbol.getName().substring(0,1).toUpperCase()
                                        +"__"+ "mutate_"+ (mut3+1) +"_"+ mut3Symbol.getName().substring(0,1).toUpperCase() +"");

                                    String label = "dependency";

                                    boolean isGoodFeature;
                                    if (newFeatureScore >= 1.0) {
                                        isGoodFeature = true;
                                    }
                                    else {
                                        isGoodFeature = false;
                                    }

                                    WeightMatrixFeature feature = new WeightMatrixFeature(
                                        label,
                                        featurePSAM,
                                        mandatoryColumns,
                                        isGoodFeature,
                                        false,
                                        newFeatureScore,
                                        modsArray,
                                        DNATools.complementTable());

                                    newFeatureMap.put(feature.getKey(), feature);

                                    ////////////////////////////////////////////////////////////////////////////////////
                                    // put it into the current feature map if greedy
                                    ////////////////////////////////////////////////////////////////////////////////////
                                    featureKeyToMultFeatureMap.put(feature.getKey(), feature);

                                }
                            }
                        }
                    }
                }
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(newFeatureMap);
    }

    public SymbolList getModSymbolList(SymbolList seedSymList, Symbol[] modsArray) {
        Symbol[] modSymArray = new Symbol[seedSymList.length()];
        for (int i = 0; i < seedSymList.length(); i++) {
            if (modsArray[i*2] == null) {
                modSymArray[i] = seedSymList.symbolAt(i+1);
            }
            else {
                modSymArray[i] = modsArray[i*2];
            }
        }
        return(new SimpleSymbolList(modSymArray, modSymArray.length, this.alphabet));
    }

// //     // creates a new multiplicative model (on top of the currentMult) based on the AdditiveModel, and currentMultiplicativeModel
// //     public LinkedHashMap<FeatureKey, WeightMatrixFeature> equilibrateMultNNDD() {
// //         LinkedHashMap<FeatureKey, WeightMatrixFeature> newFeatureMap = new LinkedHashMap<FeatureKey, WeightMatrixFeature>();

// //         try {
// //             SymbolList seedSymList = WeightMatrixTools.getMostLikelySymList(this.posStrandPWM);

// //             // get Additive Model Score(HighestAffinitySeq)
// //             double seedAddScore = addModelScore(
// //                 seedSymList,
// //                 //this.positionalWeights,
// //                 (ScoreType) null);

// //             for (WeightMatrixFeature feature : featureKeyToAddFeatureMap.values()) {


// // //                 FeatureKey newFeatureKey = new FeatureKey(feature.getModsArray());
// // //                 feature.setKey(newFeatureKey);
// // //                 newFeatureKeyToFeatureMap.put(newFeatureKey, feature);
// //                 WeightMatrixFeature featureClone = feature.clone();
// //                 SymbolList modSymbolList = getModSymbolList(seedSymList, featureClone.getModsArray());

// //                 // get Additive Model Score(featureSeq)
// //                 double addScore = addModelScore(
// //                     modSymbolList,
// //                     //this.positionalWeights,
// //                     (ScoreType)null);

// //                 // get current Mult Model Score(featureSeq) WHILE WITHHOLDING THIS FEATURE
// //                 double currentMultScore = multModelScore(
// //                     modSymbolList,
// //                     //this.positionalWeights,
// //                     featureClone.getKey());

// //                 double newFeatureScore = (addScore/currentMultScore) / seedAddScore;
// //                 //newFeatureScore = Math.pow(newFeatureScore, 1/3);
// //                 //newFeatureScore = Math.pow(newFeatureScore, 1/2);
// //                 //newFeatureScore = Math.pow(newFeatureScore, 0.8);
// // //                 newFeatureScore = Math.pow(newFeatureScore, 0.9);
// // //                 newFeatureScore *= 0.1;

// //                 //System.out.println("\naddScore="+addScore+"\tcurrentMultScore="+currentMultScore+"\tseedAddScore="+seedAddScore+"\n");

// //                 featureClone.setRelativeAffinity(newFeatureScore);

// //                 boolean isGoodFeature;
// //                 if (newFeatureScore >= 1.0) {
// //                     isGoodFeature = true;
// //                 }
// //                 else {
// //                     isGoodFeature = false;
// //                 }
// //                 featureClone.setIsGood(isGoodFeature);

// //                 newFeatureMap.put(featureClone.getKey(), featureClone);
// // //                 if (newFeatureScore < 1000) {
// // //                     newFeatureMap.put(featureClone.getKey(), featureClone);
// // //                 }
// // //                 else {
// // //                     System.out.println("Rejected Feature Score "+newFeatureScore+" for feature "+modSymbolList.seqString());
// // //                 }

// //                 ////////////////////////////////////////////////////////////////////////////////////
// //                 // put it into the current feature map if greedy
// //                 ////////////////////////////////////////////////////////////////////////////////////
// //                 featureKeyToMultFeatureMap.put(featureClone.getKey(), featureClone);

// //             }
// //         }
// //         catch (Exception ex) {
// //             ex.printStackTrace();
// //         }
// //         return(newFeatureMap);
// //     }

    // creates a new additive model (on top of the newPSAM) based on the oldPSAM, oldAdditiveModel, and newPSAM
    public LinkedHashMap<FeatureKey, WeightMatrixFeature> equilibrateAddNNDD(WeightMatrix aNewPosStrandPWM) {
        LinkedHashMap<FeatureKey, WeightMatrixFeature> newFeatureMap = new LinkedHashMap<FeatureKey, WeightMatrixFeature>();

        try {
            for (WeightMatrixFeature feature : featureKeyToAddFeatureMap.values()) {


//                 FeatureKey newFeatureKey = new FeatureKey(feature.getModsArray());
//                 feature.setKey(newFeatureKey);
//                 newFeatureKeyToFeatureMap.put(newFeatureKey, feature);
                WeightMatrixFeature featureClone = feature.clone();

                double oldPSAMproduct = getPSAMProduct(this.posStrandPWM, featureClone.getModsArray());
                double newPSAMproduct = getPSAMProduct(aNewPosStrandPWM, featureClone.getModsArray());
                featureClone.setRelativeAffinity(featureClone.getRelativeAffinity() + oldPSAMproduct - newPSAMproduct);
//                 featureClone.setRelativeAffinity(featureClone.getRelativeAffinity() - oldPSAMproduct + newPSAMproduct);

                newFeatureMap.put(featureClone.getKey(), featureClone);
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(newFeatureMap);
    }

    public double getPSAMProduct(WeightMatrix aPSAM, Symbol[] modsArray) {
        double product = 1.0;
        try {
            for (int i = 0; i < aPSAM.columns(); i++) {
                if (modsArray[i*2] != null) {
                    product *= aPSAM.getColumn(i).getWeight(modsArray[i*2]);
                }
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return(product);
    }

//     // creates a new multiplicative model (on top of the newPSAM) based on the oldPSAM, oldAdditiveModel, and newPSAM
//     public LinkedHashMap<FeatureKey, WeightMatrixFeature> equilibrateMultNNDD(WeightMatrix aNewPosStrandPWM) {
//         LinkedHashMap<FeatureKey, WeightMatrixFeature> newFeatureMap = new LinkedHashMap<FeatureKey, WeightMatrixFeature>();

//         try {
//             SymbolList seedSymList = WeightMatrixTools.getMostLikelySymList(aNewPosStrandPWM);
//             WeightMatrix aNewPosStrandPWMRevComp = WeightMatrixTools.reverseComplement(aNewPosStrandPWM, DNATools.complementTable());
//             double newRevCompSimilarity = WeightMatrixTools.similarity(aNewPosStrandPWM, aNewPosStrandPWMRevComp);

//             // get Additive Model Score(HighestAffinitySeq)
//             double seedAddScore = addModelScore(
//                 seedSymList,
//                 //this.positionalWeights,
//                 (ScoreType) null);

//             // 0-based indexing
//             for (int mut1 = 0; mut1 < aNewPosStrandPWM.columns()-1; mut1++) {
//                 for (int mut2 = mut1 + 1; mut2 <= mut1 + 1; mut2++) { // Nearest-Neighbor

//                     // iterate over every symbol at the mut1 position
//                     for (Iterator dnaIter1 = ((FiniteAlphabet)alphabet).iterator(); dnaIter1.hasNext(); ) {
//                         Symbol mut1Symbol = (Symbol)dnaIter1.next();

//                         // iterate over every symbol at the mut2 position
//                         for (Iterator dnaIter2 = ((FiniteAlphabet)alphabet).iterator(); dnaIter2.hasNext(); ) {
//                             Symbol mut2Symbol = (Symbol)dnaIter2.next();

// //                             Symbol avoidMut1Symbol = DistributionTools.getMostLikelySymbol(aNewPosStrandPWM.getColumn(mut1));
// //                             Symbol avoidMut2Symbol = DistributionTools.getMostLikelySymbol(aNewPosStrandPWM.getColumn(mut2));

// //                             if ( (mut1Symbol == avoidMut1Symbol) && (mut2Symbol == avoidMut2Symbol) ) {
// //                                 //out.println((byte)100, "\nNot including high affinity dinucleotide (("+(mut1+1)+","+mut1Symbol.getName().substring(0,1).toUpperCase()+"), ("+(mut2+1)+","+mut2Symbol.getName().substring(0,1).toUpperCase()+")) in regression.");
// //                                 continue;
// //                             }

//                             SimpleSymbolList editedSeed = new SimpleSymbolList(seedSymList);
//                             // Edit uses 1-based indexing
//                             Edit edit1 = new Edit(mut1+1, alphabet, mut1Symbol);
//                             editedSeed.edit(edit1);
//                             Edit edit2 = new Edit(mut2+1, alphabet, mut2Symbol);
//                             editedSeed.edit(edit2);

//                             // get Additive Model Score(featureSeq)
//                             double addScore = addModelScore(
//                                 editedSeed,
//                                 //this.positionalWeights,
//                                 (ScoreType)null);

//                             // get new PSAM Model Score(featureSeq)
//                             double newPsamScore = WeightMatrixTools.score(
//                                 aNewPosStrandPWM,
//                                 aNewPosStrandPWMRevComp,
//                                 editedSeed,
//                                 this.strand,
//                                 this.calc,
//                                 this.eToMu, //eToMu
//                                 this.nonSpecKa, //nonSpecKa
//                                 newRevCompSimilarity);

//                             double newFeatureScore = (addScore/newPsamScore) / seedAddScore;
//                             //System.out.println("\naddScore="+addScore+"\tnewPsamScore="+newPsamScore+"\tseedAddScore="+seedAddScore+"\n");

//                             if (newFeatureScore == 1.0) {
//                                 continue;
//                             }

//                             WeightMatrix featurePSAM = WeightMatrixTools.mutate(aNewPosStrandPWM, mut1, mut1Symbol, true);
//                             featurePSAM = WeightMatrixTools.mutate(featurePSAM, mut2, mut2Symbol, false);

//                             // name has 1-based indexing
//                             featurePSAM.setName(
//                                 "mutate_"+ (mut1+1) +"_"+ mut1Symbol.getName().substring(0,1).toUpperCase()
//                                 +"__"+ "mutate_"+ (mut2+1) +"_"+ mut2Symbol.getName().substring(0,1).toUpperCase() +"");

//                             String label = "dependency";

//                             boolean isGoodFeature;
//                             if (newFeatureScore >= 1.0) {
//                                 isGoodFeature = true;
//                             }
//                             else {
//                                 isGoodFeature = false;
//                             }

//                             boolean[] mandatoryColumns = new boolean[aNewPosStrandPWM.columns()];
//                             mandatoryColumns[mut1] = true;
//                             mandatoryColumns[mut2] = true;

//                             Symbol[] modsArray = new Symbol[(aNewPosStrandPWM.columns()*2)-1];
//                             modsArray[(mut1)*2] = mut1Symbol;
//                             modsArray[(mut2)*2] = mut2Symbol;

//                             WeightMatrixFeature feature = new WeightMatrixFeature(
//                                 label,
//                                 featurePSAM,
//                                 mandatoryColumns,
//                                 isGoodFeature,
//                                 false,
//                                 newFeatureScore,
//                                 modsArray,
//                                 DNATools.complementTable());

//                             newFeatureMap.put(feature.getKey(), feature);
//                         }
//                     }
//                 }
//             }
//         }
//         catch (Exception ex) {
//             ex.printStackTrace();
//         }
//         return(newFeatureMap);
//     }

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// New AddModelScores
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

    // // score WITH the weights
    // // use the mandatory columns of each feature
    // public double addModelScore(SymbolList posStrandSeq, double[][] weights, FeatureKey withholdFeatureKey) {

    //     this.positionalWeights = weights;

    //     //double score = WeightMatrixTools.score(
    //     double score = addModelScore(
    //         posStrandSeq,
    //         positionalWeights,
    //         strand,
    //         calc,
    //         eToMu,
    //         nonSpecKa,
    //         //this.getRevCompSimilarity(),
    //         revCompSimilarity,
    //         withholdFeatureKey);


    //     if (score < 0) {
    //         if (nonSpecKa != 0) {
    //             score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
    //         }
    //         else {
    //             if (mandatoryColumns == null) {
    //                 //mandatoryColumns = new boolean[2][this.columns];
    //                 mandatoryColumns = new boolean[2][];
    //             }
    //             // score with PSAM only
    //             score = WeightMatrixTools.score(
    //                 posStrandPWM,
    //                 negStrandPWM,
    //                 posStrandSeq,
    //                 positionalWeights,
    //                 mandatoryColumns,
    //                 strand,
    //                 calc,
    //                 eToMu,
    //                 nonSpecKa,
    //                 revCompSimilarity);

    //         }
    //     }

    //     return(score);
    // }

    // // score WITHOUT weights
    // // does not use the mandatory columns of each feature
    // public double addModelScore(SymbolList posStrandSeq, FeatureKey withholdFeatureKey, int startPos) {

    //     //double score = WeightMatrixTools.score(
    //     double score = addModelScore(
    //         posStrandSeq,
    //         strand,
    //         calc,
    //         eToMu,
    //         nonSpecKa,
    //         //this.getRevCompSimilarity(),
    //         revCompSimilarity,
    //         startPos,
    //         withholdFeatureKey);


    //     if (score < 0) {
    //         if (nonSpecKa != 0) {
    //             score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
    //         }
    //         else {
    //             if (mandatoryColumns == null) {
    //                 //mandatoryColumns = new boolean[2][this.columns];
    //                 mandatoryColumns = new boolean[2][];
    //             }
    //             // score with PSAM only
    //             score = WeightMatrixTools.score(
    //                 posStrandPWM,
    //                 negStrandPWM,
    //                 posStrandSeq,
    //                 strand,
    //                 calc,
    //                 eToMu,
    //                 nonSpecKa,
    //                 revCompSimilarity,
    //                 startPos);

    //         }
    //     }

    //     return(score);
    // }

    // // score WITHOUT weights
    // // does not use the mandatory columns of each feature
    // public double addModelScore(SymbolList posStrandSeq, ScoreType scoreType) {
    //     return(addModelScore(posStrandSeq, (FeatureKey)null));
    // }

    // // score WITHOUT weights
    // // does not use the mandatory columns of each feature
    // public double addModelScore(SymbolList posStrandSeq, FeatureKey withholdFeatureKey) {


    //     //double score = WeightMatrixTools.score(
    //     double score = addModelScore(
    //         posStrandSeq,
    //         strand,
    //         calc,
    //         eToMu,
    //         nonSpecKa,
    //         //this.getRevCompSimilarity(),
    //         revCompSimilarity,
    //         withholdFeatureKey);


    //     if (score < 0) {
    //         if (nonSpecKa != 0) {
    //             score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
    //         }
    //         else {
    //             if (mandatoryColumns == null) {
    //                 //mandatoryColumns = new boolean[2][this.columns];
    //                 mandatoryColumns = new boolean[2][];
    //             }
    //             // score with PSAM only
    //             score = WeightMatrixTools.score(
    //                 posStrandPWM,
    //                 negStrandPWM,
    //                 posStrandSeq,
    //                 strand,
    //                 calc,
    //                 eToMu,
    //                 nonSpecKa,
    //                 revCompSimilarity);

    //         }
    //     }

    //     return(score);
    // }


    // public double addModelScore(
    //     //FeaturedWeightMatrix anFWM,
    //     SymbolList posStrandSeq,
    //     WeightMatrixTools.BindingStrand strand,
    //     WeightMatrixTools.BothStrandsCalc calc,
    //     double eToMu,
    //     double nonSpecKa,
    //     double revCompSimilarity,
    //     int startPos,
    //     FeatureKey withholdFeatureKey)
    // {
    //     double scoreSum = 0;
    //     double psamScore = 0;

    //     WeightMatrix addPosStrandPWM;
    //     WeightMatrix addNegStrandPWM;
    //     if (oldAMposStrandPWM == null) {
    //         addPosStrandPWM = posStrandPWM;
    //         addNegStrandPWM = negStrandPWM;
    //     }
    //     else {
    //         addPosStrandPWM = oldAMposStrandPWM;
    //         addNegStrandPWM = oldAMnegStrandPWM;
    //     }

    //     if (includePwm) {
    //         scoreSum += WeightMatrixTools.score(
    //             addPosStrandPWM,
    //             addNegStrandPWM,
    //             posStrandSeq,
    //             strand,
    //             calc,
    //             eToMu,
    //             nonSpecKa,
    //             revCompSimilarity,
    //             startPos);
    //         psamScore = scoreSum;
    //     }

    //     // add any features
    //     if (featureKeyToMultFeatureMap != null) {

    //         // set the startColumn
    //         int startColumn = 0;
    //         if (startPos < 0) {
    //             startColumn = -1 * startPos;
    //         }

    //         // set the endColumn
    //         int seqLength = posStrandSeq.length();
    //         int endColumn = this.columns;
    //         if (startPos + endColumn > seqLength) {
    //             endColumn = seqLength - startPos;
    //         }

    //         ArrayList<WeightMatrixFeature> features = new ArrayList<WeightMatrixFeature>();
    //         features.addAll(getNNDDFeatures(posStrandSeq, strand, startPos, startColumn, endColumn));

    //         for (WeightMatrixFeature feature : features) {

    //             if ((withholdFeatureKey!=null) && (withholdFeatureKey==feature.getKey())) {
    //                 continue;
    //             }

    //             double featureRelAffinity = feature.getRelativeAffinity();

    //             double featureScore = WeightMatrixTools.score(
    //                 feature.getPosStrandPWM(),
    //                 feature.getNegStrandPWM(),
    //                 posStrandSeq,
    //                 strand,
    //                 calc,
    //                 eToMu,
    //                 0, //nonSpecKa,
    //                 feature.getRevCompSimilarity(),
    //                 startPos);

    //             scoreSum += featureRelAffinity * featureScore;
    //         }
    //     }

    //     if (scoreSum < 0) {
    //         if (nonSpecKa != 0) {
    //             scoreSum = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
    //         }
    //         else {
    //             scoreSum = psamScore;
    //         }
    //     }

    //     return(scoreSum);
    // }

    // public double addModelScore(
    //     //FeaturedWeightMatrix anFWM,
    //     SymbolList posStrandSeq,
    //     WeightMatrixTools.BindingStrand strand,
    //     WeightMatrixTools.BothStrandsCalc calc,
    //     double eToMu,
    //     double nonSpecKa,
    //     double revCompSimilarity,
    //     FeatureKey withholdFeatureKey)
    // {
    //     double scoreSum = 0;
    //     double psamScore = 0;

    //     WeightMatrix addPosStrandPWM;
    //     WeightMatrix addNegStrandPWM;
    //     if (oldAMposStrandPWM == null) {
    //         addPosStrandPWM = posStrandPWM;
    //         addNegStrandPWM = negStrandPWM;
    //     }
    //     else {
    //         addPosStrandPWM = oldAMposStrandPWM;
    //         addNegStrandPWM = oldAMnegStrandPWM;
    //     }

    //     // Sweep down the posStrandSeq scoring each window
    //     // 0-based indexing
    //     for (int startPos = 0; startPos <= posStrandSeq.length() - this.columns; startPos++) {

    //         if (includePwm) {
    //             scoreSum += WeightMatrixTools.score(
    //                 addPosStrandPWM,
    //                 addNegStrandPWM,
    //                 posStrandSeq,
    //                 strand,
    //                 calc,
    //                 eToMu,
    //                 nonSpecKa,
    //                 revCompSimilarity,
    //                 startPos);
    //             psamScore = scoreSum;
    //         }

    //         // add any features
    //         if (featureKeyToMultFeatureMap != null) {

    //             // set the startColumn
    //             int startColumn = 0;
    //             if (startPos < 0) {
    //                 startColumn = -1 * startPos;
    //             }

    //             // set the endColumn
    //             int seqLength = posStrandSeq.length();
    //             int endColumn = this.columns;
    //             if (startPos + endColumn > seqLength) {
    //                 endColumn = seqLength - startPos;
    //             }

    //             ArrayList<WeightMatrixFeature> features = new ArrayList<WeightMatrixFeature>();
    //             features.addAll(getNNDDFeatures(posStrandSeq, strand, startPos, startColumn, endColumn));

    //             for (WeightMatrixFeature feature : features) {

    //                 if ((withholdFeatureKey!=null) && (withholdFeatureKey==feature.getKey())) {
    //                     continue;
    //                 }

    //                 double featureRelAffinity = feature.getRelativeAffinity();

    //                 double featureScore = WeightMatrixTools.score(
    //                     feature.getPosStrandPWM(),
    //                     feature.getNegStrandPWM(),
    //                     posStrandSeq,
    //                     strand,
    //                     calc,
    //                     eToMu,
    //                     0, //nonSpecKa,
    //                     feature.getRevCompSimilarity(),
    //                     startPos);

    //                 scoreSum += featureRelAffinity * featureScore;
    //             }
    //         }
    //     }

    //     if (scoreSum < 0) {
    //         if (nonSpecKa != 0) {
    //             scoreSum = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
    //         }
    //         else {
    //             scoreSum = psamScore;
    //         }
    //     }

    //     return(scoreSum);
    // }

    // public double addModelScore(
    //     //FeaturedWeightMatrix anFWM,
    //     SymbolList posStrandSeq,
    //     double[][] weights,
    //     WeightMatrixTools.BindingStrand strand,
    //     WeightMatrixTools.BothStrandsCalc calc,
    //     double eToMu,
    //     double nonSpecKa,
    //     double revCompSimilarity,
    //     FeatureKey withholdFeatureKey)
    // {
    //     double scoreSum = 0;
    //     double psamScore = 0;

    //     WeightMatrix addPosStrandPWM;
    //     WeightMatrix addNegStrandPWM;
    //     if (oldAMposStrandPWM == null) {
    //         addPosStrandPWM = posStrandPWM;
    //         addNegStrandPWM = negStrandPWM;
    //     }
    //     else {
    //         addPosStrandPWM = oldAMposStrandPWM;
    //         addNegStrandPWM = oldAMnegStrandPWM;
    //     }

    //     // create an empty mandatoryColumns
    //     boolean[][] emptyMandatoryColumns = new boolean[2][this.columns];

    //     // Sweep down the posStrandSeq scoring each window
    //     // 0-based indexing
    //     for (int startPos = -1 * (this.columns - 1); startPos < posStrandSeq.length(); startPos++) {

    //         if (includePwm) {
    //             scoreSum += WeightMatrixTools.score(
    //                 addPosStrandPWM,
    //                 addNegStrandPWM,
    //                 posStrandSeq,
    //                 weights,
    //                 emptyMandatoryColumns,
    //                 strand,
    //                 calc,
    //                 eToMu,
    //                 nonSpecKa,
    //                 revCompSimilarity,
    //                 startPos);
    //             psamScore = scoreSum;
    //         }

    //         // add any features
    //         if (featureKeyToMultFeatureMap != null) {

    //             // set the startColumn
    //             int startColumn = 0;
    //             if (startPos < 0) {
    //                 startColumn = -1 * startPos;
    //             }

    //             // set the endColumn
    //             int seqLength = posStrandSeq.length();
    //             int endColumn = this.columns;
    //             if (startPos + endColumn > seqLength) {
    //                 endColumn = seqLength - startPos;
    //             }

    //             ArrayList<WeightMatrixFeature> features = new ArrayList<WeightMatrixFeature>();
    //             features.addAll(getNNDDFeatures(posStrandSeq, strand, startPos, startColumn, endColumn));

    //             for (WeightMatrixFeature feature : features) {

    //                 if ((withholdFeatureKey!=null) && (withholdFeatureKey==feature.getKey())) {
    //                     continue;
    //                 }

    //                 double featureRelAffinity = feature.getRelativeAffinity();

    //                 double featureScore = WeightMatrixTools.score(
    //                     feature.getPosStrandPWM(),
    //                     feature.getNegStrandPWM(),
    //                     posStrandSeq,
    //                     weights,
    //                     feature.getMandatoryColumns(),
    //                     strand,
    //                     calc,
    //                     eToMu,
    //                     0, //nonSpecKa
    //                     feature.getRevCompSimilarity(),
    //                     startPos);

    //                 scoreSum += featureRelAffinity * featureScore;
    //             }
    //         }
    //     }

    //     if (scoreSum < 0) {
    //         if (nonSpecKa != 0) {
    //             scoreSum = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
    //         }
    //         else {
    //             scoreSum = psamScore;
    //         }
    //     }

    //     return(scoreSum);
    // }



/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// old AddModelScores
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////
    // score WITH the weights
    // use the mandatory columns of each feature
    //////////////////////////////////////////////
    public double addModelScore(SymbolList posStrandSeq, double[][] weights, FeatureKey blah) {
        return(addModelScore(posStrandSeq, weights, null, 0));
    }

    public double addModelScore(SymbolList posStrandSeq, double[][] weights, ScoreType scoreType) {
        return(addModelScore(posStrandSeq, weights, scoreType, 0));
    }

    public double addModelScore(SymbolList posStrandSeq, double[][] weights, ScoreType scoreType, double threshold) {

        this.positionalWeights = weights;

        //System.out.println("\nthreshold = "+threshold);
        //System.out.println("\nweights = "+StringTools.toString(weights, "\n", "\t"));

        // If the PSAM is an reverse-complement palindrome, then we should look at only the positive
        // or negative strand (BUT NOT BOTH)
//         if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.BOTH)) {
//             //out.println("\n\nSince the motif is self reverse-complement, we should look only at the positive strand.");
//             strand = WeightMatrixTools.BindingStrand.POS;
//         }

        // create an empty mandatoryColumns
        boolean[][] mandatoryColumns = new boolean[2][this.columns];

        // score the posStrandSeq with the posStrandPWM
        double score = 0;
        double psamScore = 0;

        WeightMatrix addPosStrandPWM;
        WeightMatrix addNegStrandPWM;
        if (oldAMposStrandPWM == null) {
            addPosStrandPWM = posStrandPWM;
            addNegStrandPWM = negStrandPWM;
        }
        else {
            addPosStrandPWM = oldAMposStrandPWM;
            addNegStrandPWM = oldAMnegStrandPWM;
        }

        if (includePwm) {
            if (scoreType != null) {
                score += WeightMatrixTools.score(
                    addPosStrandPWM,
                    addNegStrandPWM,
                    posStrandSeq,
                    positionalWeights,
                    mandatoryColumns,
                    strand,
                    calc,
                    scoreType,
                    eToMu,
                    nonSpecKa,
                    revCompSimilarity);
            }
            else {
                score += WeightMatrixTools.score(
                    addPosStrandPWM,
                    addNegStrandPWM,
                    posStrandSeq,
                    positionalWeights,
                    mandatoryColumns,
                    strand,
                    calc,
                    eToMu,
                    nonSpecKa,
                    revCompSimilarity);
            }
            psamScore = score;
            // System.out.println("\nPSAM score = "+score);
        }

        // Iterate over each feature and add their contribution to the affinity score
        getFeatures(0, true);
        if (featuresTable != null) {

            // sort the table by the relative affinities if the table is in an unsorted state
            sortIfNeeded();

            double lastScore = score;

//             ArrayList<WeightMatrixFeature> nNDDFeatures = getNNDDFeatures(posStrandSeq);
//             for (WeightMatrixFeature feature : nNDDFeatures) {

//             ArrayList<WeightMatrixFeature> features = new ArrayList<WeightMatrixFeature>();
//             features.addAll(anFWM.getNNDDFeatures(aSeq, strand, startPos, startColumn, endColumn));

            for (Iterator iter = featuresTable.iterator(); iter.hasNext(); ) {
                List rowList = (List)iter.next();
                WeightMatrixFeature feature = (WeightMatrixFeature)rowList.get(0);

                // skip the feature if strand==BOTH and isRevComp
                //             if (strand==WeightMatrixTools.BindingStrand.BOTH && feature.isRevComp()) {
                //                 continue;
                //             }

                double relativeAffinity = feature.getRelativeAffinity();

                // if the relativeAffinity is below the threshold then stop adding feature corrections
                if (Math.abs(relativeAffinity) < threshold) {
                    System.out.println("\nRelative affinity "+relativeAffinity+" is less than the threshold "+threshold);
                    break;
                }

                // score the posStrandSeq with the featurePosStrandPWM
                double featureScore;
                if (scoreType != null) {
                    featureScore = WeightMatrixTools.score(
                        feature.getPosStrandPWM(),
                        feature.getNegStrandPWM(),
                        posStrandSeq,
                        positionalWeights,
                        feature.getMandatoryColumns(),
                        strand,
                        calc,
                        scoreType,
                        eToMu,
                        //nonSpecKa,
                        0,
                        feature.getRevCompSimilarity());
                }
                else {
                    featureScore = WeightMatrixTools.score(
                        feature.getPosStrandPWM(),
                        feature.getNegStrandPWM(),
                        posStrandSeq,
                        positionalWeights,
                        feature.getMandatoryColumns(),
                        strand,
                        calc,
                        eToMu,
                        //nonSpecKa,
                        0,
                        feature.getRevCompSimilarity());
                }

//                 if (!feature.isGood()) {
//                     relativeAffinity *= -1.0;
//                 }
                score += relativeAffinity * featureScore;

                // This feature applies in this sequence context iff the score has changed
                if (score != lastScore) {
                    // System.out.println("feature = "+feature.modsToString()+"\tupdated score = "+score);
                }
                lastScore = score;
            }
        }

        // If our motif isRevCompPalindrome and we looked at only the POS strand, then double our score
//         if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.POS)) {
//             score *= 2.0;
//         }

        if (score < 0) {
            if (nonSpecKa != 0) {
                score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
            }
            else {
                score = psamScore;
            }
        }

        return(score);
    }

    ///////////////////////////////////////
    // score with a startPos
    // use the mandatory columns of each feature
    ///////////////////////////////////////
    public double addModelScore(SymbolList posStrandSeq, FeatureKey blah, int startPos) {
        return(addModelScore(posStrandSeq, null, 0, startPos));
    }

    public double addModelScore(SymbolList posStrandSeq, ScoreType scoreType, int startPos) {
        return(addModelScore(posStrandSeq, scoreType, 0, startPos));
    }

    public double addModelScore(SymbolList posStrandSeq, ScoreType scoreType, double threshold, int startPos) {

        // If the PSAM is an reverse-complement palindrome, then we should look at only the positive
        // or negative strand (BUT NOT BOTH)
//         if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.BOTH)) {
//             //out.println("\n\nSince the motif is self reverse-complement, we should look only at the positive strand.");
//             strand = WeightMatrixTools.BindingStrand.POS;
//         }

        // score the posStrandSeq with the addPosStrandPWM
        double score = 0;
        double psamScore = 0;

        WeightMatrix addPosStrandPWM;
        WeightMatrix addNegStrandPWM;
        if (oldAMposStrandPWM == null) {
            addPosStrandPWM = posStrandPWM;
            addNegStrandPWM = negStrandPWM;
        }
        else {
            addPosStrandPWM = oldAMposStrandPWM;
            addNegStrandPWM = oldAMnegStrandPWM;
        }

        if (includePwm) {
            if (scoreType != null) {
                score += WeightMatrixTools.score(
                    addPosStrandPWM,
                    addNegStrandPWM,
                    posStrandSeq,
                    strand,
                    calc,
                    scoreType,
                    eToMu,
                    nonSpecKa,
                    revCompSimilarity,
                    startPos);
            }
            else {
                score += WeightMatrixTools.score(
                    addPosStrandPWM,
                    addNegStrandPWM,
                    posStrandSeq,
                    strand,
                    calc,
                    eToMu,
                    nonSpecKa,
                    revCompSimilarity,
                    startPos);
            }
            psamScore = score;
            // System.out.println("\nPSAM score = "+score);
        }

        // make an array of unity weights
        double[][] weights = new double[2][posStrandSeq.length() + this.columns - 1];
        Arrays.fill(weights[0], 1.0);
        Arrays.fill(weights[1], 1.0);

        // Iterate over each feature and add their contribution to the affinity score
        getFeatures(0, true);
        if (featuresTable != null) {

            // sort the table by the relative affinities if the table is in an unsorted state
            sortIfNeeded();

            double lastScore = score;

            for (Iterator iter = featuresTable.iterator(); iter.hasNext(); ) {
                List rowList = (List)iter.next();
                WeightMatrixFeature feature = (WeightMatrixFeature)rowList.get(0);

                // skip the feature if strand==BOTH and isRevComp
                //             if (strand==WeightMatrixTools.BindingStrand.BOTH && feature.isRevComp()) {
                //                 continue;
                //             }

                double relativeAffinity = feature.getRelativeAffinity();

                // if the relativeAffinity is below the threshold then stop adding feature corrections
                if (Math.abs(relativeAffinity) < threshold) {
                    System.out.println("\nRelative affinity "+relativeAffinity+" is less than the threshold "+threshold);
                    break;
                }

                // score the posStrandSeq with the featurePosStrandPWM
                double featureScore;
                if (scoreType != null) {
                    featureScore = WeightMatrixTools.score(
                        feature.getPosStrandPWM(),
                        feature.getNegStrandPWM(),
                        posStrandSeq,
                        weights,
                        feature.getMandatoryColumns(),
                        strand,
                        calc,
                        scoreType,
                        eToMu,
                        //nonSpecKa,
                        0,
                        feature.getRevCompSimilarity(),
                        startPos);
                }
                else {
                    featureScore = WeightMatrixTools.score(
                        feature.getPosStrandPWM(),
                        feature.getNegStrandPWM(),
                        posStrandSeq,
                        weights,
                        feature.getMandatoryColumns(),
                        strand,
                        calc,
                        eToMu,
                        //nonSpecKa,
                        0,
                        feature.getRevCompSimilarity(),
                        startPos);
                }

//                 if (!feature.isGood()) {
//                     relativeAffinity *= -1.0;
//                 }
                score += relativeAffinity * featureScore;

                // This feature applies in this sequence context iff the score has changed
                if (score != lastScore) {
                    // System.out.println("feature = "+feature.modsToString()+"\tupdated score = "+score);
                }
                lastScore = score;
            }
        }

        // If our motif isRevCompPalindrome and we looked at only the POS strand, then double our score
//         if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.POS)) {
//             score *= 2.0;
//         }

        if (score < 0) {
            if (nonSpecKa != 0) {
                score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
            }
            else {
                score = psamScore;
            }
        }

        return(score);
    }

    /////////////////////////////////////////////////////
    // score WITHOUT weights
    // does not use the mandatory columns of each feature
    /////////////////////////////////////////////////////
    public double addModelScore(SymbolList posStrandSeq, FeatureKey blah) {
        return(addModelScore(posStrandSeq, (ScoreType)null, 0));
    }

    public double addModelScore(SymbolList posStrandSeq, ScoreType scoreType) {
        return(addModelScore(posStrandSeq, scoreType, 0));
    }

    public double addModelScore(SymbolList posStrandSeq, ScoreType scoreType, double threshold) {

        // If the PSAM is an reverse-complement palindrome, then we should look at only the positive
        // or negative strand (BUT NOT BOTH)
//         if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.BOTH)) {
//             //out.println("\n\nSince the motif is self reverse-complement, we should look only at the positive strand.");
//             strand = WeightMatrixTools.BindingStrand.POS;
//         }

        // score the posStrandSeq with the posStrandPWM
        double score = 0;
        double psamScore = 0;

        WeightMatrix addPosStrandPWM;
        WeightMatrix addNegStrandPWM;
        if (oldAMposStrandPWM == null) {
            addPosStrandPWM = posStrandPWM;
            addNegStrandPWM = negStrandPWM;
        }
        else {
            addPosStrandPWM = oldAMposStrandPWM;
            addNegStrandPWM = oldAMnegStrandPWM;
        }

        if (includePwm) {
            if (scoreType != null) {
                score += WeightMatrixTools.score(
                    addPosStrandPWM,
                    addNegStrandPWM,
                    posStrandSeq,
                    strand,
                    calc,
                    scoreType,
                    eToMu,
                    nonSpecKa,
                    revCompSimilarity);
            }
            else {
                score += WeightMatrixTools.score(
                    addPosStrandPWM,
                    addNegStrandPWM,
                    posStrandSeq,
                    strand,
                    calc,
                    eToMu,
                    nonSpecKa,
                    revCompSimilarity);
            }
            psamScore = score;
            //System.out.println("\nPSAM score = "+score);
        }

        // Iterate over each feature and add their contribution to the affinity score
        getFeatures(0, true);
        if (featuresTable != null) {

            // sort the table by the relative affinities if the table is in an unsorted state
            sortIfNeeded();

            double lastScore = score;

            for (Iterator iter = featuresTable.iterator(); iter.hasNext(); ) {
                List rowList = (List)iter.next();
                WeightMatrixFeature feature = (WeightMatrixFeature)rowList.get(0);

                // skip the feature if strand==BOTH and isRevComp
                //             if (strand==WeightMatrixTools.BindingStrand.BOTH && feature.isRevComp()) {
                //                 continue;
                //             }

                double relativeAffinity = feature.getRelativeAffinity();

                // if the relativeAffinity is below the threshold then stop adding feature corrections
                if (Math.abs(relativeAffinity) < threshold) {
                    System.out.println("\nRelative affinity "+relativeAffinity+" is less than the threshold "+threshold);
                    break;
                }

                // score the posStrandSeq with the featurePosStrandPWM
                double featureScore;
                if (scoreType != null) {
                    featureScore = WeightMatrixTools.score(
                        feature.getPosStrandPWM(),
                        feature.getNegStrandPWM(),
                        posStrandSeq,
                        strand,
                        calc,
                        scoreType,
                        eToMu,
                        //nonSpecKa,
                        0,
                        feature.getRevCompSimilarity());
                }
                else {
                    featureScore = WeightMatrixTools.score(
                        feature.getPosStrandPWM(),
                        feature.getNegStrandPWM(),
                        posStrandSeq,
                        strand,
                        calc,
                        eToMu,
                        //nonSpecKa,
                        0,
                        feature.getRevCompSimilarity());
                }

//                 if (!feature.isGood()) {
//                     relativeAffinity *= -1.0;
//                 }
                score += relativeAffinity * featureScore;

                // This feature applies in this sequence context iff the score has changed
                if (score != lastScore) {
                    //System.out.println("feature = "+feature.modsToString()+"\tupdated score = "+score);
                }
                lastScore = score;
            }
        }

        // If our motif isRevCompPalindrome and we looked at only the POS strand, then double our score
//         if (isRevCompPalindrome && (strand == WeightMatrixTools.BindingStrand.POS)) {
//             score *= 2.0;
//         }

        if (score < 0) {
            if (nonSpecKa != 0) {
                score = nonSpecKa * (posStrandSeq.length() - this.columns + 1);
            }
            else {
                score = psamScore;
            }
        }

        return(score);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    public void setScoringKmerAffinities(double[] topPercentagePerKmer, boolean[] includeNegative) {
        if (kmerToAffinityMatrix == null) {
            return;
        }

        for (int i=0; i < kmerToAffinityMatrix.length; i++) {

            if (topPercentagePerKmer[i] >= 1) {
                continue;
            }

            int maxWord = KmerMatrix.multBy4ToPower(1, kmerLengths[i]);

            Table kmersToAffinitiesTable = new Table();
            for (int index=0; index < maxWord; index++) {

                java.util.List keyValueList = Arrays.asList(
                    new Integer(index),
                    new Double(kmerToAffinityMatrix[i][index]),
                    new Double(Math.abs(kmerToAffinityMatrix[i][index])));

                kmersToAffinitiesTable.add(keyValueList);
            }

            if (includeNegative[i]) {
                kmersToAffinitiesTable.sort(2); // sort by ascending abs(affinity)
            }
            else {
                kmersToAffinitiesTable.sort(1); // sort by ascending affinity
            }

            kmersToAffinitiesTable.reverse(); // sort by descending
            int thresh = (int) Math.ceil(topPercentagePerKmer[i] * maxWord);

            // set affinities to 0 for those kmers above thresh
            for (int index=thresh; index <kmersToAffinitiesTable.rows(); index++) {
                kmersToAffinitiesTable.setElement(index, 1, new Double(0));
            }

            kmersToAffinitiesTable.sort(0); // sort by ascending words
            java.util.List newKmerAffinities = kmersToAffinitiesTable.getColumn(1);

            // add the intercept if exists
            if (kmerToAffinityMatrix[i].length > maxWord) {
                newKmerAffinities.add(kmerToAffinityMatrix[i][maxWord]);
            }

            kmerToAffinityMatrix[i] = ArrayTools.toDoubleArray((ArrayList)newKmerAffinities);
        }
    }

    public void cleanKmerAffinities() {
        if (kmerToAffinityMatrix != null) {
            for (int i=0; i < kmerToAffinityMatrix.length; i++) {
                for (int j=0; j < kmerToAffinityMatrix[i].length; j++) {
                    if (Double.isNaN(kmerToAffinityMatrix[i][j])) {
                        kmerToAffinityMatrix[i][j] = 0;
                    }
                }
            }
        }
    }


    double kmerScore(SymbolList aSymList, boolean includePosBias) {
        double score = 0;
        for (int i=0; i < kmerLengths.length; i++) {
        //for (int i=0; i < 1; i++) {
            score += kmerScore(aSymList, i, includePosBias);
        }

        return(score);
    }

    double kmerScore(SymbolList aSymList, int kmerIndex, boolean includePosBias) {
        double[][] weights;

        if (includePosBias) {
            int trainingProbeSeqLengths = this.positionalWeights[0].length - this.columns + 1;

            //weights = this.kmerPositionalWeights[kmerIndex];

            if (aSymList.length() == trainingProbeSeqLengths) {
                weights = this.kmerPositionalWeights[kmerIndex];
            }
            else {
                weights = null;
            }
        }
        else {
            weights = null;
        }

        return(kmerScore(aSymList, kmerIndex, weights));
    }

    double kmerScore(SymbolList aSymList, int kmerIndex, double[][] weights) {
        //throws IllegalAlphabetException, IllegalSymbolException, BioException{

        double score = 0;

        try {

            char kmerLength = kmerLengths[kmerIndex];

//             // HACK!!
//             if ((kmerLength == 8) && (weights[0][0] != 1.0)) {
//                 // set negative weights to positive weights !!
//                 //weights[1] = weights[0];

//                 // use revComps!!
//                 this.includeRevComps[kmerIndex] = true;
//             }

            // if (weights == null) {
            //     weights = this.kmerPositionalWeights[kmerIndex];
            // }

            // add the intecept if it exists
            int maxWord = (int)Math.pow(4, kmerLength);
            if (kmerToAffinityMatrix[kmerIndex].length > maxWord) {
                score += kmerToAffinityMatrix[kmerIndex][maxWord];
            }

            // sum over each kmer window
            for (int startPos = 1; startPos <= aSymList.length() - kmerLength + 1; startPos++) {
                double affinity = 0;

                // Positive Strand
                // subList(start_incl, end_incl), 1-based indexing
                // SymbolList aKmerSymList = aSymList.subList(startPos, startPos+kmerLength-1);
                // char aWordIndex = (char)KmerMatrix.getWordIndex(aKmerSymList, null);
                char aWordIndex = (char)KmerMatrix.getWordIndex(aSymList, startPos, startPos+kmerLength-1);

                if (weights != null) {
                    affinity += weights[0][startPos-1] * kmerToAffinityMatrix[kmerIndex][aWordIndex];
                }
                else {
                    affinity += kmerToAffinityMatrix[kmerIndex][aWordIndex];
                }

                if (this.includeRevComps[kmerIndex]) {
                    // Negative Strand
                    char aRevCompIndex = revCompMatrix[kmerIndex][aWordIndex];


                    if (weights != null) {
                        affinity += weights[1][startPos-1] * kmerToAffinityMatrix[kmerIndex][aRevCompIndex];
                    }
                    else {
                        affinity += kmerToAffinityMatrix[kmerIndex][aRevCompIndex];
                    }

                    // self-reverse-complement words have been counted twice
                    if (aRevCompIndex == aWordIndex) {
                        affinity = affinity / 2.0;
                    }

                }

                score += affinity;

//                 if (aRevCompIndex != aWordIndex) {
//                     //double affinity = kmerToAffinityMap.get(aRevCompIndex);
//                     affinity = kmerToAffinityMatrix[kmerIndex][aRevCompIndex];
//                     score += affinity * weights[1][startPos-1];
//                 }

            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

        return(score);
    }

    //////////////////////////////////////////////////////////////////////////////
    // COMBINED SCORE = FSAM score plus Kmer score to residuals
    //////////////////////////////////////////////////////////////////////////////
    public double modelScore(SymbolList posStrandSeq, boolean includeScaler, boolean includeKmers, boolean includeFeatures, boolean includePosBias) {
        double score = 0;

        // score fsam or psam first
        if ((includeFeatures) && (featureKeyToMultFeatureMap != null)) {
            //System.out.println("\nUsing the Additive Model to score the sequence.");
            //System.out.println("\nUsing the Multiplicative Model to score the sequence.");
            if (includePosBias) {
                //score += addModelScore(posStrandSeq, this.positionalWeights, (FeatureKey)null);
                score += multModelScore(posStrandSeq, this.positionalWeights, (FeatureKey)null);
            }
            else {
                //score += addModelScore(posStrandSeq, (FeatureKey)null);
                score += multModelScore(posStrandSeq);
            }
        }
        else {
            //System.out.println("\nUsing the PSAM Model to score the sequence.");
            if (includePosBias) {
                score += psamScore(posStrandSeq, this.positionalWeights);
            }
            else {
                score += psamScore(posStrandSeq);
            }
        }

        // score kmers last
        if ((includeKmers) && (kmerToAffinityMatrix != null)) {

            if (includeScaler) {
                //if ((includeScaler) && (this.intensitiesScaler > 1)) {
                score *= this.intensitiesScaler;
                //score *= Math.abs(scaler);

                score += this.intercept;
            }

            //System.out.println("\nUsing the Kmer Model to score the sequence.");
            score += kmerScore(posStrandSeq, includePosBias);
        }

        // score can be negative if unconstrained regression was performed to learn
        // the kmer affinities!!
        if (score <= 0) {
            score = 0.0001;
        }

        //System.exit(-1);
        return(score);
    }

    public double[] modelScoreWindows(SymbolList posStrandSeq, boolean includeFeatures, double minThreshold) {
        double[] scores = null;

        // score fsam or psam first
        if ((includeFeatures) && (featureKeyToMultFeatureMap != null)) {
            //System.out.println("\nUsing the Additive Model to score the sequence.");
            //System.out.println("\nUsing the Multiplicative Model to score the sequence.");
            //score += addModelScore(posStrandSeq, (FeatureKey)null);

            scores= multModelScoreWindows(
                this,
                posStrandSeq,
                strand,
                calc,
                eToMu,
                nonSpecKa,
                //this.getRevCompSimilarity(),
                revCompSimilarity,
                null //withholdFeatureKey
                // minThreshold
                                          );

        }
        else {
            //System.out.println("\nUsing the PSAM Model to score the sequence.");
            scores = WeightMatrixTools.scoreWindows(
                posStrandPWM,
                negStrandPWM,
                posStrandSeq,
                strand,
                calc,
                eToMu,
                nonSpecKa,
                revCompSimilarity
                // minThreshold
                                                    );
        }

        // score can be negative if unconstrained regression was performed to learn
        // the kmer affinities!!
        scores = MathTools.lowerBound(scores, 0.0001);

        //System.exit(-1);
        return(scores);
    }

    public int hashCode() {
        int hc = 17; // always start with a non-zero prime

        // posStrandPWM
        for (int c = 0; c < posStrandPWM.columns(); ++c) {
            hc = (23 * hc) + posStrandPWM.getColumn(c).hashCode();
        }

        // multFeatures
        if (featureKeyToMultFeatureMap != null) {
            for (FeatureKey aFeatureKey: featureKeyToMultFeatureMap.keySet()) {
                // System.out.println("FeatureKey "+keyNum+" ="+StringTools.toString(aFeatureKey.getModsArray(), parser, "", "-")+"  ; hashcode="+aFeatureKey.hashCode());
                WeightMatrixFeature feature = (WeightMatrixFeature)featureKeyToMultFeatureMap.get(aFeatureKey);
                hc = (23 * hc) + feature.hashCode();
            }
        }

        return hc;
    }

    public boolean equals(Object o) {
        if (o instanceof FeaturedWeightMatrix) {

            FeaturedWeightMatrix anFwm = (FeaturedWeightMatrix) o;
            if (anFwm.columns() != this.columns()) {
                return false;
            }
            if (anFwm.getAlphabet() != this.getAlphabet()) {
                return false;
            }

            // posStrandPWM
            for (int c = 0; c < posStrandPWM.columns(); ++c) {
                if (! posStrandPWM.getColumn(c).equals(anFwm.getPosStrandPWM().getColumn(c))) {
                    return false;
                }
            }

            // multFeatures
            if (featureKeyToMultFeatureMap != null) {
                for (FeatureKey aFeatureKey: featureKeyToMultFeatureMap.keySet()) {
                    // System.out.println("FeatureKey "+keyNum+" ="+StringTools.toString(aFeatureKey.getModsArray(), parser, "", "-")+"  ; hashcode="+aFeatureKey.hashCode());
                    WeightMatrixFeature thisFeature = (WeightMatrixFeature)featureKeyToMultFeatureMap.get(aFeatureKey);
                    WeightMatrixFeature aFeature = (WeightMatrixFeature)anFwm.getFeature(aFeatureKey);

                    if ((aFeature == null) || (!thisFeature.equals(aFeature))) {
                        return(false);
                    }
                }
            }

            return true;
        }
        return false;
    }

}
