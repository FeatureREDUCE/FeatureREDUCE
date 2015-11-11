/*
 * WeightMatrixBuilder.java - Todd Riley
 *
 */


import org.biojava.bio.symbol.*;
import org.biojava.bio.dist.*;
import org.biojava.utils.*;
import org.biojava.bio.dp.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.db.*;
import org.biojava.bio.seq.io.*;
import java.io.*;
import java.util.*;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.math.BigInteger;

import org.biojava.bio.mydp.*;


import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.iterator.TLongIntIterator;
import gnu.trove.map.hash.TLongIntHashMap;
import gnu.trove.set.hash.TLongHashSet;

public class WeightMatrixBuilder {

	private static String tabs = "\\t+";
	private static Pattern tabsPattern = Pattern.compile(tabs);
	FiniteAlphabet alphabet = DNATools.getDNA();
	ReversibleTranslationTable complementTable = DNATools.complementTable();

	public void parseCommandLine(String[] args) throws IOException {

		try {
			String title = "WeightMatrixBuilder";
			int numMandatoryArgs = 1;

			boolean batch = false;
			int countMin = 0;
			String displayString = "Yes";
			boolean display = true;
			String scaleByInfoContentString = "Yes";
			boolean scaleByInfoContent = true;
			String outlineString = "Yes";
			boolean outline = true;
			String centeredNumberingString = "Yes";
			boolean centeredNumbering = true;
			String revCompString = "No";
			boolean revComp = false;
			double affinitySphereThresh = 0;
			String affinitySphereFileName = null;
			float offset = 0;
			String[] rearrangeStrings = null;
			WeightMatrixLogo.PositionLoc posLoc = WeightMatrixLogo.PositionLoc.BOTTOM;
			double nullWeight = 0;
			String pwm1Format = null;
			String pwm2Format = null;
			String pwm1PathName = null;
			String pwm2PathName = null;
			String saveFormat = null;
			String saveFile = null;
			String saveLogoFormat = null;
			String saveLogoFile = null;
			int start = 0;
			int motifLength = -1;
			String regexString = null;
			String convertToString = null;

			String sequence = null;
			double pseudoCountWeight = 0;

			String psamListFile = null;
			double distanceMeasure = 0;

			boolean enrichment = false;
			// int[] columns = {6, 7, 8, 9, 10, 11};
			int[] columns = { 6, 8 };

            double weight1 = -1;
            double weight2 = -1;

			WeightMatrixLogo logoPanel = null;

            System.out.println("\n"	+ title	+ " - A weight matrix building tool that reads in and parses sequences to generate, display, and save PWMs or PSAMs.\n");
            System.out.println("\t\t    - written by Todd R. Riley\n");

			if ((!(args.length == 0 && numMandatoryArgs == 0))
                && (args.length < numMandatoryArgs || args[0].equalsIgnoreCase("-help") || args[0].equalsIgnoreCase("-?")))
            {
				System.out.println("\n Usage: " + title	+ " <Counts-File1> <Counts-File2> -Option\n");
				System.out.println(" Options:");
				System.out.println("\n   -start <Int> (default " + start + ")");
				System.out.println();
				System.out.println("\n\t The start position in the sequences used to find the motif (0-based indexing).");
				System.out.println();
				System.out.println("\n   -length <Int> (default " + "max" + ")");
				System.out.println();
				System.out.println("\n\t The length of the motif to search for.");
				System.out.println();
				System.out.println("\n   -countMin <Int> (default " + countMin	+ ")");
				System.out.println();
				System.out.println("\t Sets a minimum occurrence count in order to include sequence in generation of PWM.");
				System.out.println();
				System.out.println("\n   -match <regex>");
				System.out.println();
				System.out.println("\t A regex expression for a subsubsequence within the subsequence (specified by -start and -length) that must");
				System.out.println("\t match in order to include the subsubsequence in the generation of the PWM.");
				System.out.println();
				System.out.println("\n   -revComp <Yes/No> (default " + revCompString + ")");
				System.out.println();
				System.out.println("\t Turns on/off including the reverse complement of the sequences.");
				System.out.println();
				System.out.println("\n   -displayMotifs <Yes/No> (default " + displayString + ")");
				System.out.println();
				System.out.println("\t Turns on/off the display of PWM logo.");
				System.out.println();
				System.out.println("\n   -scaleByIC <Yes/No> (default "	+ scaleByInfoContentString + ")");
				System.out.println();
				System.out.println("\t Turns on/off scaling each column in the PWM logo by its information content.");
				System.out.println();
				System.out.println("\n   -outline <Yes/No> (default " + outlineString + ")");
				System.out.println();
				System.out.println("\t Turns on/off black outline of the PWM logo.");
				System.out.println();
				System.out.println("\n   -centerNum <Yes/No> (default "	+ centeredNumberingString + ")");
				System.out.println();
				System.out.println("\t Turns on/off centered numbering of the columns in the PWM logo.");
				System.out.println();
				System.out.println("\n   -offset <Int> (default " + offset	+ ")");
				System.out.println();
				System.out.println("\t Adds an offset to the numbering of the columns in the PWM logo.");
				System.out.println();
				System.out.println("   -rearrange \t <Columns/~Columns>");
				System.out.println();
				System.out.println("\t Rearrange, copy, drop, and/or complement columns (1-based indexing).");
				System.out.println();
				System.out.println("\t examples:  reverse  5 4 3 2 1");
				System.out.println("\t\t   reverse complement  ~5 ~4 ~3 ~2 ~1");
				System.out.println();
				System.out.println("   -weightedAvg \t <weight1> <weight2>");
				System.out.println();
				System.out.println("\t Calculate the weighted average between two models.");
				System.out.println();
				System.out.println("\t example:  -load xml model1.xml -load xml model2.xml -weightedAvg 1 3");
				System.out.println();
				System.out.println("\n   -sequence <motif> <pseudoCountWeight>");
				System.out.println();
				System.out.println("\t Creates a weight matrix representation of a motif smoothed by pseudocounts.");
				System.out.println();
				System.out.println("\n   -load <ser/table/xml/mxr> <File-Name>");
				System.out.println();
				System.out.println("\t Loads a PWM saved in XML or SER(ialized) format.");
				System.out.println();
				System.out.println("\n   -convert <probability/relAffinity>");
				System.out.println();
				System.out.println("\t converts a PWM to either a probability or relative affinity PWM.");
				System.out.println();
				System.out.println("   -save <ser/xml/table> <File-Name>");
				System.out.println();
				System.out.println("\t Saves the PWM saved in XML or SER(ialized) format.");
				System.out.println();
				System.out.println("   -saveLogo <png/jpg> <File-Name>");
				System.out.println();
				System.out.println("\t Saves the PWM Logo.");
				System.out.println();
				System.out.println("   -affinitySphere <thresh> <File-Name>");
				System.out.println();
				System.out.println("\t Saves the relative affinity sphere all sequences with relative affinity >= threshold.");
				System.out.println();
				System.out.println("   -enrichment");
				System.out.println();
				System.out.println("\t Use -batch to calculate the enrichment from PSAM1 to PSAM2.");
				System.out.println();
				System.out.println("   -batch");
				System.out.println();
				System.out.println("\t Use -batch to force termination after execution and plot generation.");
				System.out.println();
				System.out.println("\n   -help, -? \t Displays this help message");
				System.out.println();
				System.exit(-1);
			}

			for (int a = 0; a < args.length; a++) {
				if (args[a].equals("-load")) {
					if (pwm1PathName == null) {
						pwm1Format = new String(args[a + 1]);
						pwm1PathName = new String(args[a + 2]);
					}
                    else {
						pwm2Format = new String(args[a + 1]);
						pwm2PathName = new String(args[a + 2]);
					}
				}
                else if (args[a].equals("-start")) {
					start = Integer.parseInt(args[a + 1]);
				}
                else if (args[a].equals("-length")) {
					motifLength = Integer.parseInt(args[a + 1]);
				}
                else if (args[a].equals("-match")) {
					regexString = new String(args[a + 1]);
				}
                else if (args[a].equals("-revComp")) {
					revCompString = new String(args[a + 1]);
					revComp = StringTools.parseBoolean(revCompString);
				}
                else if (args[a].equals("-countMin")) {
					countMin = Integer.parseInt(args[a + 1]);
				}
                else if (args[a].equals("-displayMotifs")) {
					displayString = new String(args[a + 1]);
					display = StringTools.parseBoolean(displayString);
				}
                else if (args[a].equals("-scaleByIC")) {
					scaleByInfoContentString = new String(args[a + 1]);
					scaleByInfoContent = StringTools
                        .parseBoolean(scaleByInfoContentString);
				}
                else if (args[a].equals("-outline")) {
					outlineString = new String(args[a + 1]);
					outline = StringTools.parseBoolean(outlineString);
				}
                else if (args[a].equals("-centerNum")) {
					centeredNumberingString = new String(args[a + 1]);
					centeredNumbering = StringTools.parseBoolean(centeredNumberingString);
				}
                else if (args[a].equals("-offset")) {
					offset = Float.parseFloat(args[a + 1]);
				}
                else if (args[a].equalsIgnoreCase("-rearrange")) {
					rearrangeStrings = getParams(args, a + 1);
				}
                else if (args[a].equalsIgnoreCase("-weightedAvg")) {
					weight1 = Double.parseDouble(args[a + 1]);
					weight2 = Double.parseDouble(args[a + 2]);
				}
                else if (args[a].equalsIgnoreCase("-save")) {
					saveFormat = new String(args[a + 1]);
					saveFile = new String(args[a + 2]);
				}
                else if (args[a].equalsIgnoreCase("-saveLogo")) {
					saveLogoFormat = new String(args[a + 1]);
					saveLogoFile = new String(args[a + 2]);
				}
                else if (args[a].equalsIgnoreCase("-convert")) {
					convertToString = new String(args[a + 1]);
				}
                else if (args[a].equals("-batch")) {
					batch = true;
				}
                else if (args[a].equalsIgnoreCase("-affinitySphere")) {
					affinitySphereThresh = Double.parseDouble(args[a + 1]);
					affinitySphereFileName = new String(args[a + 2]);
				}
                else if (args[a].equalsIgnoreCase("-similarities")) {
					psamListFile = new String(args[a + 1]);
					distanceMeasure = Double.parseDouble(args[a + 2]);
				}
                else if (args[a].equalsIgnoreCase("-sequence")) {
					sequence = new String(args[a + 1]);
					pseudoCountWeight = Double.parseDouble(args[a + 2]);
				}
				// ////////////////////////////////////////////////////////////////////////////
				// ////////////////////////////////////////////////////////////////////////////
				else if (args[a].startsWith("-")) {
					System.out.println("\n\nError: Unknown argument " + args[a]
                        + "\n");
					System.exit(-1);
				}
			}

			WeightMatrix pwm1 = null;
			WeightMatrix pwm2 = null;
			WeightMatrix pwm3 = null;

			if (sequence != null) {
				Distribution backgroundDist = DistributionTools
                    .createUniformDistribution(this.alphabet, false);
				SymbolList symList = DNATools.createDNA(sequence);

				pwm1 = WeightMatrixTools.getWeightMatrix(symList,
                    pseudoCountWeight, backgroundDist);
			}
            else if (pwm1PathName != null) {
				pwm1 = loadPsam(pwm1Format, pwm1PathName);
			}
            else if (psamListFile == null) {
				pwm1 = build(args[0], start, motifLength, regexString, countMin, nullWeight, revComp);
			}

			if (convertToString != null) {
				if (convertToString.equalsIgnoreCase("probability")
                    && !WeightMatrixTools.isProbabilityMatrix(pwm1)) {
					pwm1 = WeightMatrixTools.getProbabilities(pwm1);
				}
                else if (convertToString.equalsIgnoreCase("relAffinity")
                    && WeightMatrixTools.isProbabilityMatrix(pwm1)) {
					pwm1 = WeightMatrixTools.getRelativeAffinities(pwm1);
				}
			}

			// display Weight Matrix
            if (pwm1 != null) {
                if (display) {
                    logoPanel = displayPWM("PWM1", pwm1, scaleByInfoContent,
						outline, centeredNumbering, offset, posLoc);
                }
                System.out.print("\n" + StringTools.toString(pwm1));
            }

			if ((pwm2PathName != null)
                || (!args[0].startsWith("-") && !args[1].startsWith("-"))) {

				if (pwm2PathName != null) {
					pwm2 = loadPsam(pwm2Format, pwm2PathName);
				}
                else {
					pwm2 = build(args[1], start, motifLength, regexString, countMin, nullWeight, revComp);
				}

				if (convertToString != null) {
					if (convertToString.equalsIgnoreCase("probability")
                        && !WeightMatrixTools.isProbabilityMatrix(pwm2)) {
						pwm2 = WeightMatrixTools.getProbabilities(pwm2);
					}
                    else if (convertToString.equalsIgnoreCase("relAffinity")
                        && WeightMatrixTools.isProbabilityMatrix(pwm2)) {
						pwm2 = WeightMatrixTools.getRelativeAffinities(pwm2);
					}
				}

				if (display) {
					logoPanel = displayPWM("PWM2", pwm2, scaleByInfoContent,
                        outline, centeredNumbering, offset, posLoc);
				}
				System.out.println("\n" + StringTools.toString(pwm2));
			}


            if ((weight1 != -1) && (pwm1 != null) && (pwm2 != null)) {
                pwm3 = WeightMatrixTools.weightedAverage(weight1, weight2, pwm1, pwm2);
            }
			else if (enrichment) {
				pwm3 = WeightMatrixTools.getRelativeEnrichment(pwm1, pwm2);
				pwm3 = WeightMatrixTools.getProbabilities(pwm3);

				if (display) {
					logoPanel = displayPWM("Relative Enrichment from PWM1 to PWM2", pwm3, scaleByInfoContent, outline, centeredNumbering, offset, posLoc);
				}
			}
            else {
				pwm3 = pwm1;
			}


            if (rearrangeStrings != null) {
				if (pwm2 == null) {
					pwm3 = WeightMatrixTools.rearrange(pwm3, rearrangeStrings, complementTable);
				}
                else {
					pwm3 = WeightMatrixTools.rearrange(pwm1, pwm2, rearrangeStrings, complementTable);
				}

				if (display) {
					logoPanel = displayPWM("Rearranged Logo", pwm3, scaleByInfoContent, outline, centeredNumbering, offset, posLoc);
				}
			}

			// Save the PWM
			if (saveFile != null) {
				save(pwm3, saveFormat, saveFile);
			}

			// Save the PWM Logo
			if (saveLogoFile != null) {
				GraphicsTools.writeImage(logoPanel, saveLogoFile, false);
			}

			if (affinitySphereFileName != null) {
				System.out.print("\n"
                    + "Creating the Relative Affinity Sphere that contains all sequences with relative affinity >= "
                    + affinitySphereThresh + "...");

				LinkedHashMap<SymbolList, Double> relAffinitySphere = WeightMatrixTools.getRelAffinitySphere(
                    pwm3,
                    true,
                    affinitySphereThresh,
                    DNATools.complementTable(),
                    WeightMatrixTools.BindingStrand.BOTH,
                    WeightMatrixTools.BothStrandsCalc.NORMED_SUM);

				Table affinitySphereTable = new Table(relAffinitySphere);
				affinitySphereTable.sort(1);
				affinitySphereTable.reverse();
				System.out.println(" Done.");

				// FileTools.write(affinitySphereTable.toString("\t"),
				// affinitySphereFileName, false);
				System.out.print("\n"
                    + "Writing Relative Affinity Sphere to file "
                    + affinitySphereFileName + "...");
				StringBuffer stringBuffer = new StringBuffer();
				stringBuffer.append("Kmer" + "\t" + "relAffinity");

				for (Iterator rowIter = affinitySphereTable.iterator(); rowIter.hasNext();) {
					List currentRow = (List) rowIter.next();
					String kmer = ((SymbolList) currentRow.get(0)).seqString().toUpperCase();
					Double relAffinity = (Double) currentRow.get(1);
					stringBuffer.append("\n" + kmer + "\t" + relAffinity);
				}
				FileTools.write(stringBuffer.toString(), affinitySphereFileName, false);
				System.out.println(" Done.");
			}

			if (psamListFile != null) {
				String[] psamFileNames = FileTools.readStrings(psamListFile);
				WeightMatrix[] psams = new WeightMatrix[psamFileNames.length];
				for (int i = 0; i < psamFileNames.length; i++) {
					psams[i] = WeightMatrixTools.getProbabilities(WeightMatrixTools.readFromXML(psamFileNames[i]));

					String baseName = FileTools.stripPath(psamFileNames[i]);
					String[] pars = baseName.split("\\.");
					String protein = pars[0];

					psams[i].setName(protein);

				}

				System.out.println("\n" + WeightMatrixTools.similarities(psams, distanceMeasure, columns));

			}

			// if batch mode then terminate
			if (batch) {
				System.exit(0);
			}

		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}

	private String[] getParams(String[] argsList, int index) {
		int counter = 0;
		while ((index + counter < argsList.length)
            && (!argsList[index + counter].startsWith("-"))) {
			counter++;
		}
		String[] params = new String[counter];
		for (int i = 0; i < counter; i++) {
			params[i] = argsList[index + i];
		}
		// System.out.println("params_end=" + params[params.length-1]);
		return (params);
	}

	// builds a probability PWM
	public WeightMatrix build(String inputFile, int start, int motifLength,
        String regexString, int countMin, double nullWeight, boolean revComp) {

		WeightMatrix pwm = null;
		try {
			SymbolTokenization symToken = alphabet.getTokenization("token");

			String countsTableBaseName = FileTools.stripPath(inputFile);
			BufferedReader bufferedReader = new BufferedReader(new FileReader(
					inputFile));

			// 0-based indexing
			// int varRegionLength = 16;
			// int varRegionStart = 29;

			Distribution[] dists = null;
			DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();
			dtc.setNullModelWeight(nullWeight);

			Pattern regexPattern = null;
			if (regexString != null) {
				regexPattern = Pattern.compile(regexString);
			}

			// Parse the table and populate the kmerToCountMap
			System.err.print("Parsing " + countsTableBaseName
                + " and populating nucleotide counts per position....");
			int numLines = FileTools.getLineCount(inputFile);

			int lineCounter = 0;
			String lineString;
			while ((lineString = bufferedReader.readLine()) != null) {
				// String lineEntries[] = lineString.split("\t");

				lineCounter++;
				if ((lineCounter % 1000) == 0) {
					System.err.print("\n\tParsed " + lineCounter + " out of " + numLines + " sequences so far...");
				}

				String lineEntries[] = tabsPattern.split(lineString);
				if (lineEntries.length > 1) {

					String seqUpperCase = lineEntries[0].toUpperCase();
					int seqReadCount = Integer.parseInt(lineEntries[1]);

					// System.err.println("tseqUpperCase="+seqUpperCase+"\t"+"seqReadCount="+seqReadCount);

					if (seqReadCount < countMin) {
						continue;
					}

					String origRegionString = null;
					if (motifLength == -1) {
						origRegionString = seqUpperCase;
					}
                    else {
						origRegionString = seqUpperCase.substring(start, start + motifLength);
					}

					String[] regionStrings = null;
					if (revComp) {
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
						if (regexPattern != null) { // use each substring that
													// matches the regex
							Matcher matcher = regexPattern.matcher(regionString);
							while (matcher.find()) {
								SymbolList symList = DNATools.createDNA(matcher.group());
								dists = addCounts(symList, seqReadCount, dists,	dtc);
							}
						}
                        else { // use the whole regionString
							SymbolList symList = DNATools.createDNA(regionString);
							dists = addCounts(symList, seqReadCount, dists, dtc);
						}
						// System.err.println(regionString);
					}

				}
			}
			System.err.println(" Done.");

			System.err.println("Creating the nucleotide emission distributions for each position...");
			dtc.train();
			System.err.println(" Done.");

			// make a Weight Matrix
			pwm = new SimpleWeightMatrix(dists);
		}

		catch (Exception ex) {
			ex.printStackTrace();
		}

		return (pwm);
	}

	private Distribution[] addCounts(SymbolList symList, int seqReadCount,
        Distribution[] dists, DistributionTrainerContext dtc) {
		try {
			// initialize if needed
			if (dists == null) {
				dists = new Distribution[symList.length()];
				for (int j = 0; j < dists.length; j++) {
					dists[j] = DistributionFactory.DEFAULT.createDistribution(alphabet);
					dtc.registerDistribution(dists[j]);
				}
			}

			for (int i = 0; i < symList.length(); i++) {
				// 1-based indexing
				Symbol s = symList.symbolAt(i + 1);
				dtc.addCount(dists[i], s, seqReadCount);// increment the count
														// for the symbol
			}
		}

		catch (Exception ex) {
			ex.printStackTrace();
		}
		return (dists);
	}

	public static WeightMatrixLogo displayPWM(String title, WeightMatrix pwm,
        boolean scaleByInfoContent, boolean outline,
        boolean centeredNumbering, float offset,
        WeightMatrixLogo.PositionLoc posLoc) {

		GridBagLayoutFrame logoFrame = new GridBagLayoutFrame(title, true);

		WeightMatrixLogo logoPanel = new WeightMatrixLogo(pwm, scaleByInfoContent, outline, centeredNumbering, offset, posLoc, null);

		logoFrame.add(logoPanel, 10, 10, false, true);

		return (logoPanel);
	}

	public void save(WeightMatrix pwm, String format, String filePathName) {
		if (format.equalsIgnoreCase("xml")) {
			WeightMatrixTools.writeToXML(pwm, filePathName);
		}
        else if (format.equalsIgnoreCase("table")) {
			WeightMatrixTools.writeToTable(pwm, filePathName);
		}
        else {
			FileTools.writeSerializedFile((SimpleWeightMatrix)pwm, filePathName);
		}
	}

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
        else {
			pwm = (WeightMatrix) FileTools.readSerializedFile(filePathName);
		}
		pwm.setName(FileTools.stripPath(filePathName));
		return (pwm);
	}

	public static void main(String[] args) throws IOException {
		WeightMatrixBuilder builder = new WeightMatrixBuilder();
		builder.parseCommandLine(args);
	}

}
