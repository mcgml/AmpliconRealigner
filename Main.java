package nhs.genetics.cardiff;

import htsjdk.samtools.*;
import nhs.genetics.cardiff.framework.*;
import org.apache.commons.cli.*;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Program for realigning soft clipped amplicon reads (post merging and whole genome alignment)
 *
 * @author  Matt Lyon
 * @version 1.0
 * @since   2016-10-11
 */
public class Main {

    private static final Logger log = Logger.getLogger(Main.class.getName());
    private static final String program = "AmpliconRealigner";
    private static final String version = "1.1.0";

    public static void main(String[] args) {

        CommandLineParser commandLineParser = new DefaultParser();
        CommandLine commandLine = null;
        HelpFormatter formatter = new HelpFormatter();
        Options options = new Options();

        options.addOption("I", "Input", true, "Path to input BAM file");
        options.addOption("T", "Targets", true, "Path to BED file");
        options.addOption("O", "Output", true, "Path to output BAM file");
        options.addOption("R", "Reference", true, "Path to FASTA file");
        options.addOption("S", "MinScore", true, "Minimum alignment score [50]");
        options.addOption("P", "PrimerSimilarity", true, "Primer similarity [0.8]");
        options.addOption("Go", "GapOpen", true, "Read alignment gap open penalty [-14]");
        options.addOption("Ge", "GapExtend", true, "Read alignment gap extend penalty [-4]");
        options.addOption("V", "Verbose", false, "Verbose logging [false]");

        try {
            commandLine = commandLineParser.parse(options, args);

            if (!commandLine.hasOption("I") || !commandLine.hasOption("T") || ! commandLine.hasOption("O") || ! commandLine.hasOption("R")){
                throw new NullPointerException("Incorrect arguments");
            }

        } catch (ParseException | NullPointerException e){
            formatter.printHelp(program + " " + version, options);
            log.log(Level.SEVERE, e.getMessage());
            System.exit(-1);
        }

        File inputSamOrBamFile = new File(commandLine.getOptionValue("I"));
        File outputSamOrBamFile = new File(commandLine.getOptionValue("O"));
        File bedFile = new File(commandLine.getOptionValue("T"));
        File referenceFasta = new File(commandLine.getOptionValue("R"));
        File referenceFastaFai = new File(commandLine.getOptionValue("R") + ".fai");
        ArrayList<GenomicLocation> genomicLocations = null;
        int minScore = commandLine.hasOption("S") ? Integer.parseInt(commandLine.getOptionValue("S")) : 50;
        int gapOpenPenalty = commandLine.hasOption("Go") ? Integer.parseInt(commandLine.getOptionValue("Go")) : -14;
        int gapExtendPenalty = commandLine.hasOption("Ge") ? Integer.parseInt(commandLine.getOptionValue("Ge")) : -4;
        double primerSimilarity = commandLine.hasOption("P") ? Double.parseDouble(commandLine.getOptionValue("P")) : 0.8;

        log.log(Level.INFO, "Running with settings: minScore=" +  minScore + " gapOpenPenalty=" + gapOpenPenalty + " gapExtendPenalty=" + gapExtendPenalty + " primerSimilarity=" + primerSimilarity);

        log.log(Level.INFO, "Reading BED file: " + bedFile + " ...");
        try {
            genomicLocations = BEDFile.getBedFeatures(bedFile);
        } catch (IOException e){
            log.log(Level.SEVERE, "Could not read BED file: " + e.getMessage());
            System.exit(-1);
        }

        log.log(Level.INFO, "Reading BAM file: " + inputSamOrBamFile + " ...");
        try (SamReader samReader = SamReaderFactory.makeDefault().open(inputSamOrBamFile)){

            SAMFileHeader samFileHeader = samReader.getFileHeader();
            samFileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            SAMProgramRecord samProgramRecord = samFileHeader.createProgramRecord();
            samProgramRecord.setCommandLine(String.join(" ", args));
            samProgramRecord.setProgramName(program);
            samProgramRecord.setProgramVersion(version);

            log.log(Level.INFO, "Processing reads, writing to " + outputSamOrBamFile.getName() + " ...");
            try (SAMFileWriter samFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(samFileHeader, true, outputSamOrBamFile)){

                //loop over BED records
                for (GenomicLocation genomicLocation : genomicLocations){
                    if (commandLine.hasOption("V")) log.log(Level.INFO, "Inspecting region: " + genomicLocation + " ...");

                    //get reference sequence
                    ReferenceSequence referenceSequence = new ReferenceSequence(new GenomicLocation(genomicLocation.getContig(), genomicLocation.getStartPosition(), genomicLocation.getEndPosition()), referenceFasta, referenceFastaFai);
                    referenceSequence.populateReferenceSequence();

                    //get primer sequences
                    String upstreamPrimerSequence = referenceSequence.getReferenceSequence().substring(0, genomicLocation.getUpstreamPrimerLength());
                    String downstreamPrimerSequence = referenceSequence.getReferenceSequence().substring(referenceSequence.getReferenceSequence().length() - genomicLocation.getDownstreamPrimerLength());

                    if (commandLine.hasOption("V")){
                        log.log(Level.FINE, "Reference sequence: " + referenceSequence.getReferenceSequence());
                        log.log(Level.FINE, "Upstream primer: " + upstreamPrimerSequence);
                        log.log(Level.FINE, "Downstream primer: " + downstreamPrimerSequence);
                    }

                    //query alignments
                    SAMRecordIterator samRecordIterator = samReader.queryOverlapping(genomicLocation.getContig(), genomicLocation.getStartPosition(), genomicLocation.getEndPosition());

                    samRecordIterator.stream()
                            .filter(samRecord -> !samRecord.getReadUnmappedFlag())
                            .filter(samRecord -> !samRecord.getNotPrimaryAlignmentFlag())
                            .filter(samRecord -> !samRecord.getSupplementaryAlignmentFlag())
                            .filter(samRecord -> {

                                //calculate hamming distances
                                int upstreamPrimerHammingDist = Hamming.getHammingDistance(samRecord.getReadString().substring(0, upstreamPrimerSequence.length()), upstreamPrimerSequence);
                                int downstreamPrimerHammingDist = Hamming.getHammingDistance(samRecord.getReadString().substring(samRecord.getReadLength() - downstreamPrimerSequence.length()), downstreamPrimerSequence);

                                double upstreamPrimerSimilarity = (double) (upstreamPrimerSequence.length() - upstreamPrimerHammingDist) / upstreamPrimerSequence.length();
                                double downstreamPrimerSimilarity = (double) (downstreamPrimerSequence.length() - downstreamPrimerHammingDist) / downstreamPrimerSequence.length();

                                if (upstreamPrimerSimilarity > primerSimilarity && downstreamPrimerSimilarity > primerSimilarity){
                                    return true;
                                }

                                return false;

                            })
                            .forEach(samRecord -> {

                                if (samRecord.getCigar().getFirstCigarElement().getOperator().equals(CigarOperator.SOFT_CLIP) ||
                                        samRecord.getCigar().getLastCigarElement().getOperator().equals(CigarOperator.SOFT_CLIP)){

                                    try {

                                        PairwiseAligner pairwiseAligner = new PairwiseAligner(referenceSequence.getReferenceSequence(), samRecord.getReadString());
                                        pairwiseAligner.needlemanWunschAlignment(gapOpenPenalty, gapExtendPenalty);

                                        if (pairwiseAligner.getScore() > minScore){

                                            //adjust alignment
                                            String readGroup = samRecord.getStringAttribute("RG");
                                            samRecord.clearAttributes();

                                            samRecord.setAttribute("RG", readGroup);
                                            samRecord.setAttribute("AS", (int) Math.round(pairwiseAligner.getScore()));
                                            samRecord.setAttribute("CO", genomicLocation.getName());
                                            samRecord.setAttribute("XC", samRecord.getCigarString());
                                            samRecord.setCigar(pairwiseAligner.getCigar());
                                            samRecord.setAlignmentStart(genomicLocation.getStartPosition());

                                            samFileWriter.addAlignment(samRecord);
                                        }

                                    } catch (CompoundNotFoundException e){
                                        log.log(Level.SEVERE, "Could not perform pairwise alignment: " + e.getMessage());
                                        System.exit(-1);
                                    }

                                } else {
                                    String readGroup = samRecord.getStringAttribute("RG");
                                    int alignmentScore = samRecord.getIntegerAttribute("AS");
                                    samRecord.clearAttributes();

                                    samRecord.setAttribute("RG", readGroup);
                                    samRecord.setAttribute("AS", alignmentScore);
                                    samRecord.setAttribute("CO", genomicLocation.getName());

                                    samFileWriter.addAlignment(samRecord);
                                }

                            });

                    samRecordIterator.close();
                }

            }


        } catch (IOException e){
            log.log(Level.SEVERE, "Could not read BAM file: " + e.getMessage());
            System.exit(-1);
        }

    }

}
