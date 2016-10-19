package nhs.genetics.cardiff.framework;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.logging.Logger;

/**
 * A class for representing genomic regions
 *
 * @author  Matt Lyon
 * @version 1.0
 * @since   2015-02-05
 */
public class GenomicLocation {

    private static final Logger log = Logger.getLogger(GenomicLocation.class.getName());

    private String contig, name;
    private int startPosition, endPosition, strand, upstreamPrimerLength, downstreamPrimerLength;

    public GenomicLocation(String contig, int startPosition, int endPosition){
        this.contig = contig;
        this.startPosition = startPosition;
        this.endPosition = endPosition;
    }
    public GenomicLocation(String contig, int startPosition, int endPosition, String name){
        this.contig = contig;
        this.startPosition = startPosition;
        this.endPosition = endPosition;
        this.name = name;
    }
    public GenomicLocation(String contig, int startPosition, int endPosition, String name, int strand){
        this.contig = contig;
        this.startPosition = startPosition;
        this.endPosition = endPosition;
        this.name = name;
        this.strand = strand;
    }
    public GenomicLocation(String contig, int startPosition){
        this.contig = contig;
        this.startPosition = startPosition;
        this.endPosition = endPosition;
        this.name = name;
    }

    public static boolean isInteger(String s) {
        try {
            Integer.parseInt(s);
        } catch(NumberFormatException e) {
            return false;
        }
        return true;
    }

    public static ArrayList<GenomicLocation> mergeBases(ArrayList<GenomicLocation> bases) {

        ArrayList<GenomicLocation> regions = new ArrayList<>();
        int lastStartPos = 0, lastEndPos = 0;
        String lastStartChrom = "";

        //sort bases by chrom pos
        Collections.sort(bases, GenomicLocation.contigComparator);

        //make windows of overlapping regions
        for (int j = 0; j < bases.size(); ++j){ //loop over bases

            if (j == 0) { //first base

                lastStartChrom = bases.get(j).getContig();
                lastStartPos = bases.get(j).getStartPosition();
                lastEndPos = bases.get(j).getStartPosition();

            } else if (j == (bases.size() - 1)) { //last bases

                lastEndPos++;
                regions.add(new GenomicLocation(lastStartChrom, lastStartPos, lastEndPos));

            } else { //not first or last base

                if (!bases.get(j).getContig().equals(lastStartChrom)) { //different contig -- start new window

                    regions.add(new GenomicLocation(lastStartChrom, lastStartPos, lastEndPos));

                    lastStartChrom = bases.get(j).getContig();
                    lastStartPos = bases.get(j).getStartPosition();
                    lastEndPos = bases.get(j).getStartPosition();

                } else if (bases.get(j).getStartPosition() != lastEndPos + 1) { //different position -- start new window

                    regions.add(new GenomicLocation(lastStartChrom, lastStartPos, lastEndPos));

                    lastStartPos = bases.get(j).getStartPosition();
                    lastEndPos = bases.get(j).getStartPosition();

                } else {
                    lastEndPos++;
                }

            }

        }

        return regions;
    }

    public void setLeftPadding(int noBases){
        if (startPosition > noBases) {
            startPosition = (startPosition - noBases);
        } else {
            startPosition = 1;
        }
    }

    public void convertTo1Based(){
        startPosition += 1;
    }

    public void convertTo0Based(){
        startPosition -= 1;
    }

    public String getContig() {
        return contig;
    }
    public String getName() {
        return name;
    }
    public int getStartPosition() {
        return startPosition;
    }
    public int getEndPosition() {
        return endPosition;
    }
    public int getStrand() {
        return strand;
    }
    public String getConcatenatedLocation(){
        return contig + ":" + startPosition + "-" + endPosition;
    }
    public String getConcatenatedPosition(){
        return contig + ":" + startPosition;
    }
    public int getUpstreamPrimerLength() {
        return upstreamPrimerLength;
    }
    public int getDownstreamPrimerLength() {
        return downstreamPrimerLength;
    }

    public void setRightPadding(int noBases){
        endPosition = (endPosition + noBases);
    }
    public void setContig(String contig) {
        this.contig = contig;
    }
    public void setName(String name) {
        this.name = name;
    }
    public void setStartPosition(int startPosition) {
        this.startPosition = startPosition;
    }
    public void setEndPosition(int endPosition) {
        this.endPosition = endPosition;
    }
    public void setStrand(int strand) {
        this.strand = strand;
    }
    public void setUpstreamPrimerLength(int upstreamPrimerLength) {
        this.upstreamPrimerLength = upstreamPrimerLength;
    }
    public void setDownstreamPrimerLength(int downstreamPrimerLength) {
        this.downstreamPrimerLength = downstreamPrimerLength;
    }

    /*Comparator for sorting the list by contig and position*/
    public static Comparator<GenomicLocation> contigComparator = new Comparator<GenomicLocation>() {
        public int compare(GenomicLocation g1, GenomicLocation g2) {
            int returnVal = 0;

            if (isInteger(g1.getContig()) && isInteger(g2.getContig())){ //both contigs are numeric

                if (Integer.parseInt(g1.getContig()) > Integer.parseInt(g2.getContig())) {
                    returnVal = 1;
                } else if (Integer.parseInt(g1.getContig()) < Integer.parseInt(g2.getContig())) {
                    returnVal =  -1;
                } else if (g1.getStartPosition() > g2.getStartPosition()){
                    returnVal =  1;
                } else if (g1.getStartPosition() < g2.getStartPosition()){
                    returnVal =  -1;
                }

            } else {

                if (g1.getContig().compareTo(g2.getContig()) != 0){ //contigs are string
                    returnVal = g1.getContig().compareTo(g2.getContig());
                } else if (g1.getStartPosition() > g2.getStartPosition()){
                    returnVal =  1;
                } else if (g1.getStartPosition() < g2.getStartPosition()){
                    returnVal =  -1;
                }

            }

            return returnVal;
        }
    };

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        GenomicLocation that = (GenomicLocation) o;

        if (startPosition != that.startPosition) return false;
        if (endPosition != that.endPosition) return false;
        return contig.equals(that.contig);

    }

    @Override
    public int hashCode() {
        int result = contig.hashCode();
        result = 31 * result + startPosition;
        result = 31 * result + endPosition;
        return result;
    }

    @Override
    public String toString(){
        return contig + ":" + startPosition + "-" + endPosition;
    }
}
