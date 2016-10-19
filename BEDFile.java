package nhs.genetics.cardiff.framework;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import static nhs.genetics.cardiff.framework.GenomicLocation.contigComparator;

/**
 * Class for parsing BED files
 *
 * @author  Matt Lyon
 * @version 1.0
 * @since   2016-10-13
 */
public class BEDFile {
    public static ArrayList<GenomicLocation> getBedFeatures(File filePath) throws IOException {

        HashMap<String, GenomicLocation> amplicons = new HashMap<>();
        ArrayList<GenomicLocation> genomicLocations = new ArrayList<>();
        String line;

        try (BufferedReader bedReader = new BufferedReader(new FileReader(filePath))){
            while((line = bedReader.readLine())!= null) {

                if (!line.equals("")){
                    String[] fields = line.split("\t");

                    String contig = fields[0];
                    int startPosition = Integer.parseInt(fields[1]) + 1;
                    int endPosition = Integer.parseInt(fields[2]);
                    int thickStartPosition = Integer.parseInt(fields[6]) + 1;
                    int thickEndPosition = Integer.parseInt(fields[7]);
                    int upstreamPrimerLength = thickStartPosition - startPosition;
                    int downstreamPrimerLength = endPosition - thickEndPosition;
                    String name = fields[3];
                    String key = contig + ":" + startPosition + "-" + endPosition;

                    GenomicLocation genomicLocation = new GenomicLocation(contig, startPosition, endPosition, name);
                    genomicLocation.setUpstreamPrimerLength(upstreamPrimerLength);
                    genomicLocation.setDownstreamPrimerLength(downstreamPrimerLength);

                    if (amplicons.containsKey(key)){

                        //record longest primers for amplicons with the same start-stop coordinates

                        if (amplicons.get(key).getUpstreamPrimerLength() < upstreamPrimerLength){
                            amplicons.get(key).setUpstreamPrimerLength(upstreamPrimerLength);
                        }

                        if (amplicons.get(key).getDownstreamPrimerLength() < downstreamPrimerLength){
                            amplicons.get(key).setDownstreamPrimerLength(downstreamPrimerLength);
                        }

                        amplicons.get(key).setName(amplicons.get(key).getName() + "_" + name);

                    } else {
                        amplicons.put(key, genomicLocation);
                    }

                }

            }

        }

        //collect unique targets
        for (Map.Entry<String, GenomicLocation> iter : amplicons.entrySet()){
            genomicLocations.add(iter.getValue());
        }

        //sort by chromosome
        genomicLocations.sort(contigComparator);

        return genomicLocations;
    }
}
