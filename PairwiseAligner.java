package nhs.genetics.cardiff.framework;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;

import java.util.ArrayList;

/**
 * Wrapper around BioJava pairwise alignment tools
 *
 * @author  Matt Lyon
 * @version 1.0
 * @since   2016-10-13
 */
public class PairwiseAligner {

    private DNASequence referenceDNASequence, queryDNASequence;
    private PairwiseSequenceAligner<DNASequence, NucleotideCompound> alignment;

    public PairwiseAligner(String referenceString, String queryString) throws CompoundNotFoundException {
        this.referenceDNASequence = new DNASequence(referenceString, AmbiguityDNACompoundSet.getDNACompoundSet());
        this.queryDNASequence = new DNASequence(queryString, AmbiguityDNACompoundSet.getDNACompoundSet());
    }

    public void smithWatermanAlignment(int gapOpen, int gapExt){
        SimpleGapPenalty simpleGapPenalty = new SimpleGapPenalty(gapOpen, gapExt);

        alignment = Alignments.getPairwiseAligner(
                queryDNASequence,
                referenceDNASequence,
                Alignments.PairwiseSequenceAlignerType.LOCAL,
                simpleGapPenalty,
                SubstitutionMatrixHelper.getNuc4_4());
    }

    public void needlemanWunschAlignment(int gapOpen, int gapExt){
        SimpleGapPenalty simpleGapPenalty = new SimpleGapPenalty(gapOpen, gapExt);

        alignment = Alignments.getPairwiseAligner(
                queryDNASequence,
                referenceDNASequence,
                Alignments.PairwiseSequenceAlignerType.GLOBAL,
                simpleGapPenalty,
                SubstitutionMatrixHelper.getNuc4_4());
    }

    public Cigar getCigar(){
        ArrayList<CigarElement> cigarElements = new ArrayList<>();

        CigarOperator currentCigarOperator = null;
        Integer currentCigarLength = 0;

        //loop over reference alignment
        for (int t = 1; t < alignment.getPair().getTarget().getLength() + 1; t++){
            if (alignment.getPair().getQuery().getCompoundAt(t).toString().equals("-")){ //query deletion

                if (currentCigarOperator == null){
                    currentCigarLength = 1;
                    currentCigarOperator = CigarOperator.D;
                } else if (currentCigarOperator == CigarOperator.D){
                    currentCigarLength++;
                } else {
                    cigarElements.add(new CigarElement(currentCigarLength, currentCigarOperator));

                    currentCigarLength = 1;
                    currentCigarOperator = CigarOperator.D;
                }

            } else if (alignment.getPair().getTarget().getCompoundAt(t).toString().equals("-")){ //query insertion

                if (currentCigarOperator == null){
                    currentCigarLength = 1;
                    currentCigarOperator = CigarOperator.I;
                } else if (currentCigarOperator == CigarOperator.I){
                    currentCigarLength++;
                } else {
                    cigarElements.add(new CigarElement(currentCigarLength, currentCigarOperator));

                    currentCigarLength = 1;
                    currentCigarOperator = CigarOperator.I;
                }

            } else { //match or mismatch

                if (currentCigarOperator == null){
                    currentCigarLength = 1;
                    currentCigarOperator = CigarOperator.M;
                } else if (currentCigarOperator == CigarOperator.M){
                    currentCigarLength++;
                } else {
                    cigarElements.add(new CigarElement(currentCigarLength, currentCigarOperator));

                    currentCigarLength = 1;
                    currentCigarOperator = CigarOperator.M;
                }

            }
        }

        //push last cigar
        cigarElements.add(new CigarElement(currentCigarLength, currentCigarOperator));

        return new Cigar(cigarElements);
    }

    public PairwiseSequenceAligner<DNASequence, NucleotideCompound> getAlignment() {
        return alignment;
    }
    public double getDistance(){
        return alignment.getDistance();
    }
    public double getScore(){
        return alignment.getScore();
    }

}
