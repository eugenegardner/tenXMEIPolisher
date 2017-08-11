package tenX.MEIBuilder;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava3.core.sequence.DNASequence;

import MELT.MELTIllumina.classification.lineU.RunAlignment;
import jaligner.matrix.MatrixLoaderException;
import tenX.MEISorter.BXIndexReader;

public class MEIBuilder {

	public static void main (String args[]) throws IOException, MatrixLoaderException, InterruptedException, ExecutionException {
				
		/*
		 * 0 = Original BAM
		 * 1 = BXTag index
		 * 2 = List of FL-L1s to assess (still working on format)
		 * 
		 */
		
		Map<String,Long> bxIndex = BXIndexReader.buildBXIndex(new File(args[0] + ".tagSorted.sam.idx"));
		bxIndex = BXIndexReader.filterBXIndex(bxIndex, new File(args[1]));
		
		BufferedReader pacReader = new BufferedReader(new FileReader(new File(args[2])));
		File tmpDir = new File(args[0] + ".tagSorted.sam.idx.tmp");
		if (!tmpDir.exists()) {
			tmpDir.mkdir();
		}
		String line;
		Pattern endDel = Pattern.compile("(\\S+)\\|d\\d+\\-\\d+");
		Pattern begDel = Pattern.compile("d\\d+\\-\\d+\\|(\\S+)");
		
		while ((line = pacReader.readLine()) != null) {
		
			String data[] = line.split("\t");
			String chr = data[0];
			int position = Integer.parseInt(data[1]);
			int start = position - 500;
			int end = position + 500;
			String originalSeq = getOrientation(data[2].toUpperCase());
			if (originalSeq == null) { //Twin Prime -- Check both sequences
				boolean isForward = isForward(data[2].substring(data[2].length() - 100, data[2].length()));
				PairAssessment assessForward = new PairAssessment(new File(args[0]), new File(args[0] + ".bai"), chr, position, start, end, bxIndex, data[2], tmpDir);
				PairAssessment assessReverse = new PairAssessment(new File(args[0]), new File(args[0] + ".bai"), chr, position, start, end, bxIndex, new DNASequence(data[2]).getReverseComplement().getSequenceAsString(), tmpDir); //Reverse Strand
				int totalFail = assessForward.getTotalFail() + assessReverse.getTotalFail();
				int totalFinalFail = assessForward.getTotalFinalFail() + assessReverse.getTotalFinalFail();
				int max = Math.max(assessForward.getIterations(), assessReverse.getIterations());
				String finalSeq;
				String finalChanges;
				if (isForward) {
					finalSeq = new DNASequence(assessReverse.getFinalSeq()).getReverseComplement().getSequenceAsString() + assessForward.getFinalSeq();
					Matcher endMatch = endDel.matcher(assessReverse.getLineNewChanges());
					Matcher begMatch = begDel.matcher(assessForward.getLineChanges());
					if (endMatch.matches() && begMatch.matches()) {
						finalChanges = endMatch.group(1) + "|INV|" + begMatch.group(1);
					} else {
						finalChanges = null;
					}
					finalSeq = assessForward.getFinalSeq(); //THIS IS TO ONLY OUTPUT UNINVERTED SEQ FOR ALIGNMENT PURPOSES! REMOVE WHEN DONE!
				} else {
					finalSeq = new DNASequence(assessForward.getFinalSeq()).getReverseComplement().getSequenceAsString() + assessReverse.getFinalSeq();
					Matcher endMatch = endDel.matcher(assessForward.getLineNewChanges());
					Matcher begMatch = begDel.matcher(assessReverse.getLineChanges());
					if (endMatch.matches() && begMatch.matches()) {
						finalChanges = endMatch.group(1) + "|INV|" + begMatch.group(1);
					} else {
						finalChanges = null;
					}
					finalSeq = assessReverse.getFinalSeq(); //THIS IS TO ONLY OUTPUT UNINVERTED SEQ FOR ALIGNMENT PURPOSES! REMOVE WHEN DONE!
				}
				
				System.out.println(chr + "\t" + position + "\t" + totalFail + "\t" + totalFinalFail + "\t" + finalSeq.length() + "\t" + max + "\ttrue\t" + finalChanges + "\t" + finalSeq);
			} else {
				PairAssessment assess = new PairAssessment(new File(args[0]), new File(args[0] + ".bai"), chr, position, start, end, bxIndex, originalSeq, tmpDir);
				System.out.println(chr + "\t" + position + "\t" + assess.getTotalFail() + "\t" + assess.getTotalFinalFail() + "\t" + assess.getFinalSeq().length() + "\t" + assess.getIterations() + "\tfalse\t" + assess.getLineNewChanges() + "\t" + assess.getFinalSeq());
			}
				
			
		}
		
		pacReader.close();
				
		
	}
	
	private static String getOrientation(String originalSeq) throws MatrixLoaderException {
		
		String revSeq = new DNASequence(originalSeq).getReverseComplement().getSequenceAsString();
		RunAlignment forwardAlignment = new RunAlignment(originalSeq);
		RunAlignment reverseAlignment = new RunAlignment(revSeq);
		
		float forwardScore = forwardAlignment.getScore();
		float reverseScore = reverseAlignment.getScore();

		if (forwardScore > 500 && reverseScore > 500) {
			return null;
		} else if (forwardScore > reverseScore) {
			return originalSeq;
		} else {
			return revSeq;
		}
		
	}
	private static boolean isForward(String testSeq) throws MatrixLoaderException {
		
		RunAlignment forwardAlignment = new RunAlignment(testSeq);
		
		int startPos = forwardAlignment.getStart();

		if (startPos > 5900) {
			return true;
		} else {
			return false;
		}
		
	}
	
}
