package tenX.MEIBuilder;

import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Hashtable;

public class Utilities {

	public Utilities () {
		super();
	}
	
	//Check Quality String
	public Hashtable<Character,Integer> loadQualities() throws IOException {
		
		BufferedReader qualityStream = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/tenX/MEIBuilder/importFiles/qualities.txt"), "UTF-8"));
		Hashtable<Character,Integer> qualityTable = new Hashtable<Character,Integer>();
		String myLine;
		String data[];
		
		try {
			while((myLine = qualityStream.readLine()) != null) {
				data = myLine.split("\t");
				qualityTable.put(data[0].charAt(0), Integer.parseInt(data[1]));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		qualityStream.close();
		return qualityTable;
		
	}
	public double qualityString(String errorString, Hashtable<Character,Integer> qualityTable) {
		
		double quality = 0;
		double totalQuals = 0;
		char charArray[] = errorString.toCharArray();
		
		for (char itr : charArray) {
			totalQuals++;
			quality+=qualityTable.get(itr);
		}

		quality/=totalQuals;
		return quality;
	}
	
	//Check Alignment
	public AlignmentDecision alignmentMatrix(Sequence querySeq, Sequence queryRevSeq, Sequence refSeq, double forLength, double queryQual, Matrix matrix) {
		
		AlignmentDecision decision;
		float openGap = 50;
		float extendPen = 5;
		Alignment forwardAlignment = SmithWatermanGotoh.align(querySeq, refSeq, matrix, openGap, extendPen);
		double forIdentity = forwardAlignment.getIdentity();
		double forPerc = (forIdentity / forLength) * 100;
		
//		System.out.println(String.valueOf(forwardAlignment.getSequence1()) + "\n" +
//				String.valueOf(forwardAlignment.getMarkupLine()) + "\n" +
//				String.valueOf(forwardAlignment.getSequence2()) + "\n" +
//				forPerc + "\n" + querySeq.getSequence());
		
		if (forPerc >= queryQual) {
			
			decision = AlignmentDecision.FOR;
			
		} else {
			
			Alignment reverseAlignment = SmithWatermanGotoh.align(queryRevSeq, refSeq, matrix, openGap, extendPen);
			double revIdentity = reverseAlignment.getIdentity();
			double revPerc = (revIdentity / forLength) * 100;
			
//			System.out.println(String.valueOf(reverseAlignment.getSequence1()) + "\n" +
//					String.valueOf(reverseAlignment.getMarkupLine()) + "\n" +
//					String.valueOf(reverseAlignment.getSequence2()) + "\n" +
//					revPerc + "\n" + queryRevSeq.getSequence());
			
			if (revPerc >= queryQual) {

				decision = AlignmentDecision.REV;
	
			} else {
				
				decision = AlignmentDecision.NONE;
				
			}
			
		}
		
		return decision;
	
	}
	public enum AlignmentDecision {
		FOR,REV,NONE;
	}
	
	//PolA caller
	public polADecision callPolyA (Sequence querySeq) {
		
		char queries[] = {'A','T'};
		String subString;
		String revSubString;
		String querySeqString = querySeq.getSequence();
		char querySeqChar[] = querySeq.getSequence().toCharArray();	
		polADecision Decision = new polADecision();
		
		for (char query : queries) {	
			
			int Score = 5;
			int droppedHere = 1;
			boolean triggered = false;
			int currentPos = 1;
			int repeatPlus = 0;
			
			for (char base : querySeqChar) {
				
				if (base != query && triggered == false) {
					if (repeatPlus >= 2) {
						droppedHere = currentPos;
					}
					repeatPlus = 0;
					triggered = true;
					Score--;
					if (Score <= 0) {break;}
				} else if (base != query && triggered == true) {
					repeatPlus = 0;
					Score--;
					if (Score <= 0) {break;}				
				} else if (currentPos == querySeqChar.length) { 
					droppedHere = currentPos + 1;
					break;
				} else {
					repeatPlus++;
					triggered = false;
				}
	
				currentPos++;
				
			}
		
			//Check if result is at least 80% the char we are searching (eliminates things that didn't start with T/A
			subString = querySeqString.substring(0, (droppedHere - 1));
			revSubString = querySeqString.substring(droppedHere - 1);
			
			char checkArray[] = subString.toCharArray();
			
			double checkTot = 0;
			double checkSize = checkArray.length;

			for (char check : checkArray) {
				
				if (check == query) {
					
					checkTot++;
					
				}
				
			}
			
			// If matches, returns positive (1) for match, and the Alignment without PolyA Tail
			// If Doesn't match, returns negative (2) for match, and a blank substring
			
			if ((checkTot / checkSize) >= .70 && checkSize > 5) {

				Decision.setPolDec(true);
				Decision.setSubString(revSubString);
				Decision.setPolALength(checkSize);
				break;
				
			} else {
				
				Decision.setPolDec(false);
				Decision.setSubString("");
				Decision.setPolALength(0);
				
			}
			
		}	
		
		return Decision;
		
	}
	public class polADecision {
		private boolean polDec;
		private String subString;
		private double polALength;
		
		public boolean getPolDec() {
			return polDec;
		}
		private void setPolDec(boolean polDec) {
			this.polDec = polDec;
		}
		public String getSubString() {
			return subString;
		}
		private void setSubString(String subString) {
			this.subString = subString;
		}
		public double getPolALength() {
			return polALength;
		}
		public void setPolALength(double polALength) {
			this.polALength = polALength;
		}
		
		
	}
	
}
