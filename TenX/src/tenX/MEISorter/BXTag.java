package tenX.MEISorter;

import java.util.ArrayList;
import java.util.List;

public class BXTag {

	private List<Read> reads;
	
	public BXTag (String chr, int start) {
		
		reads = new ArrayList<Read>();
		addRead(chr, start);
		
	}
	public BXTag () {
		
		reads = new ArrayList<Read>();
		
	}
	
	public void addRead (String chr, int start) {
		
		Read currentReader = new Read(chr, start);
		reads.add(currentReader);
		
	}
	
	public List<Read> getReads () {
		return reads;
	}
	
	public class Read {
		
		private String chr;
		private int start;
		
		public Read (String chr, int start) {
			
			this.chr = chr;
			this.start = start;
						
		}

		public String getChr() {
			return chr;
		}

		public int getStart() {
			return start;
		}
		
	}
	
}
