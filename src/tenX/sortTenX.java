package tenX;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileHeader.SortOrder;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

public class sortTenX {

	public sortTenX() {
		super();
	}
	
	private static Map<String, Integer> bxCount;
	
	public static File makeTagBam (SamReader finalReader, File revBam) {
		
		bxCount = new HashMap<String, Integer>();
		SAMRecordIterator itr = finalReader.iterator();
		SAMFileHeader reverseHeader = finalReader.getFileHeader();
		reverseHeader.setSortOrder(SortOrder.tenX);
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory().setMaxRecordsInRam(3000000);
		SAMFileWriter samWriter = writerFactory.makeSAMWriter(reverseHeader, false, revBam);
		while (itr.hasNext()) {
			SAMRecord rec = itr.next();
			String bxTag = rec.getStringAttribute("BX");
			if (bxTag != null) {
				samWriter.addAlignment(rec);
				addTag(bxTag);
			}
		}
		samWriter.close();
		itr.close();	
		spillTags();
		return revBam;
	}
	
	private static void addTag (String tag) {
		
		if (bxCount.containsKey(tag)) {
			int total = bxCount.get(tag);
			total++;
			bxCount.put(tag, total);
		} else {
			bxCount.put(tag, 1);
		}
		
	}
	
	private static void spillTags() {
		
		for (Map.Entry<String, Integer> entry : bxCount.entrySet()) {
			
			System.out.println(entry.getKey() + "\t" + entry.getValue());
			
		}
		
	}
	
}
