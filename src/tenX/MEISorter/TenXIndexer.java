package tenX.MEISorter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import tenX.MEISorter.BXTag.Read;

public class TenXIndexer {

	public static void main (String args[]) throws IOException {
		
		SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		SamReader tenXreader = readerFactory.open(SamInputResource.of(args[0]).index(new File(args[0] + ".bai")));
		
		SAMRecordIterator samItr = tenXreader.iterator();
		
		Map<String, BXTag> bxTags = new HashMap<String, BXTag>();
		
		int totalRecords = 0;
				
		//Iterate through and build index:
		while (samItr.hasNext() == true) {
			
			SAMRecord current = samItr.next();
			totalRecords++;
			
			if (current.getReadPairedFlag()) {
				if (!current.getMateUnmappedFlag() && !current.getNotPrimaryAlignmentFlag()) {
					String bxTag = current.getStringAttribute("BX");
					if (bxTag != null) {
						if (bxTags.containsKey(bxTag)) {
							bxTags.get(bxTag).addRead(current.getReferenceName(), current.getAlignmentStart());
						} else {
							bxTags.put(bxTag, new BXTag(current.getReferenceName(), current.getAlignmentStart()));
						}
					}
				}
			}
			if ((totalRecords % 10000000) == 0) {
				System.out.println(totalRecords + " reads processed" + "\t" + System.currentTimeMillis());
			}
		}
		
		// 1.idx = main index of all coordinates
		// 2.idx = index of file pointers to specific barcodes
		
		RandomAccessFile idx1Writer = new RandomAccessFile(new File(args[0] + ".1.idx"), "rw");
		BufferedWriter idx2Writer = new BufferedWriter(new FileWriter(new File(args[0] + ".2.idx")));
		long currentFilePointer = idx1Writer.getFilePointer();
		
		for (Map.Entry<String, BXTag> entry : bxTags.entrySet()) {
			
			List<Read> reads = entry.getValue().getReads();
			if (reads.size() > 15 && reads.size() < 4000) {
				idx1Writer.writeBytes("+" + entry.getKey() + "\n");
				currentFilePointer = idx1Writer.getFilePointer();
				idx2Writer.write(entry.getKey() + "\t" + currentFilePointer + "\n");
				for (Read read : reads) {
					idx1Writer.writeBytes(read.getChr() + "\t" + read.getStart() + "\n");
				}
			}
		}
		idx1Writer.close();
		idx2Writer.close();
		
	}
		
}
