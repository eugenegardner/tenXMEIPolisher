package tenX.MEISorter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class BXIndexReader {

	public BXIndexReader() {
		super();
	}
	
	public static Map<String, Long> buildBXIndex (File bxIndexFile) throws IOException {
		
		Map<String, Long> bxIndex = new HashMap<String, Long>();
		
		BufferedReader idx2Reader = new BufferedReader(new FileReader(bxIndexFile));
		
		String line;
		String data[];
		
		while ((line = idx2Reader.readLine()) != null) {
			
			data = line.split("\t");
			bxIndex.put(data[0], Long.parseLong(data[1]));
			
		}
		
		idx2Reader.close();
		return bxIndex;
		
	}
	
	public static Map<String, Long> filterBXIndex (Map<String, Long> oldIndex, File tags) throws IOException {
		
		Map<String, Long> bxIndex = new HashMap<String, Long>();
		Map<String, Integer> filter = new HashMap<String, Integer>();
		
		DescriptiveStatistics stat = new DescriptiveStatistics();
		
		BufferedReader filterReader = new BufferedReader(new FileReader(tags));
		
		String line;
		String data[];
		double min = 15;
		
		while ((line = filterReader.readLine()) != null) {
			
			data = line.split("\t");
			int totalReads = Integer.parseInt(data[1]);
			filter.put(data[0], Integer.parseInt(data[1]));
			if (totalReads > min) {
				stat.addValue(Double.parseDouble(data[1]));
			}
			
		}
		filterReader.close();
		
		double max = stat.getPercentile(50.0) * 5;
		
		for (Map.Entry<String, Long> oldEntry : oldIndex.entrySet()) {
			
			if (filter.containsKey(oldEntry.getKey())) {
				
				int totalReads = filter.get(oldEntry.getKey());
				if (totalReads > min && totalReads < max) {
					bxIndex.put(oldEntry.getKey(), oldEntry.getValue());
				}
				
			} else {
				System.out.println("MISSING ENTRY : " + oldEntry.getKey());
			}
			
		}
		
		return bxIndex;
		
		
	}
	
}
