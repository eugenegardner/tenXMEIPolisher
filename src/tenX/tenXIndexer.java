package tenX;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class tenXIndexer {

	public static void main (String args[]) throws IOException {
		
		/* File I/O:
		 * 
		 * 0 = SAMFile w/o header
		 * 
		 * Index is printed to args[0].idx
		 * 
		 */ 
		
		//Get a random Access File to build Index
		RandomAccessFile random = new RandomAccessFile(new File(args[0]), "r");
		//Iterate through file line by line grabbing the first instance of each sequence family
		Map<String, Long> tags = new LinkedHashMap<String, Long>();
		String line;
		long lastFilePointer = 0;
		int totalRecords = 0;
		BufferedWriter indexWriter = new BufferedWriter(new FileWriter(new File(args[0] + ".idx")));
		Pattern BXpatt = Pattern.compile("[\\S\\s]+BX:Z:([ATCG]+\\-1)[\\S\\s]+");

		while ((line = random.readLine()) != null) {
			
			Matcher BXmatch = BXpatt.matcher(line);
			if (line.matches("@[\\S\\s]+")) {
				lastFilePointer = random.getFilePointer();
				continue;
			} else if (BXmatch.matches()) {
				totalRecords++;
				String bxTag = BXmatch.group(1);
				if (!tags.containsKey(bxTag)) {
					indexWriter.write(bxTag + "\t" + lastFilePointer + "\n");
					indexWriter.flush();
					tags.put(bxTag, lastFilePointer);
				}
				if ((totalRecords % 10000000) == 0) {
					System.out.println(totalRecords + " reads processed" + "\t" + System.currentTimeMillis());
				}
			} else {
				System.err.println("BXTAG DID NOT MATCH FOR " + line);
				System.exit(1);
			}
			lastFilePointer = random.getFilePointer();
						
		}
		indexWriter.close();
		random.close();
		
	}
	
}
