package tenX.MEIBuilder;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.ExecuteException;
import org.apache.commons.exec.PumpStreamHandler;

import MELT.MELTIllumina.classification.lineU.LineResults;
import MELT.MELTIllumina.classification.lineU.LineU;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMLineParser;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import jaligner.Sequence;
import jaligner.matrix.MatrixLoaderException;

public class PairAssessment {
	
	private String chr;
	private int pos;
	private int start;
	private int end;
	private File tenXBam;
	private File tenXIndex;
	private String originalSeq;
	private File tmpDir;
	
	private Sequence meiLoaded;
	
	private Map<String, Long> index;
	private SamReaderFactory readerFactory;
	private Map<Integer, Change> changes;
	
	private SamReader tenXReader;
	
	private DefaultExecutor executor;
	
	private String finalSeq;
	private String lineChanges;
	private String lineNewChanges;
	private int totalFail;
	private int totalFinalFail;
	private int iterations;
	
	//Results list:
	Map<String, Integer> positiveBXTags;
	
	public PairAssessment (File tenXBam, File tenXIndex, String chr, int pos, int start, int end, Map<String, Long> index, String originalSeq, File tmpDir) throws IOException, MatrixLoaderException, InterruptedException, ExecutionException {
		
		this.tenXBam = tenXBam;
		this.tenXIndex = tenXIndex;
		this.chr = chr;
		this.pos = pos;
		this.start = start;
		this.end = end;
		this.index = index;
		this.originalSeq = originalSeq.toUpperCase();
		this.tmpDir = tmpDir;
		
		positiveBXTags = new HashMap<String, Integer>();
		executor = new DefaultExecutor();
						
		meiLoaded = new Sequence("GGGGGAGGAGCCAAGATGGCCGAATAGGAACAGCTCCGGTCTACAGCTCCCAGCGTGAGCGACGCAGAAGACGGTGATTTCTGCATTTCCATCTGAGGTACCGGGTTCATCTCACTAGGGAGTGCCAGACAGTGGGCGCAGGCCAGTGTGTGTGCGCACCGTGCGCGAGCCGAAGCAGGGCGAGGCATTGCCTCACCTGGGAAGCGCAAGGGGTCAGGGAGTTCCCTTTCTGAGTCAAAGAAAGGGGTGACGGTCGCACCTGGAAAATCGGGTCACTCCCACCCGAATATTGCGCTTTTCAGACCGGCTTAAGAAACGGCGCACCACGAGACTATATCCCACACCTGGCTCGGAGGGTCCTACGCCCACGGAATCTCGCTGATTGCTAGCACAGCAGTCTGAGATCAAACTGCAAGGCGGCAACGAGGCTGGGGGAGGGGCGCCCGCCATTGCCCAGGCTTGCTTAGGTAAACAAAGCAGCCGGGAAGCTCGAACTGGGTGGAGCCCACCACAGCTCAAGGAGGCCTGCCTGCCTCTGTAGGCTCCACCTCTGGGGGCAGGGCACAGACAAACAAAAAGACAGCAGTAACCTCTGCAGACTTAAGTGTCCCTGTCTGACAGCTTTGAAGAGAGCAGTGGTTCTCCCAGCACGCAGCTGGAGATCTGAGAACGGGCAGACAGACTGCCTCCTCAAGTGGGTCCCTGACTCCTGACCCCCGAGCAGCCTAACTGGGAGGCACCCCCCAGCAGGGGCACACTGACACCTCACACGGCAGGGTATTCCAACAGACCTGCAGCTGAGGGTCCTGTCTGTTAGAAGGAAAACTAACAACCAGAAAGGACATCTACACCGAAAACCCATCTGTACATCACCATCATCAAAGACCAAAAGTAGATAAAACCACAAAGATGGGGAAAAAACAGAACAGAAAAACTGGAAACTCTAAAACGCAGAGCGCCTCTCCTCCTCCAAAGGAACGCAGTTCCTCACCAGCAACGGAACAAAGCTGGATGGAGAATGATTTTGACGAGCTGAGAGAAGAAGGCTTCAGACGATCAAATTACTCTGAGCTACGGGAGGACATTCAAACCAAAGGCAAAGAAGTTGAAAACTTTGAAAAAAATTTAGAAGAATGTATAACTAGAATAACCAATACAGAGAAGTGCTTAAAGGAGCTGATGGAGCTGAAAACCAAGGCTCGAGAACTACGTGAAGAATGCAGAAGCCTCAGGAGCCGATGCGATCAACTGGAAGAAAGGGTATCAGCAATGGAAGATGAAATGAATGAAATGAAGCGAGAAGGGAAGTTTAGAGAAAAAAGAATAAAAAGAAATGAGCAAAGCCTCCAAGAAATATGGGACTATGTGAAAAGACCAAATCTACGTCTGATTGGTGTACCTGAAAGTGATGTGGAGAATGGAACCAAGTTGGAAAACACTCTGCAGGATATTATCCAGGAGAACTTCCCCAATCTAGCAAGGCAGGCCAACGTTCAGATTCAGGAAATACAGAGAACGCCACAAAGATACTCCTCGAGAAGAGCAACTCCAAGACACATAATTGTCAGATTCACCAAAGTTGAAATGAAGGAAAAAATGTTAAGGGCAGCCAGAGAGAAAGGTCGGGTTACCCTCAAAGGAAAGCCCATCAGACTAACAGTGGATCTCTCGGCAGAAACCCTACAAGCCAGAAGAGAGTGGGGGCCAATATTCAACATTCTTAAAGAAAAGAATTTTCAACCCAGAATTTCATATCCAGCCAAACTAAGCTTCATAAGTGAAGGAGAAATAAAATACTTTATAGACAAGCAAATGTTGAGAGATTTTGTCACCACCAGGCCTGCCCTAAAAGAGCTCCTGAAGGAAGCGCTAAACATGGAAAGGAACAACCGGTACCAGCCGCTGCAAAATCATGCCAAAATGTAAAGACCATCGAGACTAGGAAGAAACTGCATCAACTAATGAGCAAAATCACCAGCTAACATCATAATGACAGGATCAAATTCACACATAACAATATTAACTTTAAATATAAATGGACTAAATTCTGCAATTAAAAGACACAGACTGGCAAGTTGGATAAAGAGTCAAGACCCATCAGTGTGCTGTATTCAGGAAACCCATCTCACGTGCAGAGACACACATAGGCTCAAAATAAAAGGATGGAGGAAGATCTACCAAGCCAATGGAAAACAAAAAAAGGCAGGGGTTGCAATCCTAGTCTCTGATAAAACAGACTTTAAACCAACAAAGATCAAAAGAGACAAAGAAGGCCATTACATAATGGTAAAGGGATCAATTCAACAAGAGGAGCTAACTATCCTAAATATTTATGCACCCAATACAGGAGCACCCAGATTCATAAAGCAAGTCCTCAGTGACCTACAAAGAGACTTAGACTCCCACACATTAATAATGGGAGACTTTAACACCCCACTGTCAACATTAGACAGATCAACGAGACAGAAAGTCAACAAGGATACCCAGGAATTGAACTCAGCTCTGCACCAAGCAGACCTAATAGACATCTACAGAACTCTCCACCCCAAATCAACAGAATATACCTTTTTTTCAGCACCACACCACACCTATTCCAAAATTGACCACATAGTTGGAAGTAAAGCTCTCCTCAGCAAATGTAAAAGAACAGAAATTATAACAAACTATCTCTCAGACCACAGTGCAATCAAACTAGAACTCAGGATTAAGAATCTCACTCAAAGCCGCTCAACTACATGGAAACTGAACAACCTGCTCCTGAATGACTACTGGGTACATAACGAAATGAAGGCAGAAATAAAGATGTTCTTTGAAACCAACGAGAACAAAGACACCACATACCAGAATCTCTGGGACGCATTCAAAGCAGTGTGTAGAGGGAAATTTATAGCACTAAATGCCTACAAGAGAAAGCAGGAAAGATCCAAAATTGACACCCTAACATCACAATTAAAAGAACTAGAAAAGCAAGAGCAAACACATTCAAAAGCTAGCAGAAGGCAAGAAATAACTAAAATCAGAGCAGAACTGAAGGAAATAGAGACACAAAAAACCCTTCAAAAAATCAATGAATCCAGGAGCTGGTTTTTTGAAAGGATCAACAAAATTGATAGACCGCTAGCAAGACTAATAAAGAAAAAAAGAGAGAAGAATCAAATAGACACAATAAAAAATGATAAAGGGGATATCACCACCGATCCCACAGAAATACAAACTACCATCAGAGAATACTACAAACACCTCTACGCAAATAAACTAGAAAATCTAGAAGAAATGGATACATTCCTCGACACATACACTCTCCCAAGACTAAACCAGGAAGAAGTTGAATCTCTGAATAGACCAATAACAGGCTCTGAAATTGTGGCAATAATCAATAGTTTACCAACCAAAAAGAGTCCAGGACCAGATGGATTCACAGCCGAATTCTACCAGAGGTACATGGAGGAACTGGTACCATTCCTTCTGAAACTATTCCAATCAATAGAAAAAGAGGGAATCCTCCCTAACTCATTTTATGAGGCCAGCATCATTCTGATACCAAAGCCGGGCAGAGACACAACCAAAAAAGAGAATTTTAGACCAATATCCTTGATGAACATTGATGCAAAAATCCTCAATAAAATACTGGCAAACCGAATCCAGCAGCACATCAAAAAGCTTATCCACCATGATCAAGTGGGCTTCATCCCTGGGATGCAAGGCTGGTTCAATATACGCAAATCAATAAATGTAATCCAGCATATAAACAGAGCCAAAGACAAAAACCACATGATTATCTCAATAGATGCAGAAAAAGCCTTTGACAAAATTCAACAACCCTTCATGCTAAAAACTCTCAATAAATTAGGTATTGATGGGACGTATTTCAAAATAATAAGAGCTATCTATGACAAACCCACAGCCAATATCATACTGAATGGGCAAAAACTGGAAGCATTCCCTTTGAAAACCGGCACAAGACAGGGATGCCCTCTCTCACCGCTCCTATTCAACATAGTGTTGGAAGTTCTGGCCAGGGCAATCAGGCAGGAGAAGGAAATAAAGGGTATTCAATTAGGAAAAGAGGAAGTCAAATTGTCCCTGTTTGCAGACGACATGATTGTATATCTAGAAAACCCCATCGTCTCAGCCCAAAATCTCCTTAAGCTGATAAGCAACTTCAGCAAAGTCTCAGGATACAAAATCAATGTACAAAAATCACAAGCATTCTTATACACCAACAACAGACAAACAGAGAGCCAAATCATGGGTGAACTCCCATTCGTAATTGCTTCAAAGAGAATAAAATACCTAGGAATCCAACTTACAAGGGATGTGAAGGACCTCTTCAAGGAGAACTACAAACCACTGCTCAAGGAAATAAAAGAGGACACAAACAAATGGAAGAACATTCCATGCTCATGGGTAGGAAGAATCAATATCGTGAAAATGGCCATACTGCCCAAGGTAATTTACAGATTCAATGCCATCCCCATCAAGCTACCAATGACTTTCTTCACAGAATTGGAAAAAACTACTTTAAAGTTCATATGGAACCAAAAAAGAGCCCGCATTGCCAAGTCAATCCTAAGCCAAAAGAACAAAGCTGGAGGCATCACACTACCTGACTTCAAACTATACTACAAGGCTACAGTAACCAAAACAGCATGGTACTGGTACCAAAACAGAGATATAGATCAATGGAACAGAACAGAGCCCTCAGAAATAATGCCGCATATCTACAACTATCTGATCTTTGACAAACCTGAGAAAAACAAGCAATGGGGAAAGGATTCCCTATTTAATAAATGGTGCTGGGAAAACTGGCTAGCCATATGTAGAAAGCTGAAACTGGATCCCTTCCTTACACCTTATACAAAAATCAATTCAAGATGGATTAAAGATTTAAACGTTAAACCTAAAACCATAAAAACCCTAGAAGAAAACCTAGGCATTACCATTCAGGACATAGGCGTGGGCAAGGACTTCATGTCCAAAACACCAAAAGCAATGGCAACAAAAGACAAAATTGACAAATGGGATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTACCATCAGAGTGAACAGGCAACCTACAACATGGGAGAAAATTTTCGCAACCTACTCATCTGACAAAGGGCTAATATCCAGAATCTACAATGAACTTAAACAAATTTACAAGAAAAAAACAAACAACCCCATCAAAAAGTGGGCGAAGGACATGAACAGACACTTCTCAAAAGAAGACATTTATGCAGCCAAAAAACACATGAAGAAATGCTCATCATCACTGGCCATCAGAGAAATGCAAATCAAAACCACTATGAGATATCATCTCACACCAGTTAGAATGGCAATCATTAAAAAGTCAGGAAACAACAGGTGCTGGAGAGGATGCGGAGAAATAGGAACACTTTTACACTGTTGGTGGGACTGTAAACTAGTTCAACCATTGTGGAAGTCAGTGTGGCGATTCCTCAGGGATCTAGAACTAGAAATACCATTTGACCCAGCCATCCCATTACTGGGTATATACCCAAATGAGTATAAATCATGCTGCTATAAAGACACATGCACACGTATGTTTATTGCGGCACTATTCACAATAGCAAAGACTTGGAACCAACCCAAATGTCCAACAATGATAGACTGGATTAAGAAAATGTGGCACATATACACCATGGAATACTATGCAGCCATAAAAAATGATGAGTTCATATCCTTTGTAGGGACATGGATGAAATTGGAAACCATCATTCTCAGTAAACTATCGCAAGAACAAAAAACCAAACACCGCATATTCTCACTCATAGGTGGGAATTGAACAATGAGATCACATGGACACAGGAAGGGGAATATCACACTCTGGGGACTGTGGTGGGGTCGGGGGAGGGGGGAGGGATAGCATTGGGAGATATACCTAATGCTAGATGACACATTAGTGGGTGCAGCGCACCAGCATGGCACATGTATACATATGTAACTAACCTGCACAATGTGCACATGTACCCTAAAACTTAGAGTAT");
		doWork();
		
	}
	
	//Getters for information generated during run-time
	public String getFinalSeq() {
		return finalSeq;
	}
	public String getLineChanges() {
		return lineChanges;
	}
	public String getLineNewChanges() {
		return lineNewChanges;
	}
	public int getTotalFail() {
		return totalFail;
	}
	public int getTotalFinalFail() {
		return totalFinalFail;
	}
	public int getIterations() {
		return iterations;
	}
	
	//Method for running other methods in this class
	private void doWork() throws IOException, InterruptedException, ExecutionException, MatrixLoaderException {
		readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		tenXReader = readerFactory.open(SamInputResource.of(tenXBam).index(tenXIndex));
		checkReads();
		assessFinals();
		runBowtie();
		runLineU();
		Map<Integer, Integer> depth = countDepth();
		totalFail = 99999;
		iterations = 1;
		while (totalFail > 0 && iterations < 9) { 
			totalFail = countKmer(depth);
			iterations++;
		}
		//Final iteration to consider possible cluster of errors:
		if (totalFail > 0) {
			Map<Integer, Change> finalCheck = countKmerFinal(depth);
			iterations++;
			if (finalCheck != null) {
				changes = finalCheck;
			}
		}
		
		if (totalFail < 999) {
			buildFinal();
		} else {
			buildFinal();
			//Do something here for BAD sites???
		}
		tenXReader.close();
	}

	
	//Check reads for mapping to LINE reference
	private void checkReads() {
				
		SAMRecordIterator samItr = tenXReader.queryOverlapping(chr, start, end);
		while (samItr.hasNext()) {
				
			SAMRecord currentRecord = samItr.next();
			String bxTag = currentRecord.getStringAttribute("BX");
			if (currentRecord.getReadUnmappedFlag() == true || currentRecord.getMappingQuality() == 0) {
				continue;
			}
			
			if (bxTag != null) {
				if (positiveBXTags.containsKey(bxTag)) {
					int total = positiveBXTags.get(bxTag);
					total++;
					positiveBXTags.put(bxTag, total);
				} else {
					positiveBXTags.put(bxTag, 1);
				}
			}
			
//			String cigar = currentRecord.getCigarString();
//			Matcher right = cigarCaptureRight.matcher(cigar);
//			Matcher left = cigarCaptureLeft.matcher(cigar);
//			if (currentRecord.getMateReferenceName().equals(currentRecord.getReferenceName()) == false || currentRecord.getInferredInsertSize() > 1000 || currentRecord.getInferredInsertSize() < -1000) {
//
//				if (!(right.matches() && left.matches())) {
//					if (positiveBXTags.containsKey(bxTag)) {
//						int total = positiveBXTags.get(bxTag);
//						total++;
//						positiveBXTags.put(bxTag, total);
//					} else {
//
//						positiveBXTags.put(bxTag, 1);
//					}
//				}
//				
//			} else if (right.matches() ^ left.matches()) {
//				
//				if (catalogSplit(cigar, currentRecord)) {
//					if (positiveBXTags.containsKey(bxTag)) {
//						int total = positiveBXTags.get(bxTag);
//						total++;
//						positiveBXTags.put(bxTag, total);
//					} else {
//						positiveBXTags.put(bxTag, 1);
//					}
//				}
//									
//			}
			
		}
		samItr.close();
		
	}
	private void assessFinals() throws IOException, InterruptedException, ExecutionException, MatrixLoaderException {
		
		//Collect 'pass' records:
		Map<String, Pairs> validPairs = new HashMap<String, Pairs>();
		RandomAccessFile bxSortedSAM = new RandomAccessFile(new File(tenXBam.getAbsolutePath() + ".tagSorted.sam"), "r");
		ReadFetcher fetcher = new ReadFetcher(bxSortedSAM, new SAMLineParser(new DefaultSAMRecordFactory(), ValidationStringency.SILENT, tenXReader.getFileHeader(), null, null), index, chr, pos);
				
		for (Map.Entry<String, Integer> entry : positiveBXTags.entrySet()) {
			String bxTag = entry.getKey();
			int total = entry.getValue();
			if (total >= 1) { //Add check here some time for tags with more or less total tags at site?

				if (index.containsKey(bxTag)) {
				
					List<SAMRecord> records = fetcher.getRecords(bxTag);
					
					for (SAMRecord read : records) {
					
						if (validPairs.containsKey(read.getReadName())) {
							validPairs.get(read.getReadName()).parseRead(read);
						} else {
							Pairs pair = new Pairs(read);
							validPairs.put(read.getReadName(), pair);
						}
					}
				}
			}
		}
		
		bxSortedSAM.close();
		
		//Now dump all valid pairs:
		
		//Set FQ IO:
		FastqWriter fastqWriterFirst = new FastqWriterFactory().newWriter(new File(tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".1.fq"));
		FastqWriter fastqWriterSecond = new FastqWriterFactory().newWriter(new File(tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".2.fq"));
		FastqWriter fastqWriterUnpaired = new FastqWriterFactory().newWriter(new File(tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".up.fq"));
		String fastqSeperator = "+";
		
		for (Map.Entry<String, Pairs> entry : validPairs.entrySet()) {
			
			Pairs current = entry.getValue();
			//Write FQ in interleaved format:
			if (current.isBothPairs()) {
				fastqWriterFirst.write(new FastqRecord(current.getFirstPair().getReadName(), current.getFirstPair().getReadString(), fastqSeperator, current.getFirstPair().getBaseQualityString()));
				fastqWriterSecond.write(new FastqRecord(current.getSecondPair().getReadName(), current.getSecondPair().getReadString(), fastqSeperator, current.getSecondPair().getBaseQualityString()));
			} else {
				if (current.getFirstPair() != null) {
					fastqWriterUnpaired.write(new FastqRecord(current.getFirstPair().getReadName(), current.getFirstPair().getReadString(), fastqSeperator, current.getFirstPair().getBaseQualityString()));
				} else if (current.getSecondPair() != null) {
					fastqWriterUnpaired.write(new FastqRecord(current.getSecondPair().getReadName(), current.getSecondPair().getReadString(), fastqSeperator, current.getSecondPair().getBaseQualityString()));
				}
			}
			
		}
				
		fastqWriterFirst.close();
		fastqWriterSecond.close();
		fastqWriterUnpaired.close();
		
	}
	private void runBowtie () throws ExecuteException, IOException {
		
		String command = "bowtie2 --local --quiet -q -x /local/aberdeen2rw/DEVINE/SVU41/ILLUMINA/MERGE_ANALYSIS/me_refs/LINE1 -1 " + tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".1.fq -2 " + tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".2.fq -U " + tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".up.fq -S " + tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".sam";
		executor.execute(CommandLine.parse(command));
		SamReader testReader = readerFactory.open(new File(tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".sam"));
		SAMFileHeader newHeader = testReader.getFileHeader();
		newHeader.setSortOrder(SortOrder.coordinate);
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory().setCreateIndex(true);
		SAMFileWriter bamWriter = writerFactory.makeBAMWriter(newHeader, false, new File(tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".sorted.bam"));
		FastqWriter fastqWriterKmer = new FastqWriterFactory().newWriter(new File(tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".fq"));
		SAMRecordIterator itr = testReader.iterator();
		int allowedMis = 4;
		while (itr.hasNext()) {
			
			SAMRecord currentRecord = itr.next();
			boolean isUnMapped = currentRecord.getReadUnmappedFlag();
			
			if (isUnMapped == false) {
				
				int totalMis = currentRecord.getIntegerAttribute("XM");

				if (totalMis <= allowedMis) {
				
					bamWriter.addAlignment(currentRecord);
					fastqWriterKmer.write(new FastqRecord(currentRecord.getReadName(), currentRecord.getReadString(),"+",currentRecord.getBaseQualityString()));
				
				}
					
			}
			
		}
		
		itr.close();
		fastqWriterKmer.close();
		bamWriter.close();
				
	}
	
	//Count depth at each position in 10X data
	private Map<Integer, Integer> countDepth () throws ExecuteException, IOException {
	
		Map<Integer, Integer> depth = new HashMap<Integer, Integer>();
		
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PumpStreamHandler ps = new PumpStreamHandler(os);
		executor.setStreamHandler(ps);
		executor.execute(CommandLine.parse("bedtools genomecov -ibam " + tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".sorted.bam -g /local/aberdeen2rw/DEVINE/SVU41/ILLUMINA/MERGE_ANALYSIS/me_refs/LINE1.fa.fai -d"));
		executor.setStreamHandler(ps);
		BufferedReader r = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(os.toByteArray())));
		String line;
		
		while ((line = r.readLine()) != null) {
			String data[] = line.split("\t");
			depth.put(Integer.parseInt(data[1]), Integer.parseInt(data[2]));
		}
		r.close();
		return depth;
		
	}
	
	//Run LINEU on original seq to determine changes
	private void runLineU () throws IOException, MatrixLoaderException {
		
		LineU lineU = new LineU(originalSeq);
		LineResults result = lineU.getPosResults();
		lineChanges = result.getDifferences();
		String change[] = result.getDifferences().split("\\|");
		
		changes = new ConcurrentHashMap<Integer, Change>();
		
		Pattern basePattern = Pattern.compile("[atcg](\\d+)([atcg])");
		Pattern delPattern = Pattern.compile("[di](\\d+)");
		Pattern delPatternTwo = Pattern.compile("[di](\\d+)\\-(\\d+)");
		Pattern inPattern = Pattern.compile("i(\\d+)([atcg]+)");
						
		for (String d : change) {
			Matcher baseMatcher = basePattern.matcher(d);
			Matcher delMatcher = delPattern.matcher(d);
			Matcher delMatcherTwo = delPatternTwo.matcher(d);
			Matcher inMatcher = inPattern.matcher(d);
			if (baseMatcher.matches()) {
				changes.put(Integer.parseInt(baseMatcher.group(1)) - 1, new Change(baseMatcher.group(2), ChangeType.SNP));
			} else if (delMatcher.matches()) {
				changes.put(Integer.parseInt(delMatcher.group(1)) - 1, new Change(null, ChangeType.DEL));
			} else if (delMatcherTwo.matches()) {
				int start = Integer.parseInt(delMatcherTwo.group(1));
				int stop = Integer.parseInt(delMatcherTwo.group(2));
				for (; start <= stop; start++) {
					changes.put(start - 1, new Change(null, ChangeType.DEL));
				}
			} else if (inMatcher.matches()) {
				changes.put(Integer.parseInt(inMatcher.group(1)) - 1, new Change(inMatcher.group(2), ChangeType.INS));
			}
		}
		
		
	}
	
	//Count Kmers for each of the changes determined by LINEU
	private int countKmer(Map<Integer, Integer> depth) throws ExecuteException, IOException {
		
		executor.execute(CommandLine.parse("/local/projects-t3/eugenetemp/EugeneTools/jellyfish-2.2.3/bin/jellyfish count -m 21 -s 100 -o " + tmpDir.getAbsolutePath() + "/mer_counts.jf " + tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".fq"));
		int totalFail = 0;
		//Check to see if changes are real:
		//Build a string of length 40, and then break it apart and test Kmers
		for (int start = 0; start < meiLoaded.length(); start++) {
//		for (Map.Entry<Integer, Change> change : changes.entrySet()) {
			
			String forward = "";
			String behind = "";
			if (changes.containsKey(start)) {
				
				int testPos = start;
				Change testChange = changes.get(testPos);
				if (!testChange.isVerified()) {
				
					//Go reverse (if possible):
					if (testPos < 50) {
						for (int x = 0; x < testPos; x++) {
							if (changes.containsKey(x)) {
								Change current = changes.get(x);
								if (current.getType() == ChangeType.SNP) {
									behind += current.getChange();
								} else if (current.getType() == ChangeType.INS) {
									behind += current.getChange() + meiLoaded.subsequence(x, 1);
								}
							} else {
								behind += meiLoaded.subsequence(x, 1);
							}
						}
					} else {
						for (int x = (testPos - 50); x < testPos;x++) {
							if (changes.containsKey(x)) {
								Change current = changes.get(x);
								if (current.getType() == ChangeType.SNP) {
									behind += current.getChange();
								} else if (current.getType() == ChangeType.INS) {
									behind += current.getChange() + meiLoaded.subsequence(x, 1);
								}
							} else {
								behind += meiLoaded.subsequence(x, 1);
							}
						}
					}
					
					//Go forward (if possible)
					if (testPos > 5969) {
						if (testChange.getType() == ChangeType.INS) {
							for (int x = testPos; x < 6019; x++) {
								if (changes.containsKey(x)) {
									Change current = changes.get(x);
									if (current.getType() == ChangeType.SNP) {
										forward += current.getChange();
									} else if (current.getType() == ChangeType.INS) {
										forward += current.getChange() + meiLoaded.subsequence(x, 1);
									}
								} else {
									forward += meiLoaded.subsequence(x, 1);
								}
							}
						} else {
							for (int x = testPos+1; x < 6019; x++) {
								if (changes.containsKey(x)) {
									Change current = changes.get(x);
									if (current.getType() == ChangeType.SNP) {
										forward += current.getChange();
									} else if (current.getType() == ChangeType.INS) {
										forward += current.getChange() + meiLoaded.subsequence(x, 1);
									}
								} else {
									forward += meiLoaded.subsequence(x, 1);
								}
							}
						}
					} else {
						if (testChange.getType() == ChangeType.INS) {
							for (int x = testPos; x < (testPos + 50); x++) {
								if (x != testPos && changes.containsKey(x)) {
									Change current = changes.get(x);
									if (current.getType() == ChangeType.SNP) {
										forward += current.getChange();
									} else if (current.getType() == ChangeType.INS) {
										forward += current.getChange();
									}
								} else {
									forward += meiLoaded.subsequence(x, 1);
								}
							}
						} else {
							for (int x = testPos + 1; x < (testPos + 50); x++) {
								if (changes.containsKey(x)) {
									Change current = changes.get(x);
									if (current.getType() == ChangeType.SNP) {
										forward += current.getChange();
									} else if (current.getType() == ChangeType.INS) {
										forward += current.getChange() + meiLoaded.subsequence(x, 1);
									}
								} else {
									forward += meiLoaded.subsequence(x, 1);
								}
							}
						}
						
					}
					
					String toTestBehind = "";
					String toTestForward = "";
					
					if (testChange.getType() == ChangeType.SNP || testChange.getType() == ChangeType.INS) {
						if (behind.length() < 20) {
							toTestBehind += behind;
						} else {
							toTestBehind += behind.substring((behind.length() - 20), behind.length());
						}
						if (forward.length() < 20) {
							toTestForward += forward;
						} else {
							toTestForward += forward.substring(0, 20);
						}
					} else { // DEL
						if (behind.length() < 21) {
							toTestBehind += behind;
						} else {
							toTestBehind += behind.substring((behind.length() - 21), behind.length());
						}
						if (forward.length() < 21) {
							toTestForward += forward;
						} else {
							toTestForward += forward.substring(0, 21);
						}
					}
//					System.out.println("\n----------------------------------------");
//					if (testChange.getType() == ChangeType.DEL) {
//						System.out.println(testPos + "D");
//					} else {
//						System.out.println(testPos + testChange.getChange());
//					}
//					System.out.println(toTestBehind + "-" + toTestForward);
					int minSupport;
					if (testChange.getType() == ChangeType.SNP || testChange.getType() == ChangeType.INS) {
						minSupport = testKmer(toTestBehind + testChange.getChange() + toTestForward);
					} else {
						minSupport = testKmer(toTestBehind + toTestForward);
					}
					
					if (minSupport <= 1) {
//						System.out.println("RETEST");
						int retestMin;
						//Try REF to see if we recover:
						if (testChange.getType() == ChangeType.INS || testChange.getType() == ChangeType.SNP) {
							if (testChange.getType() == ChangeType.INS) {
								retestMin = testKmer(toTestBehind + toTestForward);
							} else {
								retestMin = testKmer(toTestBehind + meiLoaded.subsequence(testPos, 1) + toTestForward);
							}
							if (retestMin > minSupport) {
								changes.remove(testPos);
							} else {
								totalFail++;
							}
						} else { // This is DEL
							//Need to test all possible bases for the del, as a del is often the result of an erroneous SNP call:
							int retestA = testKmer(toTestBehind + "A" + toTestForward);
							int retestG = testKmer(toTestBehind + "G" + toTestForward);
							int retestC = testKmer(toTestBehind + "C" + toTestForward);
							int retestT = testKmer(toTestBehind + "T" + toTestForward);
							if (retestA <= 1 && retestG <= 1 && retestC <= 1 && retestT <= 1) {
								totalFail++;
							} else {
								changes.remove(testPos);
								Change ch = new Change(resolveChange(retestA, retestG, retestC, retestT), ChangeType.SNP);
								ch.setIsVerified(true);
								changes.put(testPos, ch);
							}
						}
					} else {
						testChange.setIsVerified(true);
						changes.put(testPos, testChange);
					}
					
				}
				
			}
		
		}
		
		return totalFail;
				
	}
	//For last iteration, try and correct with reference only for all changes (should help correct clustered errors?).
	private Map<Integer,Change> countKmerFinal(Map<Integer, Integer> depth) throws ExecuteException, IOException {
		
		Map<Integer,Change> finalCorrection = new ConcurrentHashMap<Integer, Change>(changes);
		totalFinalFail = 0;
		//Check to see if changes are real:
		//Build a string of length 40, and then break it apart and test Kmers
		for (int start = 0; start < meiLoaded.length(); start++) {
//		for (Map.Entry<Integer, Change> change : changes.entrySet()) {
			
			String forward = "";
			String behind = "";
			if (finalCorrection.containsKey(start)) {
				
				int testPos = start;
				Change testChange = changes.get(testPos);
				if (!testChange.isVerified()) {
				
					//Go reverse (if possible):
					if (testPos < 50) {
						for (int x = 0; x < testPos; x++) {
							if (changes.containsKey(x)) {
								Change current = changes.get(x);
								if (current.isVerified()) {
									if (current.getType() == ChangeType.SNP) {
										behind += current.getChange();
									} else if (current.getType() == ChangeType.INS) {
										behind += current.getChange() + meiLoaded.subsequence(x, 1);
									}
								} else {
									behind += meiLoaded.subsequence(x, 1);
								}
							} else {
								behind += meiLoaded.subsequence(x, 1);
							}
						}
					} else {
						for (int x = (testPos - 50); x < testPos;x++) {
							if (changes.containsKey(x)) {
								Change current = changes.get(x);
								if (current.isVerified()) {
									if (current.getType() == ChangeType.SNP) {
										behind += current.getChange();
									} else if (current.getType() == ChangeType.INS) {
										behind += current.getChange() + meiLoaded.subsequence(x, 1);
									}
								} else {
									behind += meiLoaded.subsequence(x, 1);
								}
							} else {
								behind += meiLoaded.subsequence(x, 1);
							}
						}
					}
					
					//Go forward (if possible)
					if (testPos > 5969) {
						if (testChange.getType() == ChangeType.INS) {
							for (int x = testPos; x < 6019; x++) {
								if (changes.containsKey(x)) {
									Change current = changes.get(x);
									if (current.isVerified()) {
										if (current.getType() == ChangeType.SNP) {
											forward += current.getChange();
										} else if (current.getType() == ChangeType.INS) {
											forward += current.getChange() + meiLoaded.subsequence(x, 1);
										}
									} else {
										forward += meiLoaded.subsequence(x, 1);
									}
								} else {
									forward += meiLoaded.subsequence(x, 1);
								}
							}
						} else {
							for (int x = testPos + 1; x < 6019; x++) {
								if (changes.containsKey(x)) {
									Change current = changes.get(x);
									if (current.isVerified()) {
										if (current.getType() == ChangeType.SNP) {
											forward += current.getChange();
										} else if (current.getType() == ChangeType.INS) {
											forward += current.getChange() + meiLoaded.subsequence(x, 1);
										}
									} else {
										forward += meiLoaded.subsequence(x, 1);
									}
								} else {
									forward += meiLoaded.subsequence(x, 1);
								}
							}
						}
					} else {
						if (testChange.getType() == ChangeType.INS) {
							for (int x = testPos; x < (testPos + 50); x++) {
								if (x != testPos && changes.containsKey(x)) {
									Change current = changes.get(x);
									if (current.isVerified()) {
										if (current.getType() == ChangeType.SNP) {
											forward += current.getChange();
										} else if (current.getType() == ChangeType.INS) {
											forward += current.getChange();
										}
									} else {
										forward += meiLoaded.subsequence(x, 1);
									}
								} else {
									forward += meiLoaded.subsequence(x, 1);
								}
							}
						} else {
							for (int x = testPos + 1; x < (testPos + 50); x++) {
								if (changes.containsKey(x)) {
									Change current = changes.get(x);
									if (current.isVerified()) {
										if (current.getType() == ChangeType.SNP) {
											forward += current.getChange();
										} else if (current.getType() == ChangeType.INS) {
											forward += current.getChange() + meiLoaded.subsequence(x, 1);
										}
									} else {
										forward += meiLoaded.subsequence(x, 1);
									}
								} else {
									forward += meiLoaded.subsequence(x, 1);
								}
							}
						}
						
					}
					
					String toTestBehind = "";
					String toTestForward = "";
					
					if (behind.length() < 20) {
						toTestBehind += behind;
					} else {
						toTestBehind += behind.substring((behind.length() - 20), behind.length());
					}
					if (forward.length() < 20) {
						toTestForward += forward;
					} else {
						toTestForward += forward.substring(0, 20);
					}
					
//					System.out.println("\n----------------------------------------");
//					if (testChange.getType() == ChangeType.DEL) {
//						System.out.println(testPos + "D");
//					} else {
//						System.out.println(testPos + testChange.getChange());
//					}
//					System.out.println(toTestBehind + "-" + toTestForward);
//					System.out.println(toTestBehind + meiLoaded.subsequence(testPos, 1) + toTestForward);
					
					int minSupport;
					minSupport = testKmer(toTestBehind + meiLoaded.subsequence(testPos, 1) + toTestForward);
					if (minSupport <= 1) {
						totalFinalFail++;
					} else {
						finalCorrection.remove(testPos);
					}
					
				}
				
			}
		
		}
		
		if (totalFinalFail == 0) {
			return finalCorrection;
		} else {
			return null;
		}
				
	}
	//Utility commands for kmer counting;
	private int testKmer (String toTest) throws ExecuteException, IOException {
		int minSupport = Integer.MAX_VALUE;
		for (int x = 0; x <= (toTest.length() - 21); x++) {
			String currentSubstr = toTest.substring(x, (x + 21));
			String currentCommand = "/local/projects-t3/eugenetemp/EugeneTools/jellyfish-2.2.3/bin/jellyfish query " + tmpDir.getAbsolutePath() + "/mer_counts.jf " + currentSubstr;
			int totalSupport = runCommand(currentCommand);
			if (totalSupport < minSupport) {
				minSupport = totalSupport;
			}
//			System.out.println(currentSubstr + "\t" + (x + 1) + "\t" + (x + 22) + "\t" + totalSupport);
		}
		return minSupport;
	}
	private int runCommand (String command) throws ExecuteException, IOException {
		
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PumpStreamHandler ps = new PumpStreamHandler(os);
		executor.setStreamHandler(ps);
		executor.execute(CommandLine.parse(command));
		executor.setStreamHandler(ps);
		BufferedReader r = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(os.toByteArray())));
		String line;
		int returnable = -1;
		while ((line = r.readLine()) != null) {
			returnable = Integer.parseInt(line.split(" ")[1]);
		}
		r.close();
		return returnable;
		
	}
	private String resolveChange(int A, int G, int C, int T) {
		
		String trueChange;
		
		if (A > G && A > C && A > T) {
			trueChange = "A";
		} else if (G > A && G > C && G > T) {
			trueChange = "G";
		} else if (C > G && C > A && C > T) {
			trueChange = "C";
		} else {
			trueChange = "T";
		}
		
		return trueChange;
	
	}
	
	//Build finalSequence to return
	private void buildFinal() throws IOException, MatrixLoaderException {
		
		finalSeq = "";
		for (int start = 0; start < meiLoaded.length(); start++) {
			if (changes.containsKey(start)) {
				Change current = changes.get(start);
				if (current.isVerified()) {
					if (current.getType() == ChangeType.SNP) {
						finalSeq += current.getChange();
					} else if (current.getType() == ChangeType.INS) {
						finalSeq += current.getChange() + meiLoaded.subsequence(start, 1);
					}
				} else {
					if (current.getType() == ChangeType.SNP || current.getType() == ChangeType.DEL) {
						finalSeq += "N";
					} else if (current.getType() == ChangeType.INS) {
						int len = current.getChange().length();
						String add = "";
						for (int x = 1; x <= len; x++) {
							add += "N" + meiLoaded.subsequence(start, 1);
						}
						finalSeq += add;
					}
				}
			} else {
				finalSeq += meiLoaded.subsequence(start, 1);
			}
		}
		BufferedWriter fastaWriter = new BufferedWriter(new FileWriter(new File(tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".cons.fa"),false));
		fastaWriter.write(">original\n");
		fastaWriter.write(originalSeq);
		fastaWriter.newLine();
		fastaWriter.write(">corrected\n");
		fastaWriter.write(finalSeq);
		fastaWriter.newLine();
		fastaWriter.close();
		
		LineU lineU = new LineU(finalSeq);
		lineNewChanges = lineU.getPosResults().getDifferences();
		
		//Align to L1 REF
		executor.execute(CommandLine.parse("bowtie2 -f -x /local/aberdeen2rw/DEVINE/SVU41//ILLUMINA/MERGE_ANALYSIS/me_refs/LINE1 -U " + tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".cons.fa -S " + tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".cons.sam"));
		executor.execute(CommandLine.parse("samtools view -Sbo " + tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".cons.bam " + tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".cons.sam"));
		executor.execute(CommandLine.parse("samtools sort " + tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".cons.bam " + tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".cons.sorted"));
		executor.execute(CommandLine.parse("samtools index " + tmpDir.getAbsolutePath() + "/" + chr + "_" + start + ".cons.sorted.bam"));
		
	}
	
	//Classes for holding information on read-pairs and LINE-1 changes

	private class Pairs {
		
		private SAMRecord firstPair;
		private SAMRecord secondPair;
		private boolean bothPairs;
		
		private Pairs (SAMRecord read) {
			parseRead(read);
		}
		
		private void parseRead (SAMRecord read) {
			
			if (read.getReadPairedFlag()) {
				boolean isFirst = read.getFirstOfPairFlag();
				if (isFirst && firstPair == null) {
					setFirstPair(read);
				} else if (!isFirst && secondPair == null) {
					setSecondPair(read);
				}
			}
			
		}
		private void setFirstPair (SAMRecord read) {
			firstPair = read;
			if (firstPair != null && secondPair != null) {
				bothPairs = true;
			}
		}
		private void setSecondPair (SAMRecord read) {
			secondPair = read;
			if (firstPair != null && secondPair != null) {
				bothPairs = true;
			}
		}
		
		private boolean isBothPairs() {
			return bothPairs;
		}
		private SAMRecord getFirstPair() {
			return firstPair;
		}
		private SAMRecord getSecondPair() {
			return secondPair;
		}
		
	}
	private class Change {
		
		private String change;
		private ChangeType type;
		private boolean verified;
				
		private Change (String change, ChangeType type) {
			this.change = change;
			this.type = type;
			verified = false;
		}

		private String getChange() {
			if (change == null) {
				return null;
			} else {
				return change.toUpperCase();
			}
		}

		private ChangeType getType() {
			return type;
		}
		
		private void setIsVerified(boolean isVerified) {
			verified = true;
		}
		
		private boolean isVerified() {
			return verified;
		}
		
	}
	
	private enum ChangeType {
		SNP,INS,DEL;
	}
	
}
