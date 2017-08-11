package tenX;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.IOException;

public class tenXparser {

	public static void main(String[] args) throws IOException {

		//SET WHEN RUNNING: -Djava.io.tmpdir=/some/tmp/dir/
		
		SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		SamReader tenXreader = readerFactory.open(SamInputResource.of(args[0]).index(new File(args[0] + ".bai")));
		sortTenX.makeTagBam(tenXreader, new File(args[0] + ".tagSorted.sam"));
		
	}
		
}
