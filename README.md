# PESC
PESC is a small python script that can softclip the primers from paired sequencing data derived primer extension targeted sequencing.
It takes SAM input from stdin and a primer file and sends the modified SAM records to stdout for further processing or storage.
The tool requires python 2.7

example: samtools view -h my_test.bam | python primer_removal9_test.py CDHS-32375Z-340.primer3.txt | 
java -Xmx8G -jar picard.jar SortSam INPUT=/dev/stdin OUTPUT=clipped_test.bam SORT_ORDER=coordinate; samtools index clipped_output.bam

example of a primer file (please do not include the header!):


chromosome  position  strand  primer-sequence

chr1	115252238	0	GAGGTATCAATGTATGGAATCCCGTGCATCTTG

The position should match the last basepair of the primer sequence at the border to the actual start of the captured DNA fragment. It is important to have the complete primer sequence, as the program determines the length to be clipped based on the length of the primer. If only the length of the primer but not its sequence is known, it is instead possible to just pass a dummy sequence of the correct length.
