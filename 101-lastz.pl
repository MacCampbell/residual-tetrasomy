#! /usr/bin/perl -w
# By Mac Campbell DrMacCampbell@gmail.com
# The goal is to drive multiple lastz comparisons using
# a "myPairs.txt" file that is space-delimited and
# gnu-parallel

# Usage: ./101-lastz.pl processors genus-species path/to/genus-species-pairs
# Example: ./101-lastz.pl 2 salmo-salar ./data/salmo-salar/salmo-salar-pairs.txt
# How many simultaneous runs to do
my $processors = shift;

my $species = shift;

# Make sure output dir exists
`mkdir outputs/101/`;

my @pairs;

while(<>) {
	chomp;
	push (@pairs, $_);

}

print "The pairs I am going to compare are:\n";
print join ("\n", @pairs);
print "\n";

#Prepare file for parallel
open(OUTFILE, ">temp.txt") or die "Cannot open outfile";
	foreach my $pair (@pairs) {
	my @a=split(/\s/, $pair);
	# Example with subsetting chromos
	# printf OUTFILE "lastz ./data/$species/$species-$a[0].fasta\[1..2000000\] ./data/$species/$species-$a[1].fasta\[1..2000000\] --chain --gapped --gfextend --identity=75.0..100.0 --matchcount=100 --format=general --rdotplot=./outputs/101/$species-$a[0]-$a[1]-dotplot.txt --ambiguous=iupac --exact=20 > ./outputs/101/$species-$a[0]-$a[1].out\n";
  # Full pairs
   printf OUTFILE "lastz ./data/$species/$species-$a[0].fasta ./data/$species/$species-$a[1].fasta --chain --gapped --gfextend --identity=75.0..100.0 --matchcount=100 --format=general --rdotplot=./outputs/101/$species-$a[0]-$a[1]-dotplot.txt --ambiguous=iupac --exact=20 > ./outputs/101/$species-$a[0]-$a[1].out\n";

	}
    
close(OUTFILE);

`parallel -j $processors < temp.txt`;

# Combine outfiles while excluding lines beginning with #
`cat ./outputs/101/$species*.out | grep -v "#" > ./outputs/101/$species-101-lastz.result`;

# Clean up 
#`"yes" | rm *.out`;
`"yes" | rm temp.txt`;
