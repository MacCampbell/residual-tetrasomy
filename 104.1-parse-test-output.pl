#! /usr/bin/perl -w
# Usage
# ./104.1-parse-test-output.pl ./outputs/104/salmo-salar-eight-compare.txt
# for f in ./outputs/104/*eight-compare.txt; do echo $f; ./104.1-parse-test-output.pl $f; done;

my @pairs=();
my @results=();
while(<>){
  chomp;
  if ($_=~/^\[1\]/) {
    my @line=split(" ", $_);
    push(@pairs, "$line[1]\t$line[2]")
  } 
  elsif ($_=~/^W = /) {
    my @nextLine=split(" ", $_);
    $nextLine[2]=~ s/\,//g;
    push (@results, "$nextLine[2]\t$nextLine[5]")
  }
}

print "Chrom1\tChrom2\tW\tp-value\n";
for (my $i=0; $i< scalar(@pairs); $i++) {
 print $pairs[$i]."\t".$results[$i]."\n";

}
