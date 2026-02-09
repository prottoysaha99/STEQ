#!/lusr/bin/perl -w

#Author: Md. Shamsuzzoha Bayzid
#September 04, 2014

#use strict;
use warnings;
use Getopt::Long;

#output
sub badInput {
  my $message = "Usage: perl $0 takes the subset file obtained from DACTAL and creates seperate files containing one subset each. It also creates the species list file required to run MP-EST.
	-i=<subset_file>";
  print STDERR $message;
  die "\n";
}

GetOptions(
	"i=s"=>\my $subsets,
	"o=s"=>\my $output,
);

badInput() if not defined $subsets;


	my $subset_file_out = $output . "tmp.subsets";
	open(INFO, $subsets);		# Open the file
	my @lines = <INFO>;		# Read it into an array
	close(INFO);

	open(OUT, ">", $subset_file_out) or die "can't open $subset_file_out: $!";


	foreach my $line (@lines){
	$line =~ s/([^\[]*)//; # removing everything before the first [
	$line =~ s/[',\[\]]//g; # removing ' , [ and ]
	my $test = $line;   

		if($test){
		print OUT "$line";
		}
      }
	close(OUT);


open(INFO, $subset_file_out);		# Open the file
@lines = <INFO>;		# Read it into an array
close(INFO);


my $i = 1;
foreach my $line (@lines)
	{
		chomp($line); # to remove trailing new lines from the $line
		my @taxa = split(/ /, $line);   # so make sure there is no space after the last taxa name
		my $n_taxa = scalar (@taxa);
	
		

	$n_taxa = scalar (@taxa);  # size of the @taxa
	#print "\n no of taxa after padding: $n_taxa";
	
	my $subset_r = $output . "subsets.$i";  # for mpest.
	open(OUT1, ">", $subset_r) or die "can't open $subset_r: $!";

	
	
	
		foreach $taxa (@taxa)
		{
		print OUT1 "$taxa\n";
		}

		
		$i++; # counter for subsets
}
		  
