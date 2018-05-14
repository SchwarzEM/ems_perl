#!/usr/bin/perl -w
# 
# Select selection_number rows randomly from a list
#
# Usage: random_select.pl listfile selection_number
#
# Ross Hall 2010
#

use strict;

if (@ARGV != 2) {
	print STDERR "Usage: random_select.pl listfile selection_number\n";
	exit(1);
}

my $infile = shift;
my $snumber = shift;


open(FFILE,$infile) || &ErrorMessage("Cannot open file ".$infile);
my @larray = <FFILE>;


my $seqnum = @larray;


print STDERR "$seqnum rows in the file: $infile\n";	
my @randarray;
my %randhash;

my $range = $seqnum-1;
my $random_number;

my $count = 0;

while ($count < $snumber) {
	$random_number = int(rand($range));
	while (exists($randhash{$random_number})) {
		$random_number = int(rand($range));
	}
	$randhash{$random_number} = 1;
	$randarray[$count] = $random_number;
	$count++;
}
	
print STDERR "Created $snumber random list numbers\n";

			
my @sortednums = sort { $a <=> $b } @randarray; # numerical sort 


my $seqcount = 0;
my $randcurrent = 0;
my $printflag = 0;
my $maxrandindex = $range-2;

print STDERR "Printing the random samples sequences.........\n";

foreach  my $x (@sortednums) {
	print $larray[$x];
}		
			
			
 	
							
sub ErrorMessage {
	my $msg = shift;
	print "Fatal error: $msg\n";
	exit(1);	
}
		

