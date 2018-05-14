#!/usr/bin/perl -w
#
# Script to join two file togethers on common column one.
# Output has all of file 1 and only matching of file two.
#
# Note: if no matching from file 2, place a "#" character in the 
# empty field locations
#
# 
# Usage: superjoin.pl field_separator file1 file2
#
# Ross Hall 2010

use strict;

if (@ARGV != 3) {
	print STDERR "superjoin.pl field_separator file1 file2\n";
	exit(1);
}


my $sep = shift;
my $file1 = shift;
my $file2 = shift;


open(FILEONE,$file1) || &ErrorMessage("Cannot open file ". $file1);
my @f1array = <FILEONE>;


open(FILETWO,$file2) || &ErrorMessage("Cannot open file ".$file2);
my @f2array = <FILETWO>;


my %f1hash;
my %f2hash;

foreach my $line (@f2array) {
	chomp($line);
	
	my @fa = split(/$sep/,$line);
	
	my $key = $fa[0];
	if (exists($f2hash{$key})) {
		my $arrayptr = $f2hash{$key};
		push(@$arrayptr,$line);
	} else {
		my @array =	();
		push(@array,$line);
		$f2hash{$key} = \@array;
	}			
}


foreach my $line (@f1array) {
	chomp($line);
	
	my @fa = split(/$sep/,$line);
	
	my $key = $fa[0];
	if (exists($f2hash{$key})) {
		my $aptr = $f2hash{$key};
		foreach my $x (@$aptr) {
			print $line . $sep . $x . "\n";
		}
	} else {
		print $line . $sep . "#" . $sep . "#\n";
 
	}	
}			
	
	





sub ErrorMessage {
	my $msg = shift;
	print STDERR "Fatal error: $msg\n";
	exit(1);	
}	

