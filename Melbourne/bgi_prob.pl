#!/usr/bin/perl -w

# Script to read in a list of the 5 params required
# for the BGI differential gene expression prob. calculation
#
# Datafile format: name length x y N1 N2
# separated by tabs
#
# Output: 
# Contig Length x y N1 N2 BGIpsum p(y|x) RPKM_x RPKM_y RPKM_log2ratio Log2x/Log2y Log2y/Log2x FDR
# 
# Usage: bgiprob.pl replacement_value datafile
# replacement_value is the number to replace zero if x=0 or y=0
#
# Ross Hall 2010


use lib '/home/est/bin';
use lib '/home/ross/Lib/Perl';
use strict;
use ProbBGI;



if (@ARGV != 3) {
	print STDERR "Usage: bgi_prob.pl [P BGIPsum] replacement_value datafile\n";
	print STDERR "FDR calculation option, P: use single p value or BGIPsum: use the summed p value as per BGI\n";
	print STDERR "replacement_value is the number to replace zero if x=0 or y=0\n";
	print STDERR "Example: bgi_prob.pl P 0.0001 mydatafile\n";
	exit(1);
}

my $popt = shift;
my $replace = shift;			
my $infile = shift;


my $PI = 3.1415926535897932384626433832795;

open(INFILE,$infile) || &ErrorMessage("Cannot open file ".$infile);
my @darray = <INFILE>;

if ($popt eq "BGIPsum") {
	print join("\t","#Contig","Length","x","y","N1","N2","BGIpsum","RPKM_x","RPKM_y","RPKM_log2ratio","Log2x/Log2y","Log2y/Log2x","FDR") . "\n";
} elsif ($popt eq "P") {
	print join("\t","#Contig","Length","x","y","N1","N2","p(y|x)","RPKM_x","RPKM_y","RPKM_log2ratio","Log2x/Log2y","Log2y/Log2x","FDR") . "\n";
} else {
	&ErrorMessage("Unkown P option: " . $popt);
}			


my $bgiobj;
my @data_array;
foreach my $line (@darray) {
	chomp($line);
	next if ($line =~ /^#/);
	 
	my @fa = split(/\t/,$line);
	my $name = $fa[0];
	my $length = $fa[1];
	my $x = $fa[2];
	my $y = $fa[3];
	my $n1 = $fa[4];
	my $n2 = $fa[5];
	
	if ($x < 0.000001) { $x = $replace;}
	if ($y < 0.000001) { $y = $replace;}
	
	# print STDERR join("\t",$x,$y,$n1,$n2) . "\n";
	
	my @parray;
	
	my $singlep = &bgi_prob($x,$y,$n1,$n2);
	my $singlestr = sprintf("%3.5e",$singlep);
	
	my $psum = 0;
	for (my $i=0;$i<=$y;$i++) {
		my $pval = &bgi_prob($x,$i,$n1,$n2);
		$psum += $pval;
		$parray[$i] = $pval;
	}
	if ($psum > 0.5) {
		$psum = 1 - $psum;
	}
	$psum *= 2;	
	
	
	if ($x == 1) { $x = 1.0000001;}
	if ($y == 1) { $y = 1.0000001;}
	
	my $rpkm_x = (1000000*$x)/($n1*$length/1000);
	my $rpkm_y = (1000000*$y)/($n2*$length/1000);
	
	my $rpkm_x_str = sprintf("%3.5f",$rpkm_x);
	my $rpkm_y_str = sprintf("%3.5f",$rpkm_y);
	# print STDERR "$x\tlog2(x): " . &log2($x) . "\n";
	my $log2ratio = sprintf("%3.5f",&log2($rpkm_x/$rpkm_y));
	my $log2x_y = sprintf("%3.5f",&log2($x)/&log2($y));
	my $log2y_x = sprintf("%3.5f",&log2($y)/&log2($x));
					
	my $numstr = sprintf("%3.5e",$psum);
	
	
	$bgiobj = ProbBGI->new(contig => $name, contiglength => $length, x => $x, y => $y, n1 => $n1,
				n2 => $n2, bgisum =>$psum, bgipval => $singlep, rpkm_x => $rpkm_x, rpkm_y => $rpkm_y,
				rpkm_ratio => &log2($rpkm_x/$rpkm_y), log2xy => &log2($x)/&log2($y),
				log2yx => &log2($y)/&log2($x), fdr => -1);
							
	push(@data_array,$bgiobj);
		
}

my @sorted_array;
if ($popt eq "BGIPsum") {
	@sorted_array = sort {$a->bgisum <=> $b->bgisum} @data_array;
} 
elsif ($popt eq "P") {
	@sorted_array = sort {$a->bgipval <=> $b->bgipval} @data_array; 
} 

my $asize = scalar(@sorted_array);
for (my $i=1;$i<=$asize;$i++) {
	my $bobj = $sorted_array[$i-1];

	if ($popt eq "BGIPsum") {
		$bobj->setFDR($i*$sorted_array[$i-1]->bgisum/$asize);
		$bobj->printBGI();
	} else {
		$bobj->setFDR($i*$sorted_array[$i-1]->bgipval/$asize);
		$bobj->printP();
	}		
}

exit(0);	



sub log2 {
	my $x = shift;
	return(log($x)/log(2));
}
	

sub log10 {

  	my $n = shift;
  	return log($n)/log(10);
}



sub bgi_prob {
	my $x = shift;
	my $y = shift;
	my $n1 = shift;
	my $n2 = shift;	
	
	my $p = $y*&log10($n2/$n1) + &lnfact($x+$y) - &lnfact($x) - &lnfact($y) - ($x+$y+1)*&log10(1+($n2/$n1));
	
	return(10**$p);
	
}
	
	
	
sub lnfact {
	my $x = shift;
	
	if ($x == 0) { return(0); }
	
	my $lnf1 = $x*log($x)  - $x + (1/2)*log(2*$PI*$x) + 1/(12*$x);
	my $inf2 = 1/(360*exp(3*log($x)));
	my $inf3 = 1/(1260*exp(5*log($x)));
	my $lnf = $lnf1 - $inf2 + $inf3; 
	
	return($lnf/log(10));
}


		
sub ErrorMessage {
	my $msg = shift;
	print STDERR "Fatal error: $msg\n";
	exit(1);	
}	
