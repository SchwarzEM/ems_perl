#!/gsc/bin/perl

use strict;

my $in = $ARGV[0];

open(IN, "<$in");

my %gene_hash;
my %cds_hash;

# get hash

foreach(<IN>) {
    next if ($_ =~ /^\#/);
    my @f = split(/\t/, $_);

    if($f[2] eq 'gene') {
	my ($gene) = $f[8] =~ /Name\=(\S+)\;/;
	my ($alias) = $f[8] =~ /Alias\=(\S+)/;
	$gene_hash{$gene}{"alias"} = $alias;
    }
    elsif($f[2] eq 'mRNA') {
	my ($gene) = $f[8] =~ /ID\=mRNA\:(\S+)/;
	my ($parent) = $f[8] =~ /Parent\=Gene\:(\S+)\;/;
	my $contig = $f[0]; 
	$contig =~ s/^crem3_c/C/;
	$cds_hash{$gene}{"contig"} = $contig;
	$cds_hash{$gene}{"parent"} = $parent;
	$cds_hash{$gene}{"strand"} = $f[6];
	$cds_hash{$gene}{"start"} = $f[3];
	$cds_hash{$gene}{"stop"} = $f[4];
    }
    elsif($f[2] eq 'CDS') {
	my ($gene) = $f[8] =~ /Parent\=mRNA\:(\S+)/;
#	print "CDS gene = $gene\n";
	my @array = ($f[3], $f[4]);
	if($cds_hash{$gene}{internal}[0]) {
	    my @array1 = @{$cds_hash{$gene}{internal}};
	    push(@array1, @array);
	    $cds_hash{$gene}{"internal"} = \@array1;
	}
	else {
	    $cds_hash{$gene}{"internal"} = \@array;
#	    print "@{$cds_hash{$gene}{internal}}\n";
	}

    }



}



# print ace format

foreach(sort keys %cds_hash) {

    my $gene = $_;

    print "// printing $gene in ace format\n";
    print "// $gene, $cds_hash{$gene}{contig}, $gene_hash{$cds_hash{$gene}{parent}}{alias}, $cds_hash{$gene}{strand}, @{$cds_hash{$gene}{internal}}\n";

    &print_as_ace($gene, $cds_hash{$gene}{contig}, $gene_hash{$cds_hash{$gene}{parent}}{alias}, $cds_hash{$gene}{strand}, $cds_hash{$gene}{internal});


}




sub print_as_ace {

    my ($gene, $contig, $alias, $strand, $coords) = @_;
    my @array = @$coords;
    my @scoords = sort {$a <=> $b} (@array);

    print "// gene : $gene\n";
    print "// @scoords\n\n";

    if($strand eq '+') {
	print "Sequence : \"$contig\"\n";
	print "CDS_child $gene $scoords[0] $scoords[$#scoords]\n\n";

	print "CDS : \"$gene\"\n";
	for(my $ii = 0; $ii <= $#scoords; $ii+=2) {
	    my $c1 = $scoords[$ii] - $scoords[0] + 1;
	    my $c2 = $scoords[$ii+1] - $scoords[0] + 1;
	    print "Source_exons $c1 $c2\n";
	}
	print "CDS\n";
	print "Method \"washu_remanei\"\n";
	print "Remark \"WBGene id : $alias\"\n";
	print "\n";

    }

    else {

	print "Sequence : \"$contig\"\n";
	print "CDS_child $gene $scoords[$#scoords] $scoords[0]\n\n";

	print "CDS : \"$gene\"\n";
	for(my $ii = $#scoords; $ii >= 0; $ii-=2) {
	    my $c1 = $scoords[$#scoords] - $scoords[$ii] + 1;
	    my $c2 = $scoords[$#scoords] - $scoords[$ii-1] + 1;
	    print "Source_exons $c1 $c2\n";
	}
	print "CDS\n";
	print "Method \"washu_remanei\"\n";
	print "Remark \"WBGene id : $alias\"\n";
	print "\n";



    }


}
