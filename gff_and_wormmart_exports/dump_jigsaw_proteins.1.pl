#!/gsc/bin/perl

use strict;
use Ace;
use Bio::SeqIO;
use Bio::Seq;
use lib "$ENV{WORMBASE_SCRIPT_DISK}/lib/WormBase";
use Database;

my $species = $ARGV[0];


my %subsequencehash;
my %parenthash;
my %dnahash;

my $db;
my @sequences;

if($species eq 'briggsae') {
    $db = Database::connect_to_database('brigace');
    @sequences = $db->fetch(-query => "Find Sequence chr*");
}
elsif($species eq 'brenneri') {
    $db = Database::connect_to_database('brenace');
    @sequences = $db->fetch(-query => "Find Sequence Contig*");
}
elsif($species eq 'remanei') {
    $db = Database::connect_to_database('remanei_acedb');
    @sequences = $db->fetch(-query => "Find Sequence Contig*");
}
elsif($species eq 'brugia') {
    $db = Database::connect_to_database('brugia_acedb');
    @sequences = $db->fetch(-query => "Find Sequence Bmal*");
}


foreach (@sequences) {
    my $clone = $_;
    my $obj = $db->fetch(Sequence => $clone);
    my @subsequences = $obj->CDS_child(1);
	
    my $dna = $obj->asDNA;
    $dna =~ s/^\>\S+//;
    $dna =~ s/\s+//g;
    $dnahash{$clone} = $dna;

# add subsequences to %subsequencehash         #
    foreach(@subsequences) {
	if($species eq 'brugia') {
	    next unless($_ =~ /^Bmal/);
	}
	else {
	    next unless ($_ =~ /^Contig/); 
	}
	my ($a, $b, $c) = $_->row;
	my @ends = ($b, $c);
	push(@{$subsequencehash{$a}}, @ends);
	$parenthash{$a} = $clone;
    }
}


my @cds;

if($species eq 'brugia') {
    @cds = $db->fetch(-query=>"Find CDS Bmal*");

}
else {
    @cds = $db->fetch(-query=>"Find CDS Contig*");
}


foreach(keys %subsequencehash) {
    my $gene = $_;
    my $obj = $db->fetch(CDS => $gene);
    my @exons = $obj->Source_exons(1);
    
    foreach(@exons) {
	my ($a, $b) = $_->row;
	my @ends = ($a, $b);
	push(@{$subsequencehash{$gene}}, @ends);
    }
}



my $currentparent = "";
my $clonedna;
my $revcompclonedna;
    
foreach (sort keys %subsequencehash) {
    my $subclone = $_;
#    my ($parentclone) = $subclone =~ /^(\w+)\./;
    my $parentclone = $parenthash{$subclone};
#	print "checking $subclone (currentparent = $currentparent parentclone = $parentclone)\n";
    
#### check for Start_not_found and End_not_found tags
    my $obj = $db->fetch(CDS => $subclone);
    
# get the dna we want to work with for this next set of genes (ie, until the child changes -- example: B0205.1, (use B0205s dna) then .2, then .3, then .4, then B0280.1 --> change to B0280s dna here)  
    unless($currentparent eq $parentclone) {   
	$currentparent = $parentclone;         
	foreach(keys %dnahash) {
	    if(uc($_) eq uc($parentclone)) {
		$clonedna = $dnahash{$_}; # access the sequence
		my $revcompcloneobj = Bio::Seq->new('-seq'=>"$clonedna");
		my $revcompclone = $revcompcloneobj->revcom();
		$revcompclonedna = $revcompclone->seq();
	    }
	}
    }
    my $totalseq = 0;
    my $totalsubseq = 0;
	
# stitch together spliced dna sequence using substr()        #
    
    my @total = @{$subsequencehash{$subclone}};
#    print "total $#total\n";
    my $start = $subsequencehash{$subclone}->[0];
    my $end = $subsequencehash{$subclone}->[1];
    
    my $subclonedna = "";
    
    for(my $count=2; $count<=$#total; $count+=2) {
	if($start < $end) {
	    my $exonstart = $start + $subsequencehash{$subclone}->[$count] - 2;
	    my $exonend = $start + $subsequencehash{$subclone}->[$count+1] - 2;
	    my $length = $exonend - $exonstart + 1;
	    my $thisseg = substr($clonedna, $exonstart, $length);
	    $subclonedna .= $thisseg;
#	    print "$exonstart $exonend\n";
#	    print "substr: $thisseg\n";
	}
	else {
	    my $exonstart = $start - $subsequencehash{$subclone}->[$count] + 1;
	    my $exonend = $start - $subsequencehash{$subclone}->[$count+1] + 1;
	    my $newexonstart = length($clonedna) - $exonstart;
	    my $newexonend = length($clonedna) - $exonend;
	    my $length = $newexonend - $newexonstart + 1;
	    my $thisseg = substr($revcompclonedna, $newexonstart, $length);
	    $subclonedna .= $thisseg;
#	    print "$exonstart $exonend ";
#	    print "$thisseg\n";
	}
    }
#	print "$subclonedna\n";

# use Bioperl module Seq.pm to translate sequence       #
    unless($subclonedna =~ /\S/) {
	print "couldn't get a sequence for $subclone\n";
#	    die;
	
    }
    
#	if($subclone =~ /W06H8/) {
#	    print "trying to translate $subclone ($subclonedna)\n";
#	}
    my $thisseqobj = Bio::Seq->new('-seq'=>"$subclonedna");
    my $tranthisseq = $thisseqobj->translate();
    my $translation = $tranthisseq->seq();

    $translation =~ s/\*$//;

# write as fasta format

    print ">$subclone\n";
    print &as_fasta($translation);


}


    
sub as_fasta {
    my ($str) = @_;
    my $fastastr = "";
    my $length = length($str);
    my $linelength = 60;
    for(my $start = 0; $start <= $length; $start+=$linelength) {
        my $this_seg = substr($str, $start, $linelength);
        $fastastr .= "$this_seg\n";
    }
    return $fastastr;
}







