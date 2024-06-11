#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $gff3 = q{};
my $pep  = q{};

$gff3 = $ARGV[0] if $ARGV[0];
$pep  = $ARGV[1] if $ARGV[1];

if ( (! $gff3 ) or (! $pep ) ) {
    die "Format: make_gff3.lift.filt_11jan2024.pl [GFF3] [proteome] > [GFF3 filtering script]\n";
}

my $pre_filt_gff3 = $gff3;
my $filt_gff3     = $gff3;

$pre_filt_gff3 =~ s/\.gff3\z/.pre-filt.gff3/;
$filt_gff3     =~ s/\.gff3\z/.filt.gff3/;

$pre_filt_gff3 = safename($pre_filt_gff3); 
$filt_gff3     = safename($filt_gff3);

my $genes_of_gff3     = 'genes_of_gff3.txt';
my $genes_of_proteome = 'genes_of_proteome.txt';

$genes_of_gff3     = safename($genes_of_gff3); 
$genes_of_proteome = safename($genes_of_proteome); 

my $mRNAs_of_gff3      = 'mRNAs_of_gff3.txt';
my $mRNAs_of_proteome  = 'mRNAs_of_proteome.txt';
my $mRNAs_of_gff3_only = 'mRNAs_of_gff3_only.txt';

$mRNAs_of_gff3      = safename($mRNAs_of_gff3);
$mRNAs_of_proteome  = safename($mRNAs_of_proteome);
$mRNAs_of_gff3_only = safename($mRNAs_of_gff3_only);

print "\n";

print "$gff3";
print ' | perl -ne \' if ( /\A \S+ \t \S+ \t gene \t .* \t ID=(\S+?);/xms ) { print "$1\n" ; }\' | sort | uniq >';
print " $genes_of_gff3 ;\n";

print "\n";

print "$pep";
print ' | grep \'>\' | perl -ne \' if ( /gene_id=(\S+?);/xms ) { print "$1\n"; } \' | sort | uniq >';
print " $genes_of_proteome ;\n";

print "\n";

print "$gff3";
print ' | perl -ne \' if ( /\A \S+ \t \S+ \t mRNA \t .* \t ID=(\S+?);/xms ) { print "$1\n" ; }\' | sort | uniq >';
print " $mRNAs_of_gff3 ;\n";

print "\n";

print "$pep";
print ' | grep \'>\' | perl -ne \' if ( /\A[>](\S+)/xms ) { print "$1\n"; } \' | sort | uniq >';
print " $mRNAs_of_proteome ;\n";

print "\n";

print "comm -23 $mRNAs_of_gff3 $mRNAs_of_proteome > $mRNAs_of_gff3_only ;\n";

print "\n";

print '$PROJECT/ems_perl/Hco_Acey/select_gff3_genes_28may2024.pl';
print " $genes_of_proteome $gff3 1>$pre_filt_gff3 2>test1.err ;\n";

print "\n";

print '$PROJECT/ems_perl/Hco_Acey/reject_gff3_txs_28may2024.pl'; 
print " $mRNAs_of_gff3_only $pre_filt_gff3 1>$filt_gff3 2>test2.err ;\n";

print "\n";

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

