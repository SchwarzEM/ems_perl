#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $infile = q{};
my $header = q{};

my $data_ref;

$infile = $ARGV[0] if $ARGV[0];

if (! -r $infile ) {
    die "Format: select_acey_annots_22oct2016.pl [infile] ;\n";
}

my $int_outfile         = $infile;
$int_outfile =~ s/\.tsv\.txt\z//;

my $int_sec_outfile     = $int_outfile;
my $non_int_outfile     = $int_outfile;
my $non_int_sec_outfile = $int_outfile;
my $non_spec_outfile    = $int_outfile;

$int_outfile         = "$int_outfile.intestine.tsv.txt";
$int_sec_outfile     = "$int_sec_outfile.intestine_secreted.tsv.txt";
$non_int_outfile     = "$non_int_outfile.non_intestine.tsv.txt";
$non_int_sec_outfile = "$non_int_sec_outfile.non_intestine_secreted.tsv.txt";
$non_spec_outfile    = "$non_spec_outfile.non_specific.tsv.txt";

$int_outfile         = safename($int_outfile);
$int_sec_outfile     = safename($int_sec_outfile);
$non_int_outfile     = safename($non_int_outfile);
$non_int_sec_outfile = safename($non_int_sec_outfile);
$non_spec_outfile    = safename($non_spec_outfile);

my @non_spec_genes = ();

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A Gene \t /xms ) {
        $header = $input;
    }
    elsif ( $input =~ /\A ([^\t]+) \t (?: [^\t]* \t){4} ([^\t]*) \t (?: [^\t]* \t){19} ([^\t]*) \t ([^\t]*) \t .* \z/xms ) { 
        my $gene    = $1;
        my $phobius = $2;   
        my $log2fc  = $3;
        my $fdr     = $4;

        if ( looks_like_number($log2fc) and looks_like_number($fdr) and ( $fdr <= 0.1 ) ) {
            if ( $log2fc > 0 ) {
                if ( $phobius =~ /SigP/xms ) { 
                    $data_ref->{'int_sec'}->{$fdr}->{'gene'}->{$gene}->{'annot'} = $input;
                }
                $data_ref->{'int'}->{$fdr}->{'gene'}->{$gene}->{'annot'} = $input;
            }
            elsif ( $log2fc < 0 ) {
                if ( $phobius =~ /SigP/xms ) {
                    $data_ref->{'non_int_sec'}->{$fdr}->{'gene'}->{$gene}->{'annot'} = $input;
                }
                $data_ref->{'non_int'}->{$fdr}->{'gene'}->{$gene}->{'annot'} = $input;
            } 
        }
        else {
            push @non_spec_genes, $gene;
            $data_ref->{'non_spec'}->{'gene'}->{$gene}->{'annot'} = $input;
        }

    }
    else { 
        die "Cannot parse input line: $input\n";
    }
}
close $INFILE;

open my $INT_OUTFILE, '>', $int_outfile;
print $INT_OUTFILE "$header\n";
my @int_fdrs = sort { $a <=> $b } keys %{ $data_ref->{'int'} };
foreach my $int_fdr (@int_fdrs) {
    my @int_genes = sort keys %{ $data_ref->{'int'}->{$int_fdr}->{'gene'} };
    foreach my $int_gene (@int_genes) {
        my $annot = $data_ref->{'int'}->{$int_fdr}->{'gene'}->{$int_gene}->{'annot'};
        print $INT_OUTFILE "$annot\n";
    }
}
close $INT_OUTFILE;

open my $INT_SEC_OUTFILE, '>', $int_sec_outfile;
print $INT_SEC_OUTFILE "$header\n";
my @int_sec_fdrs = sort { $a <=> $b } keys %{ $data_ref->{'int_sec'} };
foreach my $int_sec_fdr (@int_sec_fdrs) {
    my @int_sec_genes = sort keys %{ $data_ref->{'int_sec'}->{$int_sec_fdr}->{'gene'} };
    foreach my $int_sec_gene (@int_sec_genes) {
        my $annot = $data_ref->{'int_sec'}->{$int_sec_fdr}->{'gene'}->{$int_sec_gene}->{'annot'};
        print $INT_SEC_OUTFILE "$annot\n";
    }
}
close $INT_SEC_OUTFILE;

open my $NON_INT_OUTFILE, '>', $non_int_outfile;
print $NON_INT_OUTFILE "$header\n";
my @non_int_fdrs = sort { $a <=> $b } keys %{ $data_ref->{'non_int'} };
foreach my $non_int_fdr (@non_int_fdrs) {
    my @non_int_genes = sort keys %{ $data_ref->{'non_int'}->{$non_int_fdr}->{'gene'} };
    foreach my $non_int_gene (@non_int_genes) {
        my $annot = $data_ref->{'non_int'}->{$non_int_fdr}->{'gene'}->{$non_int_gene}->{'annot'};
        print $NON_INT_OUTFILE "$annot\n";
    }
}
close $NON_INT_OUTFILE;
 
open my $NON_INT_SEC_OUTFILE, '>', $non_int_sec_outfile;
print $NON_INT_SEC_OUTFILE "$header\n";
my @non_int_sec_fdrs = sort { $a <=> $b } keys %{ $data_ref->{'non_int_sec'} };
foreach my $non_int_sec_fdr (@non_int_sec_fdrs) {
    my @non_int_sec_genes = sort keys %{ $data_ref->{'non_int_sec'}->{$non_int_sec_fdr}->{'gene'} };
    foreach my $non_int_sec_gene (@non_int_sec_genes) {
        my $annot = $data_ref->{'non_int_sec'}->{$non_int_sec_fdr}->{'gene'}->{$non_int_sec_gene}->{'annot'};
        print $NON_INT_SEC_OUTFILE "$annot\n";
    }
}
close $NON_INT_SEC_OUTFILE;

open my $NON_SPEC_OUTFILE, '>', $non_spec_outfile;
print $NON_SPEC_OUTFILE "$header\n";
foreach my $gene (@non_spec_genes) {
    my $annot = $data_ref->{'non_spec'}->{'gene'}->{$gene}->{'annot'};
    print $NON_SPEC_OUTFILE "$annot\n";
}
close $NON_SPEC_OUTFILE;

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


