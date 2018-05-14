#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my @infiles       = ();
my $log2fc_thresh = q{};
my $fdr_thresh    = q{};
my $basename      = q{};

my $header = 'Infile' . "\t" . 'log2FC_thresh' . "\t" . 'FDR_thresh' . "\t" . 'Upreg_genes' . "\t" . 'Downreg_genes';

my $data_ref;

my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'log2fc=s'     => \$log2fc_thresh,
             'fdr=s'        => \$fdr_thresh,
             'help'         => \$help,   );

if (    $help 
     or (! @infiles) 
     or (! looks_like_number($log2fc_thresh) ) 
     or (! looks_like_number($fdr_thresh) ) 
     or ( $log2fc_thresh <= 0 ) 
     or ( $fdr_thresh < 0 ) 
     or ( $fdr_thresh > 1 ) 
   ) {
    die "Format: summarize_edgeR_data_07aug2017.pl \n",
        "    --infile|-i  <input stream/files>\n",
        "    --log2fc|-l  [log2 Fold Change threshold; selects greater/less or *equal* values; must be positive number > 0]\n",
        "    --fdr|-f     [False Discovery Rate threshold; selects smaller or *equal* FDR; must be number >= 0, <= 1]\n",
        "    --help|-h    [print this message]\n",
        ;
}

foreach my $infile (@infiles) { 
    if (! -r $infile) {
        die "Cannot read input file: $infile\n";
    }
    $basename = basename($infile);

    my $INPUT_FILE;
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
        $infile = basename $infile if $basename;
    }
    while (my $input = <$INPUT_FILE>) { 
        chomp $input;
        if ( ( $input !~ /\A Gene \t /xms ) and ( $input =~ /\A \S+ \t (\S+) \t (\S+) /xms ) ) { 
            my $log2fc = $1;
            my $fdr    = $2;

            if ( (! looks_like_number($log2fc) ) or (! looks_like_number($fdr) ) ) {
                die "For input file $infile, cannot parse numerical values in: $input\n";
            }

            if ( $fdr <= $fdr_thresh ) { 
                if ( $log2fc >= $log2fc_thresh ) { 
                    $data_ref->{'basename'}->{$basename}->{'upreg'}++;
                }
                elsif ( $log2fc <= (-1*$log2fc_thresh) ) {
                    $data_ref->{'basename'}->{$basename}->{'downreg'}++;
                }
            }
        }
    }
    close $INPUT_FILE;
}

my @base_infiles = sort keys %{ $data_ref->{'basename'} };

foreach my $base_infile (@base_infiles) {
    my $upreg_gene_count   = $data_ref->{'basename'}->{$base_infile}->{'upreg'};
    my $downreg_gene_count = $data_ref->{'basename'}->{$base_infile}->{'downreg'}++;

    $upreg_gene_count   = commify($upreg_gene_count);
    $downreg_gene_count = commify($downreg_gene_count);

    print "$header\n" if $header;
    $header = q{};

    print $base_infile, "\t", $log2fc_thresh, "\t", $fdr_thresh, "\t", $upreg_gene_count, "\t", $downreg_gene_count, "\n";

}

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}


