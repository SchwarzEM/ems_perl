#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $orthologs = q{};
my @taxa      = ();
my @rseqs     = ();
my $log2fc    = q{};
my $fdr       = q{};
my $help;

GetOptions ( 'orthologs=s' => \$orthologs,
             'taxa=s{,}'   => \@taxa,
             'rseqs=s{,}'  => \@rseqs,
             'log2fc=s'    => \$log2fc,
             'fdr=s'       => \$fdr,
             'help'        => \$help,   );

if ( $help 
         or (! $orthologs) or (! @taxa) or (! @rseqs) or (! $log2fc) or (! $fdr)
         or (! looks_like_number($log2fc) ) or (! looks_like_number($fdr) ) 
         or ( $log2fc <= 0 ) or ( $fdr < 0 ) or ( $fdr > 1 ) 
    ) {
    die "Format: strict_ortholog_data_12feb2017.pl\n",
        "            --orthologs|-o [OrthoMCL-formatted data, as file or '-']\n",
        "            --taxa|-t      [list of 2+ taxa to extract as strict orthologs and RNA-seq data]\n",
        "            --rseqs|-r     [list of 2+ RNA-seq files from which to extract log2FC and FDR; must match taxa list]\n",
        "            --log2fc|-l    [minimum threshold for bias; must be positive real number]\n",
        "            --fdr|-f       [maximum allowable FDR; must be real number between 0 and 1]\n",
        "            --help|-h      [print this message and exit]\n",
}

my $taxon_count = @taxa;
my $rseq_count  = @rseqs;

if ( ( $taxon_count != $rseq_count ) or ( $taxon_count <= 1 ) ) {
    die "Taxon count must equal rseq count; both must be >= 2\n";
}

# Deal with zero-based indexing of arrays.
$taxon_count--;

# For each taxon, record from its corresponding Rseq data file whether genes are up-biased or down-biased.
# Keep records of the actual expression values so that they can be read later.
foreach my $i (0..$taxon_count) {
    my $taxon = $taxa[$i];
    my $rseq  = $rseqs[$i];

    # Later, I will want fast way to see if a taxon is key.
    $data_ref->{'taxon'}->{$taxon}->{'key'} = 1;

    open my $RSEQ, '<', $rseq;
    while (my $input = <$RSEQ>) {
        chomp $input;
        if ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \z/xms ) {
            my $gene       = $1;
            my $log2fc_val = $2;
            my $fdr_val    = $3;
            # Note that we want *either* up-biased *or* down-biased, so we compare abs(log2FC), not log2FC, to our threshold.
            if ( ( $gene ne 'Gene' ) and ( abs($log2fc_val) >= $log2fc ) and ( $fdr_val <= $fdr ) ) {
                $data_ref->{'taxon'}->{$taxon}->{'gene'}->{$gene}->{'log2FC'} = $log2fc_val;
                $data_ref->{'taxon'}->{$taxon}->{'gene'}->{$gene}->{'FDR'} = $fdr_val;
                if ( $log2fc_val > 0 ) {
                    $data_ref->{'taxon'}->{$taxon}->{'gene'}->{$gene}->{'status'} = 'up_biased';
                }
                elsif ( $log2fc_val < 0 ) {
                    $data_ref->{'taxon'}->{$taxon}->{'gene'}->{$gene}->{'status'} = 'down_biased';
                }
            }
        }
        else {
            die "From rseq file \"$rseq\" of taxon \"$taxon\", cannot parse: $input\n";
        }
    }
}

# Build a header line for my outputs.
my $header = "Orthology_group";
foreach my $taxon (@taxa) {
    $header = $header . "\t". "$taxon.taxon\t$taxon.gene\t$taxon.status\t$taxon.log2FC\t$taxon.FDR";
}

my $ORTHOLOGS;
if ($orthologs eq '-') {
    # Special case: get the stdin handle
    $ORTHOLOGS = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $ORTHOLOGS, '<', $orthologs;
}
while (my $input = <$ORTHOLOGS>) { 
    # obligatory parsing of the input line
    if ($input =~ / \A
                    ([A-Z]+\d+)                   # $1 -> $orth_group -- no longer demands 'ORTHOMCL' as only permissible starting text 
                    \( \d+\sgenes,\d+\staxa \)    # just (punctuation)
                    : \s+ \b                      # Enforce word boundary just before start of $oseqs_line.

                    (\S.+\S)                      # $2 -> $oseqs_line
                    \s* 
                    \z 
                  /xms) {

        # Note that we refer to '$oseqs_line', not '$oprots_line', 
        #     because the OrthoMCL line may have either proteins or genes.

        my $orth_group = $1;
        my $oseqs_line = $2;
        my @orthoseqs  = split /\s+/, $oseqs_line;

        # Enforce only *one* instance of key taxon per line
        my %taxon2gene  = ();

        foreach my $o_seq (@orthoseqs) { 
            # Obligatory parsing of the $o_seq set:
            if ( $o_seq =~ / \A ([^\s\(\)]+) \( ( [^\s\(\)]+ ) \) \z /xms ) { 
                my $gene  = $1;
                my $taxon = $2;
                if ( ( exists $data_ref->{'taxon'}->{$taxon}->{'key'} ) and ( exists $taxon2gene{$taxon} ) ) {
                    die "Non-strict orthology for taxon \"$taxon\" in: $input\n";
                }
                $taxon2gene{$taxon} = $gene;
            }
            # Enforce correct format of $o_seq set:
            else {
                die "Can't parse sequence $o_seq in: $input\n";
            }
        }

        # Build output line:
        my $output = $orth_group;

        foreach my $taxon (@taxa) {
            my $gene   = $taxon2gene{$taxon};
            my $status = q{};
            my $log2fc = q{};
            my $fdr    = q{};
            if ( exists $data_ref->{'taxon'}->{$taxon}->{'gene'}->{$gene}->{'status'} ) {
                $status = $data_ref->{'taxon'}->{$taxon}->{'gene'}->{$gene}->{'status'};
                $log2fc = $data_ref->{'taxon'}->{$taxon}->{'gene'}->{$gene}->{'log2FC'};
                $fdr    = $data_ref->{'taxon'}->{$taxon}->{'gene'}->{$gene}->{'FDR'};
            }
            $output = $output . "\t" . "$taxon\t$gene\t$status\t$log2fc\t$fdr";
       }

       print "$header\n" if $header;
       $header = q{};

       print "$output\n";
    }
    else {
        die "From orthology file $orthologs, cannot parse: $input\n";
    }
}
close $ORTHOLOGS;
