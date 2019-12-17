#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);
use Scalar::Util qw(looks_like_number);

my $ofinds        = q{};
my $degs          = q{};
my $taxon         = q{};
my $expr          = q{};
my $log_fc_thresh = q{};
my $fdr_thresh    = q{};

my $header = "Gene\tOrthol_DEG\tHomol_DEG";

my $data_ref;

$ofinds        = $ARGV[0] if $ARGV[0];
$degs          = $ARGV[1] if $ARGV[1];
$taxon         = $ARGV[2] if $ARGV[2];
$expr          = $ARGV[3] if $ARGV[3];
$log_fc_thresh = $ARGV[4] if $ARGV[4];
$fdr_thresh    = $ARGV[5] if $ARGV[5];

if ( ! $fdr_thresh ) {
    &help_message;
    die "No FDR threshold specified\n";
}

if ( ( $expr ne 'up' ) and ( $expr ne 'down' ) ) {
    &help_message;
    die	"Cannot interpret 'expr' argument, which must be 'up' or 'down': $expr\n";
}

if (! looks_like_number($log_fc_thresh) ) {
    &help_message;
    die "log_FC_threshold not numerical: $log_fc_thresh\n";
}

if (! looks_like_number($fdr_thresh) ) {
    &help_message;
    die "FDR threshold not numerical: $fdr_thresh\n";
}

if ( ( $fdr_thresh < 0 ) or ( $fdr_thresh > 1 ) ) {
    &help_message;
    die "FDR threshold outside 0-1 range: $fdr_thresh\n";
}

# Basic idea: put the data-filtering into the data-collection.
# Step 1: only record homologs whose DEG data meet our thresholds.

open my $DEGS, '<', $degs;
while (my $input = <$DEGS>) {
    chomp $input;
    if ( ( $input !~ /\A Gene /xms ) and ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \z/xms ) ) {
        my $gene   = $1;
        my $log_fc = $2;
        my $fdr    = $3;
        if ( (! looks_like_number($log_fc) ) or (! looks_like_number($fdr) ) ) {
            die "In DEG file $degs, non-numeric data: $input\n";
        }
	if ( $fdr <= $fdr_thresh ) {
            if (    ( ( $expr eq 'up' )   and ( $log_fc >= $log_fc_thresh ) )
                 or ( ( $expr eq 'down' ) and ( $log_fc <= $log_fc_thresh ) ) )	{
       	       	$data_ref->{'homolog'}->{$gene}->{'log_fc'} = $log_fc;
       	       	$data_ref->{'homolog'}->{$gene}->{'fdr'}    = $fdr;
       	    }
       	}
    }
}
close $DEGS;

# Step 2: only list *genes* if they have homologs we have recorded.

open my $OFINDS, '<', $ofinds;
while (my $input = <$OFINDS>) {
    chomp $input;
    if ( ( $input !~ /\A Gene /xms ) and ( $input =~ /\A (\S+) \t ([^\t]+) \z/xms ) ) {
        my $gene  = $1;
        my $ofind = $2;
        if ( $ofind =~ /(\S+)\($taxon\)/xms ) {
            while ( $ofind =~ /(\S+)\($taxon\)/xmsg ) {
                my $homolog = $1;
                if ( ( exists $data_ref->{'homolog'}->{$homolog}->{'log_fc'} ) and ( exists $data_ref->{'homolog'}->{$homolog}->{'fdr'} ) ) {
                    $data_ref->{'gene'}->{$gene}->{'homolog'}->{$homolog} = 1;
                    $data_ref->{'gene'}->{$gene}->{'ofind'} = $ofind;
                }
            }
        }
    }
}
close $OFINDS;

# Having done both filters: print qualified genes with qualified orthologies and DEG values.

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) {
    my @homologs  = sort keys %{ $data_ref->{'gene'}->{$gene}->{'homolog'} };
    my $ofind     = $data_ref->{'gene'}->{$gene}->{'ofind'};
    my $hom_data  = q{};
    my @hom_stats = ();

    foreach my $homolog (@homologs) {
        my $log_fc   = $data_ref->{'homolog'}->{$homolog}->{'log_fc'};
        my $fdr      = $data_ref->{'homolog'}->{$homolog}->{'fdr'};
        my $hom_stat = "$homolog [logFC $log_fc; FDR $fdr]";
        push @hom_stats, $hom_stat;
    }

    @hom_stats = sort @hom_stats;
    @hom_stats = uniq(@hom_stats);
    $hom_data  = join '; ', @hom_stats;

    print "$header\n" if $header;
    $header = q{};

    print "$gene\t$ofind\t$hom_data\n";
}

sub help_message {
    warn "Format: hco.expr_onto_acey.ofind_17dec2019.pl\n",
         "            [gene-Orthofinder groups]\n",
         "            [homolog-logFC, FDR]\n",
         "            [orthology group taxon being searched]\n",
         "            ['up' or 'down' expression]\n",
         "            [logFC threshold: for 'up', select all x >= logFC thr.; for 'down', select all x <= logFC thr.]\n",
         "            [FDR threshold: 0 <= x <= 1]\n",
         ;
}
