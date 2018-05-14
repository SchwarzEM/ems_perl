#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my $txs  = q{};
my $orth = q{};
my $gff3 = q{};
my $pref = q{};

$txs  = $ARGV[0] if $ARGV[0];
$orth = $ARGV[1] if $ARGV[1];
$gff3 = $ARGV[2] if $ARGV[2];
$pref = $ARGV[3] if $ARGV[3];

my $header = "Gene\t$pref" . "_Tx\t$pref" . "_exon_sum\t$pref" . "_intron_sum";

if ( (! -r $txs) or (! -r $orth) or (! -r $gff3) ) { 
    die "Format: exon.intron_sizes_of_chosen_txs.pl [ok_tx_list] [orthofinder_tx_table] [gff3] [data prefix] > [output]\n";
}

my $data_ref;

open my $TXS, '<', $txs;
while (my $input = <$TXS>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) { 
        $data_ref->{'tx'}->{$input}->{'OK'} = 1;
    }
    else {
        die "In transcript list $txs, cannot parse: $input\n";
    }
}
close $TXS;

open my $ORTH, '<', $orth;
while (my $input = <$ORTH>) {
    chomp $input;
    if ( $input =~ /\A (OG\d+) [:] \s+ (\S .+ \S) \s* \z/xms ) {
        my $ofind   = $1;
        my $tx_text = $2;
        my @tx_list = split /\s+/, $tx_text;
        foreach my $tx (@tx_list) {
            if ( exists $data_ref->{'tx'}->{$tx}->{'OK'} ) {
                $data_ref->{'ofind'}->{$ofind}->{'tx'} = $tx;
                $data_ref->{'tx'}->{$tx}->{'ofind'}    = $ofind;
            }
        }
    }
    else {
        die "In OrthoFind transcript file $orth, cannot parse: $input\n";
    }
}
close $ORTH;

open my $GFF3, '<', $gff3;
while (my $input = <$GFF3>) {
    chomp $input;
    if ( $input =~ /\A \S+ 
                       \t [^\t]* 
                       \t (CDS|intron) 
                       \t (\d+) 
                       \t (\d+) 
                       \t [^\t]* \t [^\t]* \t [^\t]* 
                       \t [^\t]* Parent= (?: Transcript[:]){0,1} (\S+) /xms ) { 
        my $seq_type = $1;
        my $start_nt = $2;
        my $end_nt   = $3;
        my $tx       = $4;

        if ( $end_nt < $start_nt ) {
            die "In GFF3 file $gff3, cannot parse nt interval $start_nt", '-', "$end_nt: $input\n";
        }

        if ( exists $data_ref->{'tx'}->{$tx}->{'OK'} ) {
            if (! exists $data_ref->{'tx'}->{$tx}->{'ofind'} ) {
                die "Cannot assign transcript $tx to an OrthoFinder group\n";
            }
            my $size = (($end_nt - $start_nt) + 1);
            push @{ $data_ref->{'annot_tx'}->{$tx}->{'seq_type'}->{$seq_type}->{'size'} }, $size;
        }
    }
}
close $GFF3;

my @ofind_groups    = sort keys %{ $data_ref->{'ofind'} };
my @exon_sum_lens   = ();
my @intron_sum_lens = ();

foreach my $ofind (@ofind_groups) { 
    my $annot_tx = $data_ref->{'ofind'}->{$ofind}->{'tx'};
    my @indiv_exon_lengths = ();
    if ( exists $data_ref->{'annot_tx'}->{$annot_tx}->{'seq_type'}->{'CDS'}->{'size'} ) {
        @indiv_exon_lengths = @{ $data_ref->{'annot_tx'}->{$annot_tx}->{'seq_type'}->{'CDS'}->{'size'} };
    }

    my @indiv_intron_lengths = ();
    if ( exists $data_ref->{'annot_tx'}->{$annot_tx}->{'seq_type'}->{'intron'}->{'size'} ) {
        @indiv_intron_lengths = @{ $data_ref->{'annot_tx'}->{$annot_tx}->{'seq_type'}->{'intron'}->{'size'} };
    }

    my $exon_sum_length   = 0;
    my $intron_sum_length = 0;

    if (@indiv_exon_lengths) {
        $exon_sum_length = sum(@indiv_exon_lengths);
    }

    if (@indiv_intron_lengths) {
        $intron_sum_length = sum(@indiv_intron_lengths);
    }

    push @exon_sum_lens, $exon_sum_length;
    push @intron_sum_lens, $intron_sum_length;

    print "$header\n" if $header;
    $header = q{};

    # Do *not* commify $exon_sum_length or $intron_sum_length, to avoid confounding downstream regression analyses.

    print "$ofind\t$annot_tx\t$exon_sum_length\t$intron_sum_length\n";
}

my $total_exon_len   = sum(@exon_sum_lens);
my $total_intron_len = sum(@intron_sum_lens);

# Commify the overall totals, for readability.

$total_exon_len   = commify($total_exon_len);
$total_intron_len = commify($total_intron_len);

print "Total\tAll_txs\t$total_exon_len\t$total_intron_len\n";


# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

