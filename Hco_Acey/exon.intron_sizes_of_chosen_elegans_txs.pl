#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Statistics::Descriptive;

my $txs  = q{};
my $gff3 = q{};

$txs  = $ARGV[0] if $ARGV[0];
$gff3 = $ARGV[1] if $ARGV[1];

if ( (! -r $txs) or (! -r $gff3) ) { 
    die "Format: exon.intron_sizes_of_chosen_aug_txs.pl [ok_tx_list] [gff3] > [output]\n";
}

my $data_ref;

open my $TXS, '<', $txs;
while (my $input = <$TXS>) {
    chomp $input;
    # One *very* weirdly protein-coding gene in C. elegans is E_BE45912.2 -- WBGene00017153!
    # Another gene: CE7X_3.1
    if ( $input =~ /\A ([a-zA-Z0-9_]+\.[0-9]+[a-z]*) \z/xms ) { 
        $data_ref->{'tx'}->{$input}->{'OK'} = 1;
    }
    else {
        die "In transcript list $txs, cannot parse: $input\n";
    }
}
close $TXS;

open my $GFF3, '<', $gff3;
while (my $input = <$GFF3>) {
    chomp $input;

    if ( $input =~ /\A \S+ 
                       \t WormBase
                       \t CDS
                       \t (\d+) 
                       \t (\d+) 
                       \t [^\t]* \t [^\t]* \t [^\t]* 
                       \t ID=CDS: ([a-zA-Z0-9_]+\.[0-9]+[a-z]*) (?:;|\b) /xms ) {
        my $start_nt = $1;
        my $end_nt   = $2;
        my $tx       = $3;

        if ( $end_nt < $start_nt ) {
            warn "In GFF3 file $gff3, cannot parse nt interval $start_nt", '-', "$end_nt\n";
            die "Offending line: $input\n";
        }

        if ( exists $data_ref->{'tx'}->{$tx}->{'OK'} ) {
            my $size = (($end_nt - $start_nt) + 1);
            push @{ $data_ref->{'annot_tx'}->{$tx}->{'CDS'}->{'size'} }, $size;
            $data_ref->{'annot_tx'}->{$tx}->{'CDS'}->{'count'}++;

            # This is a work-around for the failure to have 'intron' data in the remanei-Fierst GFF3;
            # first, the nt coordinates of all the neighboring nt of all the exons will be recorded; 
            # then they will be sorted; then they will be edited to yield pairs of intron coordinates.

            my $pre_start_nt = ($start_nt - 1);
            my $post_end_nt  = ($end_nt + 1);
            push @{ $data_ref->{'annot_tx'}->{$tx}->{'intron_coords'} }, ($pre_start_nt, $post_end_nt);
        }
    }

    # Enforce strict ID of CDS name for each putative WormBase CDS line
    elsif ( $input =~ /\A \S+ \t WormBase \t CDS \t /xms ) { 
        die "Cannot parse putative WormBase CDS GFF3 line: $input\n";
    }
}
close $GFF3;

my @annotated_txs   = sort keys %{ $data_ref->{'annot_tx'} };
my @exon_lengths    = (); 
my @intron_lengths  = ();
my @exon_sum_lens   = ();
my @intron_sum_lens = ();

foreach my $annot_tx (@annotated_txs) {
    my @indiv_exon_lengths = ();
    if ( exists $data_ref->{'annot_tx'}->{$annot_tx}->{'CDS'}->{'size'} ) {
        @indiv_exon_lengths = @{ $data_ref->{'annot_tx'}->{$annot_tx}->{'CDS'}->{'size'} };
    }

    my @indiv_intron_lengths = ();

    if ( $data_ref->{'annot_tx'}->{$annot_tx}->{'intron_coords'} ) {
        my @proto_intron_coords = @{ $data_ref->{'annot_tx'}->{$annot_tx}->{'intron_coords'} };
        @proto_intron_coords = sort { $a <=> $b } @proto_intron_coords;

        my $proto_intron_coord_text = join "\t", @proto_intron_coords;

        if ( $proto_intron_coord_text =~ /\A \d+ \t \d+ \z/xms ) { 
            $proto_intron_coord_text = q{};
        }

        if ( $proto_intron_coord_text =~ /\A \d+ \t .+ \d+ \z/xms ) {
            # Leave a leading-edge '\t', to make successive deletions easy
            $proto_intron_coord_text =~ s/\A\d+//;
            $proto_intron_coord_text =~ s/\t\d+\z//;

            while ( $proto_intron_coord_text =~ /\A \t (\d+) \t (\d+) (.*?) \z/xmsg ) {
                my $start_nt             = $1;
                my $end_nt               = $2;
                $proto_intron_coord_text = $3;

                my $size = (($end_nt - $start_nt) + 1);
                push @indiv_intron_lengths, $size;
            }
            # ensure that the intron coordinates are all gone
            if ( $proto_intron_coord_text !~ /\A \s* \z/xms ) {
                die "Failed to parse full intron coordinates in tx. $annot_tx: still have \"$proto_intron_coord_text\"\n";
            }
        }
    }

    my $exon_sum_length   = 0;
    my $intron_sum_length = 0;

    my $exon_count   = $data_ref->{'annot_tx'}->{$annot_tx}->{'CDS'}->{'count'};
    my $intron_count = ($exon_count-1);

    if (@indiv_exon_lengths) {
        $exon_sum_length = sum(@indiv_exon_lengths);
    }

    if (@indiv_intron_lengths) {
        $intron_sum_length = sum(@indiv_intron_lengths);
    }

    # sanity checks
    if ( ( $exon_sum_length < $exon_count ) or ( $exon_sum_length == 0 ) ) {
        die "Transcript $annot_tx putatively has exon sum length of $exon_sum_length despite having $exon_count exons\n";
    } 
    if ( $intron_sum_length < $intron_count ) {
        die "Transcript $annot_tx putatively has intron sum length of $intron_sum_length despite having $intron_count introns\n";
    }

    push @exon_lengths, @indiv_exon_lengths;
    push @intron_lengths, @indiv_intron_lengths;

    push @exon_sum_lens, $exon_sum_length;
    push @intron_sum_lens, $intron_sum_length;
}

my $exon_len_stat   = Statistics::Descriptive::Full->new();
my $intron_len_stat = Statistics::Descriptive::Full->new();

my $exon_sum_len_stat = Statistics::Descriptive::Full->new();
my $intron_sum_len_stat = Statistics::Descriptive::Full->new();

$exon_len_stat->add_data(@exon_lengths);
$intron_len_stat->add_data(@intron_lengths);

$exon_sum_len_stat->add_data(@exon_sum_lens);
$intron_sum_len_stat->add_data(@intron_sum_lens);

# round means to the nearest integer; commify everything.

my $annot_tx_count = @annotated_txs;
$annot_tx_count = commify($annot_tx_count);

my $exon_count     = @exon_lengths;
$exon_count = commify($exon_count);

my $intron_count   = @intron_lengths;
$intron_count = commify($intron_count);

my $mean_exon_len = $exon_len_stat->mean();
$mean_exon_len    = int(($mean_exon_len + 0.5));
$mean_exon_len = commify($mean_exon_len);

my $median_exon_len = $exon_len_stat->median();
$median_exon_len = commify($median_exon_len);

my $mean_intron_len = $intron_len_stat->mean();
$mean_intron_len    = int(($mean_intron_len + 0.5));
$mean_intron_len = commify($mean_intron_len);

my $median_intron_len = $intron_len_stat->median();
$median_intron_len = commify($median_intron_len);

my $mean_exon_sum_len   = $exon_sum_len_stat->mean();
$mean_exon_sum_len = int(($mean_exon_sum_len + 0.5));
$mean_exon_sum_len = commify($mean_exon_sum_len);

my $median_exon_sum_len = $exon_sum_len_stat->median();
$median_exon_sum_len = commify($median_exon_sum_len);

my $mean_intron_sum_len   = $intron_sum_len_stat->mean();
$mean_intron_sum_len = int(($mean_intron_sum_len + 0.5));
$mean_intron_sum_len = commify($mean_intron_sum_len);

my $median_intron_sum_len = $intron_sum_len_stat->median();
$median_intron_sum_len = commify($median_intron_sum_len);

print "Tx_count",

      "\tExon_count",
      "\tIntron_count",

      "\tMean_exon_sum_len",
      "\tMean_intron_sum_len",

      "\tMedian_exon_sum_len",  
      "\tMedian_intron_sum_len",

      "\tMean_exon_len",
      "\tMean_intron_len",

      "\tMedian_exon_len",
      "\tMedian_intron_len",

      "\n",
      ;

print "$annot_tx_count",

      "\t$exon_count",
      "\t$intron_count",

      "\t$mean_exon_sum_len",
      "\t$mean_intron_sum_len",

      "\t$median_exon_sum_len",
      "\t$median_intron_sum_len",

      "\t$mean_exon_len",
      "\t$mean_intron_len",

      "\t$median_exon_len",
      "\t$median_intron_len",

      "\n",
      ;

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

