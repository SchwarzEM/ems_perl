#!/usr/bin/env perl

# get_cdhit_est_clusters.pl -- Erich Schwarz <ems394@cornell.edu>, 11/8/2013.
# Purpose: given a FASTA sequence file subjected to CD-HIT-EST, and the *.clustr output from that run, extract (and, as needed, revcomp) sequences of the clusters.

# TODO: make this thing print each cluster to a unique sequence file.

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $fasta   = q{};
my $clustr  = q{};
my $number = 0;

my $help;

my $seqname     = q{};
my $comment     = q{};
my $sequence    = q{};
my $orientation = q{};

my $key_seq   = q{};
my $other_seq = q{};
my $clust_no  = q{};

my $data_ref;

GetOptions ( 'fasta=s'    => \$fasta,
             'cluster=s'  => \$clustr,
             'number=i'   => \$number,
             'help'       => \$help,   );

if ( $help or (! $fasta) or (! $clustr) or (! looks_like_number($number) ) or ( $number != int($number) ) or ( $number < 2 ) ) {
    die "Format: get_cdhit_est_clusters.pl\n",
        "            --fasta|-f    [input FASTA file]\n",
        "            --cluster|-c  [input *.clustr file]\n",
        "            --number|-n   [minimum number of sequences in a cluster (must be integer, and at least 2)]\n",
        "            --help|-h     [print this message]\n",
        ;
}

open my $FASTA, '<', $fasta or die "Can't open FASTA file $fasta: $!";
while (my $input = <$FASTA>) {
    chomp $input;
    if ( $input =~ /\A > (\S+) (.*) \z/xms ) { 
        $seqname = $1;
        $comment = $2;
        if ( exists $data_ref->{'seqname'}->{$seqname} ) {
            die "In FASTA file $fasta, redundant sequence name $seqname\n";
        }
        $data_ref->{'seqname'}->{$seqname}->{'comment'} = $comment;
    }
    elsif ( $input =~ /\A > /xms ) { 
        die "Can't parse header line from FASTA file $fasta: $input\n";
    }
    elsif ( $input =~ / \S /xms ) { 
        $input =~ s/\s//g;
        $data_ref->{'seqname'}->{$seqname}->{'sequence'} .= $input;
    }
    else { 
        die "From FASTA file $fasta, can't parse text line $input\n";
    }
}        
close $FASTA or die "Can't close filehandle to FASTA file $fasta: $!";

open my $CLUSTR, '<', $clustr or die "Can't open *.clustr file $clustr: $!";
while (my $input = <$CLUSTR>) {
    chomp $input;  
    if ( $input =~ /\A > Cluster [ ] (\d+) \s* \z/xms ) {
        $clust_no = $1;
    }
    elsif ( $input =~ /\A \d+ \s+ \d+ nt , [ ] > (\S+) \.\.\. [ ] [*] \s* \z/xms ) { 
        $key_seq = $1;
        if (! exists $data_ref->{'seqname'}->{$key_seq} ) { 
            die "FASTA file $fasta does not contain key sequence $key_seq of cluster number $clust_no\n";
        }
        $data_ref->{'clust_no'}->{$clust_no}->{'key_seq'} = $key_seq;

        # By default, '+' -- but we note this explicitly, so that later code will Just Work all the time.
        $data_ref->{'seqname'}->{$key_seq}->{'orientation'} = '+';
    }
    elsif ( $input =~ /\A \d+ \s+ \d+ nt , [ ] > (\S+) \.\.\. [ ] at [ ] ([+]|[-]) \/ \d.+\d [%] \s* \z/xms ) {
        $other_seq   = $1;
        $orientation = $2;
        if (! exists $data_ref->{'seqname'}->{$other_seq} ) {
            die "FASTA file $fasta does not contain other sequence $other_seq of cluster number $clust_no\n";
        }
        $data_ref->{'clust_no'}->{$clust_no}->{'other_seq'}->{$other_seq} = 1;
        $data_ref->{'seqname'}->{$other_seq}->{'orientation'}             = $orientation;
    }
    else { 
        die "From *.clustr file $clustr, can't parse text line: $input\n";
    }
}
close $CLUSTR or die "Can't close filehandle to *.clustr file $clustr: $!";

my @clusters = grep { ( exists $data_ref->{'clust_no'}->{$_}->{'other_seq'} ) } sort { $a <=> $b } keys %{ $data_ref->{'clust_no'} };
foreach my $cluster (@clusters) { 
    $key_seq = $data_ref->{'clust_no'}->{$cluster}->{'key_seq'};
    my @other_seqs = sort keys %{ $data_ref->{'clust_no'}->{$cluster}->{'other_seq'} };
    my $seqcount = @other_seqs;
    $seqcount++;
    if ( $seqcount >= $number ) { 
        my $outfile = $clustr . '_Cluster_' . "$cluster.fa";
        $outfile    = safename($outfile);
        open my $OUTFILE, '>', $outfile or die "Can't open output cluster file $outfile: $!";
        my @seqs_to_print = ();
        push @seqs_to_print, $key_seq;
        push @seqs_to_print, @other_seqs;
        foreach my $seq_to_print (@seqs_to_print) { 
            my $final_seq_to_print = $seq_to_print;
            $comment = $data_ref->{'seqname'}->{$seq_to_print}->{'comment'};
            $sequence = $data_ref->{'seqname'}->{$seq_to_print}->{'sequence'};
            if ( $data_ref->{'seqname'}->{$seq_to_print}->{'orientation'} eq '-' ) { 
                $sequence = revcomp($sequence);
                $comment = $comment . ' [reverse complement of original]';
                $final_seq_to_print = $final_seq_to_print . '_RC';
            }
            my @output_lines = ();
            print $OUTFILE '>', $final_seq_to_print, $comment, "\n", ;
            @output_lines = unpack("a60" x (length($sequence)/60 + 1), $sequence);
            foreach my $output_line (@output_lines) {
                if ($output_line =~ /\S/) {
                    print $OUTFILE "$output_line\n";
                }
            }
        }
        close $OUTFILE or die "Can't close filehandle to output cluster file $outfile: $!";
    }
}

sub revcomp { 
    my $input_sequence = $_[0];
    $input_sequence = reverse($input_sequence);
    if ( $input_sequence =~ /[^acgtnACGTN]/xms ) { 
        die "Not currently designed to parse any letters but a, c, g, t, n, A, C, G, T, or N!\n";
    }
    $input_sequence =~ tr/acgtACGT/tgcaTGCA/;
    return $input_sequence;
}

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

