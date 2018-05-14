#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Getopt::Long;

my $data_ref;

my @cons_data  = ();
my $motif_data = q{};

my $data_set   = q{};
my $motif      = q{};
my $consensus  = q{};

my $added_head = "Consensus\tRev_comp_consensus";

my $help;

GetOptions ( 'cons_data=s{,}' => \@cons_data,
             'mot_data=s'     => \$motif_data,
             'help'           => \$help,        );

if ($help or (! @cons_data) or (! $motif_data) ) {
    die "Format: append_motif_consensi_12mar2016.pl\n",
        "            -c|--cons_data  [1+ MEME files, w\ full filenames, identifying data types and consensi]\n",
        "            -m|--mot_data   [one complex motif data table to be expanded]\n",
        "            -h|--help       [print this message]\n",
        "            [print complex annotation table to <STDOUT>]\n",
        ;
}

foreach my $cons_data_file (@cons_data) { 
    # Enforce identification of data set from reasonably full filename.
    # E.g.: ../adrienne_list1_08mar2016/adrienne_list1_08mar2016_500trans_meme_2016.03.07/meme.txt
    if ( $cons_data_file =~ / ([^\/\s]+) _500trans_meme_2016\.03\.07 \/meme\.txt \z/xms ) {
        $data_set = $1;
    }
    else {
        die "Cannot extract data set from MEME file name: $cons_data_file\n";
    }

    open my $CONS_DAT, '<', $cons_data_file;
    while (my $input = <$CONS_DAT>) {
        chomp $input;

        # Sample data:
        # MOTIF  1 MEME   width =  11  sites =  10  llr = 117  E-value = 2.9e-003
        # [...]
        # Multilevel           TCACCTACCTA
        # consensus            C    A  GC 
        # sequence                        

        if ( $input =~ /\A MOTIF \s+ (\S+) \s+ MEME \s+ width \s+ [=] \s+ \d+ \s+ 
                           sites \s+ [=] \s+ \d+ \s+ llr \s+ [=] \s+ \d+ \s+ E-value \s+ [=] \s+ \S+ \s* \z/xms ) { 
            $motif      = $1;
        }
        elsif ( $data_set and $motif and ( $input =~ /\A Multilevel \s+ ([ACGT]+) \s* \z/xms ) ) {
            $consensus  = $1;
            if ( exists $data_ref->{'data_set'}->{$data_set}->{'motif'}->{$motif}->{'consensus'} ) {
                warn "In consensus data file $cons_data_file, for data set \"$data_set\" and motif \"$motif\", redundant consensi:\n";
                warn "$data_ref->{'data_set'}->{$data_set}->{'motif'}->{$motif}->{'consensus'}\n";
                warn "vs. $consensus\n";
                die  "in input line: $input\n";
            }
            else {
                $data_ref->{'data_set'}->{$data_set}->{'motif'}->{$motif}->{'consensus'} = $consensus;
            }
            $motif      = q{};
            $consensus  = q{};
        }
    }
    close $CONS_DAT;
}

open my $MOT_DAT, '<', $motif_data;
while (my $input = <$MOT_DAT>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) .* \z/xms ) {
        my $data_set1 = $1;
        my $motif1    = $2;

        if ( $data_set1 eq 'Data' ) {
            print "$input\t$added_head\n";
        }
        elsif (! $data_ref->{'data_set'}->{$data_set1}->{'motif'}->{$motif1}->{'consensus'} ) { 
            die "Do not have consensus data for motif in: $input\n";
        }
        else {
            my $consensus     = $data_ref->{'data_set'}->{$data_set1}->{'motif'}->{$motif1}->{'consensus'};
            my $rev_comp_cons = revcomp($consensus);
            print "$input\t$consensus\t$rev_comp_cons\n";
        }
    }
    else {
        die "From motif data table $motif_data, cannot parse: $input\n";
    }
}
close $MOT_DAT;


sub revcomp { 
    my $input_sequence = $_[0];
    $input_sequence = reverse($input_sequence);
    if ( $input_sequence =~ /[^acgtACGT]/xms ) { 
        die "Not currently designed to parse any",
            " letters but a, c, g, t, A, C, G, or T!\n",
            ;
    }
    $input_sequence =~ tr/acgtACGT/tgcaTGCA/;
    return $input_sequence;
}

