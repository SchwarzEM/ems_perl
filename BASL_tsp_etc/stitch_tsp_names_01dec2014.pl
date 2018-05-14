#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Getopt::Long;

my @genes2uniprots = ();
my @new2orig_names = ();
my @fastas         = ();
my $help;

my $data_ref;

GetOptions ( 'genes2uniprots=s{,}'  => \@genes2uniprots,
             'new2orig_names=s{,}'  => \@new2orig_names,
             'fastas=s{,}'          => \@fastas,
             'help'                 => \$help, 
);


if ( $help or (! @genes2uniprots) or (! @new2orig_names) or (! @fastas) ) { 
    die "Format: stitch_tsp_names_01dec2014.pl\n",
        "    --genes2uniprots|-g   [1+ precomputed UniProt/PFAM tables from PFAM]\n",
        "    --new2orig_names|-n   [1+ new names linked to original header text files; if using several, make sure no ambiguous short names!]\n",
        "    --fastas|-f           [1+ FASTA files (can be plain sequences, or aligned sequences]\n",
        "    --help|-h             [print this message]\n",
        ;
}

foreach my $genes2uniprot (@genes2uniprots) {
    open my $GENES, '<', $genes2uniprot;
    while (my $input = <$GENES>) {
        chomp $input;
        # Sample input for @genes2uniprots:
        # LBM   Q24188
        # TM4SF O76137
        # TSP26A        Q9VMJ6
        chomp $input;
        if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
            my $gene    = $1;
            my $uniprot = $2;
            if ( ( exists $data_ref->{'uniprot'}->{$uniprot}->{'gene'} ) and ( $gene ne $data_ref->{'uniprot'}->{$uniprot}->{'gene'} ) ) {
                die "In genes-to-Uniprot file $genes2uniprot, UniProt ID $uniprot mapped to two genes, $data_ref->{'uniprot'}->{$uniprot}->{'gene'} and $gene\n";
            }
            $data_ref->{'uniprot'}->{$uniprot}->{'gene'} = $gene;
        }
        else {
            die "In genes-to-Uniprot file $genes2uniprot, can't parse input line: $input\n";
        }
    }
    close $GENES;
}

foreach my $new2orig_name (@new2orig_names) {
    open my $ORIG, '<', $new2orig_name;
    while (my $input = <$ORIG>) {
        chomp $input;
        # Sample input for @new2orig_names:
        # >A1Z6T8_DROME_8-219	tr|A1Z6T8|A1Z6T8_DROME/8-219
        # >A1Z6U3_DROME_7-203	tr|A1Z6U3|A1Z6U3_DROME/7-203
        # >A1Z6U7_DROME_8-214	tr|A1Z6U7|A1Z6U7_DROME/8-214
        if ( $input =~ /\A > /xms ) { 
            if ( $input =~ /\A > ( \S+ _ [A-Z]+ _ \d+ [-] \d+) \s+ [a-z]+ \| ([^\|\s]+) \| /xms ) {
                my $seqname   = $1;
                my $uniprot   = $2;
                if ( ( exists $data_ref->{'seqname'}->{$seqname}->{'uniprot'} ) and ( $data_ref->{'seqname'}->{$seqname}->{'uniprot'} ne $uniprot ) ) {
                    die "In new-to-orig names file $new2orig_name, sequence $seqname is mapped to two UniProt IDs, $data_ref->{'seqname'}->{$seqname}->{'uniprot'} and $uniprot\n";
                }
                $data_ref->{'seqname'}->{$seqname}->{'uniprot'} = $uniprot;
            }
            else {
                die "In new-to-orig names file $new2orig_name, cannot parse: $input\n";
            }
        }
    }
    close $ORIG;
}

foreach my $fasta (@fastas) {
    open my $FASTA, '<', $fasta;
    while (my $input = <$FASTA>) {
        chomp $input;
        # Sample input for @fastas:
        # >LBM_DROME_14-197
        # >A1Z6T8_DROME_8-219
        # >A1Z6U3_DROME_7-203
        if ( $input =~ /\A > /xms ) {
            if ( $input =~ /\A > ( \S+ (_ [A-Z]+ _ \d+ [-] \d+)) (.*) \z/xms ) {
                my $seqname   = $1;
                my $seqsuffix = $2;
                my $comments  = $3;
                if (! exists $data_ref->{'seqname'}->{$seqname}->{'uniprot'} ) {
                    print "$input\n";
                }
                else {
                    my $uniprot = $data_ref->{'seqname'}->{$seqname}->{'uniprot'};
                    if (! exists $data_ref->{'uniprot'}->{$uniprot}->{'gene'} ) {
                        print "$input\n";
                    }
                    else { 
                        my $gene = $data_ref->{'uniprot'}->{$uniprot}->{'gene'};
                        my $new_header = '>' . $gene . $seqsuffix . $comments;
                        print "$new_header\n";
                    }
                }
            }
            else {
                die "In FASTA file $fasta, cannot parse: $input\n";
            }
        }
        else {
            print "$input\n";
        }
    }   
    close $FASTA;
}

