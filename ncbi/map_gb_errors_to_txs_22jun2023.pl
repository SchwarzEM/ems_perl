#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $gb_report = q{};
my $genbank   = q{};

my %bad_peps  = ();
my %results   = ();

$gb_report = $ARGV[0] if $ARGV[0];
$genbank   = $ARGV[1] if $ARGV[1];

my $read_mRNAs = 0;
my @tx_names   = ();

my $read_peps = 0;
my @pep_names = ();

if ( (! $gb_report ) or (! $genbank ) ) {
   die "Format: map_gb_errors_to_txs_22jun2023.pl [GenBank error report] [GenBank of attempted *.sqn] > [rejected transcripts]\n";
}

open my $REPORT, '<', $gb_report;
while ( my $input = <$REPORT> ) {
    chomp $input;
    if ( $input =~ /\A ERROR: [ ] .+ \[ gnl \| (\S+) \| (cds-\d+) \] \s* \z/xms ) {
        my $locus_prefix = $1;
        my $cds_id       = $2;
        my $pep = "$locus_prefix:$cds_id";
        $bad_peps{$pep} = 1;
    }
}
close $REPORT;

open my $GENBANK, '<', $genbank;
while ( my $input = <$GENBANK> ) {
    chomp $input;
    # Add check for '//', because the very last gene in the file needs to be parsed too!
    if ( ( $input =~ /\A \s+ gene \s+ /xms ) or ( $input =~ /\A \/\/ /xms ) ) {
        # Record anything taken from a previous gene, if it matches a bad peptide.
        if (@pep_names) {
            my $count = @pep_names;
            $count--;
            for my $i (0..$count) {
                my $pep = $pep_names[$i];
                my $tx  = $tx_names[$i];
                if ( exists $bad_peps{$pep} ) {
                    my $data = "$pep\t$tx";
                    $results{$data} = 1;
                    delete $bad_peps{$pep};
                }
            }
        }
        # Then zero these variables out (again).
        $read_mRNAs = 0;
        @tx_names   = ();

        $read_peps = 0;
        @pep_names = ();
    }
    elsif ( $input =~ /\A \s+ mRNA \s+ /xms ) {
        $read_mRNAs = 1;
        $read_peps  = 0;
    }
    elsif ( $input =~ /\A \s+ CDS \s+ /xms ) {
        $read_mRNAs = 0;
        $read_peps  = 1;
    }
    elsif ( $read_mRNAs and ( $input =~ / \/ note = \" (\S+) \" /xms )  ) {
        my $tx_id  = $1;
        push @tx_names, $tx_id;
    }
    elsif ( $read_peps and ( $input =~ / \/ protein_id = \" (\S+) \" /xms ) ) {
        my $pep_id = $1;
        push @pep_names, $pep_id;
    }
}
close $GENBANK;

# Print final mapping of bad translations to transcript names:
my @final_results = sort keys %results;
foreach my $final_result (@final_results) {
    print "$final_result\n";
}

# Point out if this left any bad translations unaccounted for:
my @remaining_bad_peps = sort keys %bad_peps;
foreach my $remaining_bad_pep (@remaining_bad_peps) {
    print "Failed to account for: $remaining_bad_pep\n";
}

