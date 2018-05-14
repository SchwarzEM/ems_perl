#!/usr/bin/env perl

# parse_listed_rsem_summs_10feb2016.pl -- Erich Schwarz <ems394@cornell.edu>, 2/10/2016.
# Purpose: given table listing RNA-seq data types and corresponding RSEM readmap summaries, open and parse each report w/ note of data type.

use strict;
use warnings;
use autodie;

use File::Basename;

my $data_ref;

my @types = ();

my $header = "Data\tInput reads\tUnmapped\tUnmapped percent\tUniquely mapped\tUniq. percent\tMultiply mapped\tMulti. percent\tTotal percent";

while (my $report = <>) {
    chomp $report;
    my $basename = basename($report);
    if ( $basename =~ /\A job_rsem_sepals_ (\S+_rep\d) \.trim_exact_3nt_2016.01.21.\d{2}.e\d{8} \z/xms ) {
        my $type = $1;
        if ( exists $data_ref->{'type'}->{$type}->{'report'} ) {
            die "Data type $type given two different reports: $data_ref->{'type'}->{$type}->{'report'} and $report\n";
        }
        if (! -r $report ) {
            die "For data type $type, unreadable putative report file: $report\n";
        }
        push @types, $type;
        $data_ref->{'type'}->{$type}->{'report'} = $report;
    }
    else {
        die "Cannot parse input line: $report\n";
    }
}

foreach my $type (@types) {
    my $report = $data_ref->{'type'}->{$type}->{'report'};
    open my $REPORT, '<', $report;
    my $input_reads   = q{};
    my $zero_mapped   = q{};
    my $zero_percent  = q{};
    my $uniq_mapped   = q{};
    my $uniq_percent  = q{};
    my $multi_mapped  = q{};
    my $multi_percent = q{};
    my $total_percent = q{};
    while (my $input = <$REPORT>) {
        chomp $input;

        # Sample input line:
        # 33223373 reads; of these:
        if ( $input =~ /\A (\d+) \s+ reads; \s+ of \s+ these: \s* \z/xms ) { 
            $input_reads = $1;
            $input_reads = commify($input_reads);
        }

        # Sample input line:
        #     4437975 (13.36%) aligned 0 times
        elsif ( $input =~ /\A \s+ (\d+) \s+ \( (\S+) \) \s+ aligned \s+ 0 \s+ times \s* \z/xms ) {
            $zero_mapped  = $1;
            $zero_percent = $2;
            $zero_mapped  = commify($zero_mapped);
        }

        # Sample input line:
        #     16564126 (49.86%) aligned exactly 1 time
        elsif ( $input =~ /\A \s+ (\d+) \s+ \( (\S+) \) \s+ aligned \s+ exactly \s+ 1 \s+ time \s* \z/xms ) {
            $uniq_mapped  = $1;
            $uniq_percent = $2;
            $uniq_mapped  = commify($uniq_mapped);
        }

        # Sample input line:
        #     12221272 (36.79%) aligned >1 times
        elsif ( $input =~ /\A \s+ (\d+) \s+ \( (\S+) \) \s+ aligned \s+ [>] 1 \s+ times \s* \z/xms ) {
            $multi_mapped  = $1;
            $multi_percent = $2;
            $multi_mapped  = commify($multi_mapped);
        }

        # Sample input line:
        # 86.64% overall alignment rate
        elsif ( $input =~ /\A (\S+) \s+ overall \s+ alignment \s+ rate \s* \z/xms ) {
            $total_percent = $1;

            # Print header exactly once at start of output:
            print "$header\n" if $header;
            $header = q{};

            # Print RSEM readmapping summary information for data type, but only when all data have been gathered:
            print "$type\t$input_reads\t$zero_mapped\t$zero_percent\t$uniq_mapped\t$uniq_percent\t$multi_mapped\t$multi_percent\t$total_percent\n";

            # Zero out all the values.
            $input_reads   = q{};
            $zero_mapped   = q{};
            $zero_percent  = q{};
            $uniq_mapped   = q{};
            $uniq_percent  = q{};
            $multi_mapped  = q{};
            $multi_percent = q{};
            $total_percent = q{};
        }
    }
    close $REPORT;
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

