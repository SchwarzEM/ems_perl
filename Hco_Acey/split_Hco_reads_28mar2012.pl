#!/usr/bin/env perl

# split_Hco_reads_28mar2012.pl -- Erich Schwarz <emsch@caltech.edu>, 3/29/2012.
# Purpose: given many jumbled Hco reads, split them into defined groups based on their read headers (for Sujai filtering).

use strict;
use warnings;

my $data_ref;
 
my $target = q{};

my %header2library = ( 'HWUSI-EAS627:3'       => '300nt_42FC3AAXX.3',
                       'HWUSI-EAS627_0001:6'  => '300nt_42YPHAAXX.6',  
                       'HWUSI-EAS627_0001:1'  => '300nt_61J8HAAXX.1',
                       'HWUSI-EAS627_0012'    => '300nt_62636AAXX.1',
                       'HWUSI-EAS627_0017'    => '500nt_625LPAAXX.4',
                       'HWUSI-EAS627_0014'    => '500nt_6263DAAXX.8',
                       'ILLUMINA-33A494_0009' => '500nt_62DJDAAXX.8',
                       'HWI-ST0787_0069'      => '500nt_B08VPABXX.3', 
                       'A80PGJABXX:5'         => '2kb_Hco_02aug2011a',
                       'A812R5ABXX:7'         => '2kb_Hco_02aug2011b',
                       'A819JDABXX:6'         => '2kb_Hco_02aug2011c',
                       'A819JDABXX:7'         => '2kb_Hco_02aug2011d',
                       'B810L6ABXX:5'         => '2kb_Hco_02aug2011e',
                       'A812R5ABXX:6'         => '5kb_Hco_02aug2011',
                       'B810K5ABXX:1'         => '5kb_Hco_23jan2011',
                       'B80PGKABXX:5'         => '10kb_Hco_23jan2011',
                       'HWUSI-EAS627:8'       => 'single_30LNFAAXX.8',
                       'HWUSI-EAS627:7'       => 'single_311B7AAXX.7', 
);

my @read_labels = sort keys %header2library;

while (my $input = <> ) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) .* \z/xms ) {
        my $read_name = $1;
        my $assigned  = 0;
        foreach my $read_label (@read_labels) { 
            if ( (! $assigned ) and ( $read_name =~ /\A $read_label /xms ) ) { 
                $target = $read_label;
                push @{ $data_ref->{'reads'}->{$target} }, $input;
                $assigned = 1;
            }
        }
        if (! $assigned) { 
            warn "Can't parse read header: $input\n";
            $target = {};
        }
    }
    else { 
        if ($target) { 
            push @{ $data_ref->{'reads'}->{$target} }, $input;
        }
    }
}

foreach my $read_label (@read_labels) { 
    if ( exists $data_ref->{'reads'}->{$read_label} ) { 
        my $output_tag = $header2library{$read_label};
        my $output     = $output_tag . '.split_Hco_28mar2012.fa';
        $output        = safename($output);
        open my $OUTPUT, '>', $output or die "Can't print to output file $output: $!";
        foreach my $seq_line ( @{ $data_ref->{'reads'}->{$read_label} } ) { 
            print $OUTPUT "$seq_line\n";
        }
        close $OUTPUT or die "Can't close filehandle to to output file $output: $!";
    }
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

