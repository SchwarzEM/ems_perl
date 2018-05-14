#!/usr/bin/env perl

use strict;
use warnings;

my $input_list        = $ARGV[0];
my $stat_file         = q{};
my $refinement_script = q{};

my $stat_summ = '08dec2010.stat.summ.txt';
my $line_comm = '08dec2010.line.comm.sh';

$stat_summ = safename($stat_summ);
$line_comm = safename($line_comm);

my %ontology2GOid = ( molecular_function => '[GO:0003674]',
                      biological_process => '[GO:0008150]',
                      cellular_component => '[GO:0005575]', );

open my $STATS, '>', $stat_summ or die "Can't open statistics summary $stat_summ: $!\n";

open my $LINES, '>', $line_comm or die "Can't open line-command file $line_comm: $!\n";
print $LINES '#/bin/bash', "\n\n";

open my $INPUT_LIST, '<', $input_list or die "Can't open input list file $input_list: $!";
while (my $infile = <$INPUT_LIST>) { 
    chomp $infile;
    open my $INFILE, '<', $infile or die "Can't open putative input file $infile: $!\n";
    while (my $input1 = <$INFILE>) { 
        chomp $input1;
        if ( $input1 =~ / \A \s+ \- \s+ (\/ home \/ schwarz \/ func_work \/ [^\/\s]+ \/ statistics.txt) \s* /xms ) { 
            $stat_file = $1;
            open my $STAT_FILE, '<', $stat_file or die "Can't open statistics file $stat_file: $!";
            print $STATS "\nSummary from $stat_file:\n\n";
            my $ontology = q{};
            my $ready    = 0;
            my $high_p   = q{};
            while (my $input2 = <$STAT_FILE>) { 
                chomp $input2;
                if ($ontology) { 
                    if ( ( $ready == 2 ) and ( $input2 =~ / \A \d\S* \s+ (\d\S*) \s* \z /xms ) ) { 
                        $high_p = $1;
                        print $STATS "$ontology $ontology2GOid{$ontology}: $high_p\n";

                        # Immediately zero everything out again:
                        $ontology = q{};
                        $ready    = 0;   
                        $high_p   = q{};
                    }
                    elsif ( ( $ready == 1 ) and ( $input2 =~ / \b underrepresentation \s+ overrepresentation \b /xms ) ) {
                        $ready = 2;
                    }
                    elsif ( ( $ready == 0 ) and ( $input2 =~ / \A global-test-statistics: \s* \z /xms ) ) { 
                        $ready = 1;
                    }
                }
                elsif (! $ontology) { 
                    if ( $input2 =~ /\A (molecular_function|biological_process|cellular_component) \s* \z /xms ) {
                        $ontology = $1;
                        $ready = 0;   
                    }
                }
            }
            print $STATS "\n";
            close $STAT_FILE or die "Can't close filehandle to statistics file $stat_file: $!";
        }
        elsif ( $input1 =~ / \A Run \s+ script \s+  (\/ home \/ schwarz \/ func_work \/ [^\/\s]+ \.allGO\.dir \/ refin-2010-12-08.sh) \s+ for \s+ refinement /xms ) { 
            $refinement_script = $1;
            print $LINES "$refinement_script 0.001 0.001 1 ;\n";
        }
    }
    close $INFILE or die "Can't close filehandle to putative input file $infile: $!\n";
}

close $INPUT_LIST or die "Can't close filehandle to input list file $input_list: $!";

print $LINES "\n";
close $LINES or die "Can't close filehandle to line-command file $line_comm: $!\n"; 
system 'chmod', '+x', $line_comm;

close $STATS or die "Can't close filehandle to statistics summary $stat_summ: $!\n";

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

