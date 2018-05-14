#!/usr/bin/env perl

# renumber_wb_motifs.pl -- Erich Schwarz <emsch@caltech.edu>, 10/9/2012.
# Purpose: given a file or stream of WBMotif text in .ace format, and a designated starting number N, renumber the motifs starting from N (default is N=1).

use strict;
use warnings;
use Getopt::Long;

my @infiles = ();
my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'help'         => \$help, );

if ( (! @infiles) or $help ) { 
    die "Format: stormo2pfm_10oct2012.pl\n",
        "    --infile|-i   <input stream/files>\n",
        "    --help|-h     [print this message]\n",
        ;
}

my $i            = 1;
my $name         = q{};
my $consensus    = q{};
my $a_vals       = q{};
my $c_vals       = q{};
my $g_vals       = q{};
my $t_vals       = q{};
my $cluster_size = q{};

print "\n";

foreach my $infile (@infiles) { 
    open my $INPUT, '<', $infile or die "Can't open input file $infile. $!\n";
    while (my $input = <$INPUT>) { 
        chomp $input;
        if ( (! $name ) and ( $input =~ /\A (\S+) \z/xms ) ) { 
            $name = $1;
        }
        elsif ( $name and (! $consensus ) and ( $input =~ /\A Consensus \s+ [=] \s+ (\S+) \s*/xms ) ) { 
            $consensus = $1;
        }
        elsif ( $name and (! $a_vals ) and ( $input =~ /\A A \s+ \| \s+ \d .* \z/xms ) ) { 
            $a_vals = $input;
            $a_vals =~ s/\|//g;
        }
        elsif ( $name and (! $c_vals ) and ( $input =~ /\A C \s+ \| \s+ \d .* \z/xms ) ) {
            $c_vals = $input;
            $c_vals =~ s/\|//g;
        }
        elsif ( $name and (! $g_vals ) and ( $input =~ /\A G \s+ \| \s+ \d .* \z/xms ) ) {
            $g_vals = $input;
            $g_vals =~ s/\|//g;
        }
        elsif ( $name and (! $t_vals ) and ( $input =~ /\A T \s+ \| \s+ \d .* \z/xms ) ) {
            $t_vals = $input;
            $t_vals =~ s/\|//g;
        }
        # Gene Cluster size = 260683:
        elsif ( $name and (! $cluster_size ) and ( $input =~ /\A Gene [ ] Cluster [ ] size [ ] = [ ] (\d+) : \s* \z/xms ) ) {
            $cluster_size = $1;
        }
    }
    close $INPUT;
    if (    (! $name )
         or (! $consensus )
         or (! $a_vals )
         or (! $c_vals )
         or (! $g_vals )
         or (! $t_vals )
        or (! $cluster_size ) ) {
        die "Failed to extract all needed data from $infile\n";
    }
    my $wb_index = sprintf "%08i", $i;
    $wb_index    = 'WBPmat' . $wb_index;

    print "Position_Matrix : \"$wb_index\"\n";
    print "Brief_id          \"$name", "_Stormo_2012.pfm\"\n";
    print 'Description',
          " \"Conserved regulatory element identified by genome-wide search of",
          " C. elegans, C. briggsae and C. remanei genomes.\"",
          " Paper_evidence \"WBPaper00040751\" \/\/ pmid22540038 \n",
          ;
    print "Type           Frequency\n";
    print "Consensus      \"$consensus\"\n";
    print "Site_values    $a_vals\n";
    print "Site_values    $c_vals\n";
    print "Site_values    $g_vals\n";
    print "Site_values    $t_vals\n";
    print "Sites_used     $cluster_size\n";
    print "Remark \"From a genome-wide comparative analysis in WBPaper00040751\/pmid22540038; its C. elegans genome sequence was from WS170;",
          " its C. briggsae genome sequence was from Stein et al., 2003 (WBPaper00006396\/pmid14624247);",
          " and its C. remanei genome sequence was from genome.wustl.edu.\"\n",
          ;
    print "\n";

    # After printing out this file's matrix data, clear all the data variables to zero so that the next file's data can be read.
    $name         = q{};
    $consensus    = q{};
    $a_vals       = q{};
    $c_vals       = q{};
    $g_vals       = q{};
    $t_vals       = q{};
    $cluster_size = q{};

    # Count one more matrix as having been read and printed.
    $i++;
}

