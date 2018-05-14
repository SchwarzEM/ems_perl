#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $output_file = 'repbase_spp_dnas.fa';
$output_file    = safename($output_file);

while (my $input = <>) { 
    chomp $input;
    # repbase/RepBase19.02.embl/angrep.ref
    if ( $input =~ /\A repbase\/RepBase19\.02\.embl\/([a-z]+)rep\.ref \z/xms ) { 
        my $prefix = $1;
        my $target_file = 'repbase/RepBase19.02.fasta/' . $prefix . 'rep.ref';
        if ( (! -e $input ) or (! -e $target_file ) ) { 
            die "Can't find one of these two files:\n$input\n$target_file\n";
        }
        my $grabbed_species_info = qx/grep 'OS ' $input | head --lines=1;/;
        chomp $grabbed_species_info;
        if ( $grabbed_species_info =~ /\A OS \s+ ([A-Z]) [a-z]+ \s+ ([a-z]{4}) [a-z]* \z/xms ) {
            my $prefix_start = $1;
            my $prefix_end   = $2;
            my $prefix       = $prefix_start . $prefix_end . '.RB_';
            system "tag_FASTA_names.pl -i $target_file -p $prefix >> $output_file ;";
        }
        elsif ( $grabbed_species_info =~ /\A OS \s+ ([A-Z][a-z]{4}) [a-z]* \z/xms ) {
            my $prefix = $1;
            $prefix    = $prefix . '.RB_';
            system "tag_FASTA_names.pl -i $target_file -p $prefix >> $output_file ;";
        }
        # For dealing with really wacky 'species' like "Mammalian expression vector pCMV-Script"
        elsif ( $grabbed_species_info =~ /\A OS \s+ .+ \s+ (p[A-Z]+\S+) \z/xms ) {
            my $prefix = $1;
            $prefix    = $prefix . '.RB_';
            system "tag_FASTA_names.pl -i $target_file -p $prefix >> $output_file ;";
        }
        # For "Murine leukemia virus":
        elsif ( $grabbed_species_info =~ /\A OS \s+ ([A-Z]) [a-z]+ \s+ ([a-z]{3}) [a-z]* \s+ ([a-z]{3}) [a-z]* \z/xms ) {
            my $prefix_start  = $1;
            my $prefix_middle = $2;
            my $prefix_end    = $3;
            my $prefix        = $prefix_start . $prefix_middle . $prefix_end . '.RB_';
            system "tag_FASTA_names.pl -i $target_file -p $prefix >> $output_file ;";
        }
        else { 
            die "Can't parse captured species info: $grabbed_species_info\n";
        }
    }
    else { 
        die "Can't parse input: $input\n";
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


