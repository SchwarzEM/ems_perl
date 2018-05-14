#!/usr/bin/env perl

# remove_aug_genes.pl -- Erich Schwarz <emsch@caltech.edu>, 10/26/2012.
# Purpose: given a list of genes to be censored and an *.aug output from AUGUSTUS 2.6.1, print an *.aug without the unwanted genes.

use strict;
use warnings;
use Getopt::Long;

my $data_ref;

my $print    = 1;
my $censored = q{};
my $augustus = q{};
my $help;

GetOptions ( 'censored=s' => \$censored,
             'augustus=s' => \$augustus,
             'help'       => \$help, );

if ( (! $censored) or (! $augustus) or $help ) { 
    die "Format: remove_aug_genes.pl\n",
        "    --censored|-c  [list of genes to censor]\n",
        "    --augustus|-a  [AUGUSTUS .aug file to censor]\n",
        "    --help|-h      [print this message]\n",
        ;
}

open my $CENSOR, '<', $censored or die "Can't open file of genes to be censored: $censored\n";
while (my $input = <$CENSOR>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) { 
        my $unwanted = $1;
        $data_ref->{'unwanted'}->{$unwanted} = 1;
    }
    else { 
        die "From censored file $censored, can't parse: $input\n";
    }
}
close $CENSOR or die "Can't close filehandle to file of genes to be censored: $censored\n";

# Sample input:
# end gene Acey_2012.08.05_0001.g10
###
# start gene Acey_2012.08.05_0001.g11

open my $AUG, '<', $augustus or die "Can't open AUGUSTUS file to be censored: $augustus\n";
while (my $input = <$AUG>) {
    chomp $input;
    if ( $input =~ /\A [#] [ ] start [ ] gene [ ] (\S+) /xms ) {
        my $gene = $1;
        if ( exists $data_ref->{'unwanted'}->{$gene} ) { 
            $print = 0;
        }
    }
    elsif ( $input =~ /\A [#] [ ] end [ ] gene [ ] (\S+) /xms ) {
        my $gene = $1;
        if ( exists $data_ref->{'unwanted'}->{$gene} ) {
            $print = -1;   
        }
    }
    if ( $print == 1 ) { 
        print "$input\n";
    }
    elsif ( $print == -1 ) { 
        $print = 1;
    }
}
close $AUG or die "Can't close filehandle to AUGUSTUS file to be censored: $augustus\n";

