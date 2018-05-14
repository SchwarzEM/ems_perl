#!/usr/bin/env perl

# plot_paired_blastns_24aug2011.pl -- Erich Schwarz <emsch@caltech.edu>, 8/24/2011.
# Purpose: given two vintage-August BlastN reports in output format '6', extract % identities for the reports for any sequences having hits in both.

use strict;
use warnings;
use Getopt::Long;

my $sp1    = q{};
my $sp2    = q{};
my $label1 = q{};
my $label2 = q{};

my $data_ref;
my $help;

GetOptions ( 'sp1=s'    => \$sp1,
             'sp2=s'    => \$sp2,
             'label1=s' => \$label1,
             'label2=s' => \$label2,
             'help'     => \$help,   );

if ( $label1 !~ /\S/xms ) { 
    $label1 = 'species_1';
}
if ( $label2 !~ /\S/xms ) {  
    $label2 = 'species_2';
}

if ($label1 eq $label2) {
    die "Need different labels (not $label1 and $label2).\n";
}

if ($help or (! $sp1) or (! $sp2)) { 
    die "Format: plot_paired_blastns_24aug2011.pl",
        " --sp1 [BlastN for sp. 1] --sp2 [BlastN for sp. 2]",
        " --label1 [optional label for sp. 1] --label2 [optional label for sp. 2]",
        " --help",
        "\n",
        ;
}

open my $SP1, '<', $sp1 or die "Can't open BlastN report $sp1 for $label1: $!";
while (my $input = <$SP1>) { 
    chomp $input;
    my $seq = q{};
    my $id  = q{};
    if ( $input =~ /\A ([^\t]+) \t [^\t]* \t (\d{2,3}\.\d{2}) \t /xms ) { 
        $seq = $1;
        $id  = $2;
        $data_ref->{'seq'}->{$seq}->{'species'}->{$label1} = $id;
    }
    else { 
        die "Can't parse input from $sp1: $input\n";
    }
}
close $SP1 or die "Can't close filehandle to BlastN report $sp1 for $label1: $!";

open my $SP2, '<', $sp2 or die "Can't open BlastN report $sp2 for $label2: $!";
while (my $input = <$SP2>) {
    chomp $input;
    my $seq = q{};
    my $id  = q{};
    if ( $input =~ /\A ([^\t]+) \t [^\t]* \t (\d{2,3}\.\d{2}) \t /xms ) {
        $seq = $1;
        $id  = $2;
        $data_ref->{'seq'}->{$seq}->{'species'}->{$label2} = $id;
    }
    else {
        die "Can't parse input from $sp2: $input\n";
    }
}
close $SP2 or die "Can't close filehandle to BlastN report $sp2 for $label1: $!";

print 'Gene', "\t", $label1, "\t", $label2, "\n", ;

foreach my $sequence (sort keys %{ $data_ref->{'seq'} } ) { 
    if ( ( exists $data_ref->{'seq'}->{$sequence}->{'species'}->{$label1} ) and ( exists $data_ref->{'seq'}->{$sequence}->{'species'}->{$label2} ) ) { 
        print $sequence,
              "\t",
              $data_ref->{'seq'}->{$sequence}->{'species'}->{$label1},
              "\t",
              $data_ref->{'seq'}->{$sequence}->{'species'}->{$label2},
              "\n",
              ;        
    }
}

