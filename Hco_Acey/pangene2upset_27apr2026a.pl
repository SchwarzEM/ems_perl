#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;

my $pangenes  = q{};
my @requireds = ();
my $data_ref;
my $help;

GetOptions ( 'pangenes=s'    => \$pangenes,
             'required=s{,}' => \@requireds,
             'help'          => \$help,   );

if ( $help or (! $pangenes) ) {
    die "Format: pangene2upset_27apr2026a.pl\n",
        "    --pangenes|-p   [pangenes file from GENESPACE]\n",
        "    --required|-r   [OPTIONAL: one or more required taxa]\n",
        "    --help|-h       [print this message]\n",
        ;
}

open my $PANGENES, '<', $pangenes;
while ( my $input = <$PANGENES> ) {
    chomp $input;

    # Sample inputs:
    # 0_368   1       Necator_chrI    1       0_368   Aroian  35669   PASS    Necator_chrI.g137       Necator_chrI    13892   14485   1
    # 1_8011  1       Necator_chrI    1       0_368   Baylor  14610   NSOrtho Anhui_chrII.g9236       Anhui_chrII     28747090        28747683        7206

    if (     ( $input !~ /\A ofID \t/xms ) 
         and ( $input =~ /\A \S+ \t (\S+) \t (?: \S+ \t){3} (\S+) \t (?: \S+ \t){4} \d+ \t \d+ \t \d+ \z/xms ) 
       ) {
        my $pangene = $1;
        my $taxon   = $2;

        $data_ref->{'taxon'}->{$taxon} = 1;
        $data_ref->{'pangene'}->{$pangene}->{'taxa'}->{$taxon} = 1;
    }
}
close $PANGENES;

my @all_pangenes = sort keys %{ $data_ref->{'pangene'} };
my @pangenes     = ();

if (@requireds) {
    foreach my $pangene (@all_pangenes) {
        my $ok = 1;
        foreach my $required (@requireds) {
            if (! exists $data_ref->{'pangene'}->{$pangene}->{'taxa'}->{$required} ) {
                $ok = 0;
            }
        }
        if ($ok) {
            push @pangenes, $pangene;
        }
    }
}
else {
    @pangenes = @all_pangenes;
}

my @taxa = sort keys %{ $data_ref->{'taxon'} };
my $head = join "\t", @taxa;
$head    = "Gene\t$head";

print "$head\n";

foreach my $pangene (@pangenes) {
    print "$pangene";
    foreach my $taxon (@taxa) {
        my $i = 0;
        if ( exists $data_ref->{'pangene'}->{$pangene}->{'taxa'}->{$taxon} ) {
            $i = 1;
        }
        print "\t$i";
    }
    print "\n";
}

