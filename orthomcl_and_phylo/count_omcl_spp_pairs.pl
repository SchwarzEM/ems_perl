#!/usr/bin/env perl

# count_omcl_spp_pairs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/4/2010.
# Purpose: given any OrthoMCL file/stream, count and list all the pairwise species co-occurrences.

# (Very) basic algorithm:
# 1.  Read a line.
# 2.  Get the nonredundant species list for that line.
# 3.  For *each* species, ++ 1 for $omcl_data_ref->{$sp1}->{$sp2}, skipping self-cooccurrences.
# 4.  Having gone through the entire file: 
# 4a.     Start reading off the hashref keys, and for each one, get a count of pairings.
# 4b.     Print each count: species 1; species 2; number of co-occurrences.
# 4b.     As you do this, for each sp1-sp2 pair, *delete* the opposite-but-equal sp2-sp1 pair!
# 4c.     Sooner or later there are no live pairs left.  Stop.

use strict;
use warnings;

my $omcl_data_ref;

my $data_text = q{};
my $species   = q{};
my %obs_spp   = ();

while (my $input = <>) { 
    chomp $input;
    $data_text = q{};
    $species   = q{};
    %obs_spp   = ();

    # Sample input line:
    # ORTHOMCL0(467 genes,6 taxa):      Bm1_00320(b_malayi) [etc.]

    if ( $input =~ / \A 
                     ORTHOMCL\d+ 
                     \( 
                     \d+ [ ] genes,   # no space between this ',' and the next "\d+"!
                     \d+ [ ] taxa \): 

                     # For some absurd reason, Perl can't just work here with
                     #     a single parenthesis pair -- this double-set is needed:

                     ( (?: \s+ [^\s\(\)]+ [(] [^\W\(\)]+ [)] ){2,} )

                     \z /xms ) { 
        $data_text = $1;
        while ( $data_text =~ / \( ( [^\W\(\)]+ ) \) /xmsg ) { 
            $species = $1;
            $obs_spp{$species} = 1;
        }
        foreach my $sp1 (sort keys %obs_spp) { 
            foreach my $sp2 ( grep { $_ ne $sp1 } sort keys %obs_spp ) { 
                $omcl_data_ref->{$sp1}->{$sp2}++;
            }
        }
    }

    # Enforce successful parsing.
    else { 
        die "Can't parse input line: $input\n";
    }
}

my $length = 0;
foreach my $sp ( sort keys %{ $omcl_data_ref } ) { 
    my $len1 = length($sp);
    if ( $len1 > $length ) { 
        $length = $len1;
    }
}

$length = $length * 2;
$length = $length + 5;

foreach my $sp1 ( sort keys %{ $omcl_data_ref } ) { 
    foreach my $sp2 ( sort keys %{  $omcl_data_ref->{$sp1} } ) { 
        my $format = '%' . $length . 's:%7s';
        my $paircount = $omcl_data_ref->{$sp1}->{$sp2};
        $paircount = commify($paircount);
        my $output = sprintf( "$format",
                              "$sp1 + $sp2",
                               $paircount, 
                            );
        print "$output\n";
        delete $omcl_data_ref->{$sp2}->{$sp1};
    }
}

sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}


