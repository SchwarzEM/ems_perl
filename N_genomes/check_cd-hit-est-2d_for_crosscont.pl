#!/usr/bin/env perl

# check_cd-hit-est-2d_for_crosscont.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/22/2011.
# Purpose: given input cd-hit-est-2d file, check for cross contaminants.  Presumes that the first database compared is itself non-redundant.

use strict;
use warnings;

my $cluster_count    = 0;
my $cont_clust_count = 0;

my $new_clust_seen   = 0;
my $new_cont_seen    = 0;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > Cluster \s+ \d+ /xms ) { 
        $new_clust_seen = 1;
        $new_cont_seen  = 0;
    }
    elsif ( $new_clust_seen and ( $input =~ / [*] \s* \z/xms ) ) {
        $cluster_count++;
        $new_clust_seen = 0;
    }
    elsif ( ( $input =~ / \S \s* \z/xms ) and ( $input !~ / [*] \s* \z/xms ) ) { 
        if ( ! $new_cont_seen ) { 
            $new_cont_seen = 1; 
            $cont_clust_count++;
        }
    }
    else { 
        die "Can't parse input: $input\n";
    }
}

my $print_clust_count = commify($cluster_count);
my $print_cont_count  = commify($cont_clust_count);

print $print_clust_count;
print ' cluster';
print 's' if ( $cluster_count >= 2 );
print '; ';
print $print_cont_count;
print ' cross-contaminated cluster';
print 's' if ( $cont_clust_count >= 2 );
print ".\n";

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

