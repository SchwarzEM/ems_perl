#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $data_ref;

my $header = "library(stats)";

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A ([^\t]+)      # Data_set
                       \t ([^\t]+)   # Motif
                       \t (\S+)  # Orig_search_genes
                       \t \S+    # Orig_positive_genes 
                       \t (\S+)  # Redisc_positive_genes
                       \t \S+    # All_genes
                       \t (\S+)  # All_sepal_genes
                       \t \S+    # Genomewide_hits
                       \t (\S+)  # Sepalwide_hits 
                       \t \S+ \t \S+  # Redisc_orig_query_ratio, Redisc_orig_pos_ratio
                       \t (\S+)  # Gwide_hit_ratio
                       \t \S+ \t \S+ \t \S+ # Sepalwide_hit_ratio     Query_enrichment Sepal_enrichment
                    \z/xms ) { 
        my $data_set        = $1;
        my $motif           = $2;
        my $all_data_genes  = $3;
        my $setwide_hits    = $4;
        my $all_sepal_genes = $5;
        my $sepalwide_hits  = $6;
        my $gwide_hit_ratio = $7;
        if ( $data_set ne 'Data' ) { 
            # Strip out any commas place into the number for readability:
            $all_data_genes  =~ s/[,]//g;
            $setwide_hits    =~ s/[,]//g;
            $all_sepal_genes =~ s/[,]//g;
            $sepalwide_hits  =~ s/[,]//g;

            # Enforce numerical validity:
            if (    (! looks_like_number($all_data_genes)  )
                 or (! looks_like_number($setwide_hits)    )
                 or (! looks_like_number($sepalwide_hits)  ) 
                 or (! looks_like_number($all_sepal_genes) ) 
                 or (! looks_like_number($gwide_hit_ratio) ) ) {
                die "Failed to detect valid numbers in: $input\n";
            }

            # Print the header only once, at the start:
            print "$header\n" if $header;
            $header = q{};

            # Print one stanza of lines for each relevant data set and genomewide comparison:
            print "# Data set: \"$data_set\"; motif \"$motif\"; p-value for freq. in all_sepal_genes\n";
            print 'binom.test('
                  . $sepalwide_hits
                  . q{,}
                  . $all_sepal_genes
                  . ',p='
                  . $gwide_hit_ratio
                  . ',alternative=c("two.sided"),conf.level=0.99)'
                  . "\n"
                  ;

            print "# Data set: \"$data_set\"; motif \"$motif\"; p-value for freq. in all_query_genes\n";
            print 'binom.test('
                  . $setwide_hits 
                  . q{,}
                  . $all_data_genes
                  . ',p='
                  . $gwide_hit_ratio
                  . ',alternative=c("two.sided"),conf.level=0.99)'
                  . "\n"
                  ;
        }
    }
    else { 
        die "Cannot parse input line: $input\n";
    }
}

# Print the final line only once, *if* we have actually printed data:
print "q()\n" if (! $header);

