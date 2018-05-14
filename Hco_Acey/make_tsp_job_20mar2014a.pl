#!/usr/bin/env perl

use strict;
use warnings;

my $header = '#!/bin/bash' . "\n\n";

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A more_seqs\/(\S+\.pep)\.all\.fa \z/xms ) {
        my $stem   = $1;
        my $prefix = q{};
        if ( $stem =~ /\A ([A-Z])[a-z]* [_] ([a-z]{3})[a-z]* /xms ) { 
            my $prefix1 = $1;
            my $prefix2 = $2;
            $prefix = $prefix1 . $prefix2 . q{_};
        }
        else { 
            die "Can't parse a prefix from stem $stem\n";
        }
        print $header if $header;
        $header = q{};
        print 'get_largest_isoforms.pl -t ens -i ',
              "$input | prettify_ensembl_proteome_headers_18apr2014.pl $prefix - > nice_seqs/",
              "$stem.max_isos.pretty.fa ;\n\n",
              ;
    }
    else { 
        die "Can't parse input: $input\n";
    }
}

# Print empty footer line, if this worked:
if (! $header) { 
    print "\n";
}

