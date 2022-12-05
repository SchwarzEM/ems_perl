#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    # cut -f 1,22-24,32,34-36
    if ( $input =~ /\A ([^\t]*) \t 
                       (?: [^\t]* \t){20} 
                       ([^\t]*) \t 
                       ([^\t]*) \t 
                       ([^\t]*) \t 
                       (?: [^\t]* \t){7} 
                       ([^\t]*) \t 
                       [^\t]* \t 
                       ([^\t]*) \t 
                       ([^\t]*) \t 
                       ([^\t]*) \t/xms ) {
        print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\n";
    }
    else {
        die "From input, cannot parse: $input\n";
    }
}

