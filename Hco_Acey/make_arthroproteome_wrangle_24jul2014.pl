#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @seqfiles = ();

my $header = '#!/bin/bash' . "\n\n";

while (my $input = <>) { 
    chomp $input;
    $input =~ s/[#].+\z//;
    my $orig_command = $input;
    $input =~ s/\A\s+//;
    $input =~ s/\s+\z//;
    if ( $input =~ /\A .+ \/ ([^\/\s]+\.fa).gz \s+ [;] \z/xms ) { 
        my $seqfile = $1;
        push @seqfiles, $seqfile;
    }
    else {
        die "Can't parse input line: $input\n";
    }

    # If there is anything to print, then print header, but do it just once.
    print $header if $header;
    $header = q{};

    print "$orig_command\n";
}

print "\n";
print "    gunzip *.fa.gz ;\n";
print "\n";

@seqfiles = sort @seqfiles;

foreach my $seqfile (@seqfiles) { 
    if ( $seqfile =~ /\A ([A-Z][a-z]) [a-z]* _ ([a-z]{3}) [a-z]* /xms ) { 
        my $prefix1 = $1;
        my $prefix2 = $2;
        my $prefix = $prefix1 . $prefix2 . q{_};
        print "    mv -i $seqfile $seqfile.orig ;\n";
        print "    tag_FASTA_names.pl -p $prefix -i $seqfile.orig > $seqfile ;\n";
        print "\n";
    }
    else { 
        die "Can't parse sequence file name: $seqfile\n";
    }
}

