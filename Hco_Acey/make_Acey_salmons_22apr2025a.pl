#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $read_1 = <> ) {
    chomp $read_1;
    my $read_2 = $read_1;
    $read_2 =~ s/_1\.fq\.gz\z/_2.fq.gz/;

    if ( $read_1 eq $read_2 ) {
        die "Failed to readpair $read_1\n";
    }
    if (! -e $read_2 ) {
        die "Read_2 does not exist: $read_2\n";
    }

    my $tag = q{};
    if ( $read_1 =~ / \S+ \/ (\S+?_rep\d)_\d\.fq.gz \z/xms ) {
        $tag = $1;
    }
    else {
        die "Cannot extract tag from $read_1\n";
    }

    print 'salmon --no-version-check quant --threads 32 --seqBias --gcBias --posBias --discardOrphansQuasi';
    print ' --index $PROJECT/Acey/2025.04.01/dbs/Acey.v2.1_WBPS19_gentrome_index';
    print ' --libType A';
    print " --mates1 $read_1";
    print " --mates2 $read_2";
    print " --output $tag.salmon";
    print ' --geneMap $PROJECT/Acey/2025.04.01/annots/Acey_v2.1.2022.11.14.01.cds2gene.tsv.txt ;';
    print "\n";
}

