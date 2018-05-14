#!/usr/bin/env perl

# URLgff3_2_genecount.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/11/2008.
# Purpose: given list of one or more URLs of gff3 files, download and get gene counts, e.g., for nGASP.

use strict;
use warnings;

my $url_list = $ARGV[0];
my @urls     = ();

if ( (! $ARGV[0]) or (! -r $url_list ) ) { 
    die "Format: ./gff3_2_genecount.pl [URL list file]\n";
}

open my $URL_LIST, '<', $url_list 
    or die "Can't open URL list $url_list: $!";

while (my $input = <$URL_LIST> ) { 
    chomp $input;
    $input =~ s/\A\s*//;
    $input =~ s/\s*\z//;
    if ( $input =~ /\A\S+\z/xms ) {
        push @urls, $input;
    }
}

foreach my $url (@urls) { 
    my $gene       = q{}; 
    my %genes_seen = ();

    system "ncftpget -V $url";

    my $gz_filename = $url;
    $gz_filename =~ s{ \A .+ /}{}xms;

    my $filename = $gz_filename;
    $filename =~ s/\.gz\z//;

    if ($gz_filename =~ /\.gz\z/xms) { 
        system "gunzip $gz_filename";
    }

    open my $FILE, '<', $filename 
        or die "Can't open file $filename: $!";

    while ( my $input = <$FILE> ) { 
        if ( $input =~ / \A 
                         .+
                         \t gene \t 
                         .+ \t
                         (?: \S ;)* 
                         ID= 
                         ([^;]+)
                         /xms ) { 
            $gene = $1;
            $genes_seen{$gene} = 1;
        }
    }

    close $FILE or die "Can't close handle of file $filename: $!";

    system "rm $filename";

    my $gene_count = commify(scalar(keys %genes_seen));
    print "$gene_count genes in $filename\n";
}

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

