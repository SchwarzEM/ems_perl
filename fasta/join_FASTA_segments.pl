#!/usr/bin/perl

# basic_fasta_clean.pl, Erich Schwarz <emsch@its.caltech.edu>, 8/4/05.
# Purpose: catenate alignment regions in FASTA format with COMMON_NAME/x-y names (produced by JalView), keeping input order.

$i = 0;

while (<>) { 
    chomp ($input = $_);
    if ($input =~ /^>((.+?)[\/\s]+.+)/) { 
        $name = $2;
        $info{$name} .= "[" . $1 . "] ";
        unless ($ident{$name}) { $ident{$name} = ++$i; }
    }
    unless ($input =~ /^>/) { 
        $store{$name} .= $input;
    }
}
foreach $name (sort {$ident{$a} <=> $ident{$b}} keys %ident) { 
    print ">$name $info{$name}\n";
    while ($store{$name}) { 
        $store{$name} =~ m/^(.{0,60})(.*)/;
        print "$1\n";
        $store{$name} = $2;
    }
}
