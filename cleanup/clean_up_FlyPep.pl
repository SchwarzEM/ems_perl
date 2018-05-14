#!/usr/bin/perl

# clean_up_FlyPep.pl:  Erich Schwarz <emsch@its.caltech.edu>, 4/26/05.
# Purpose: give FlyPep concise headers, while splitting off hairy FlyPep headers to separate file.

$headerfile = $ARGV[0] . ".full_headers";
open (HEADERS, ">$headerfile") || die;

while (<>)
{
    $input = $_;
    chomp ($input);
    if ($input =~ /^>(\S+)\s+.+\[(.+ symbol\S+)\s+.*\s+(FB\D+\d+)\D+/) { 
        while ($store) {
            $store =~ m/^(.{0,60})(.*)/;
            print "$1\n";
            $store = $2;
        }
        print ">$1    $2    $3\n";
        print HEADERS "$input\n\n";
    }
    elsif ($input =~ /^>(\S+)\s+.+\[(.+ symbol\S+)\s+(FB\D+\d+)\D+/) { 
        while ($store) {
            $store =~ m/^(.{0,60})(.*)/;
            print "$1\n";
            $store = $2;
        }
        print ">$1    $2    $3\n";
        print HEADERS "$input\n\n";
    }
    elsif ($input =~ /^>(\S+)\s+.+\[(.+ symbol\S+)\s+/) { 
        while ($store) {
            $store =~ m/^(.{0,60})(.*)/;
            print "$1\n";
            $store = $2;
        }
        print ">$1    $2\n";
        print HEADERS "$input\n\n";
    }
    elsif ($input =~ /^>(\S+)\s+/) { 
        while ($store) {
            $store =~ m/^(.{0,60})(.*)/;
            print "$1\n";
            $store = $2;
        }
        print ">$1\n";
        print HEADERS "$input\n\n";
    }
    unless ($input =~ /^>/) {
        $input =~ s/[^a-zA-Z]//g;
        $store .= $input;
    }
}
while ($store) {
    $store =~ m/^(.{0,60})(.*)/;
    print "$1\n";
    $store = $2;
}

close HEADERS;

