#!/usr/bin/perl -w

# clean_up_IPI_Pep.pl:  Erich Schwarz <emsch@its.caltech.edu>, 4/26/05.
# Purpose: Reformat IPI HumanPep and MousePep headers so that simple protein names are clearly used.

# E.g., IPI HumanPep ==  ftp://ftp.ebi.ac.uk/pub/databases/IPI/current/ipi.HUMAN.fasta.gz

while (<>) {
    chomp ($input = $_);
    if ($input =~ />(IPI:([^|]+)\|.*)/) { 
        $header = $1;
        $name   = $2;
        $header =~ s/\|/  /g;
        $header =~ s/;/; /g;
        while ($store) {
            $store =~ m/^(.{0,60})(.*)/;
            print "$1\n";
            $store = $2;
        }
        print ">$name    $header\n";
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
