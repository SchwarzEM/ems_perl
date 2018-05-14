#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $file = q{};
my $md5  = q{};

my @files    = ();
my %file2md5 = ();

# Sample input --
# Sample: jj142_GACTCA
# File: 7493_1506_41050_HCHK5BGXY_jj142_GACTCA_R1.fastq.gz
# Size 892892924 bytes, MD5: 20de0d8015fd7ac67f7899fd4ebafb92
# Link: http://cbsuapps.tc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=1754321425&refid=167231

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A File: [ ] (\S+) \s* /xms ) { 
        $file = $1;
    }
    elsif ( ( $file =~ /\S/xms ) and ( $input =~ /\A Size [ ] \d+ [ ] bytes, [ ] MD5: [ ] (\S+) \s* /xms ) ) {     
        $md5 = $1;
        $file2md5{$file} = $md5;
        push @files, $file;
        $file = q{};
        $md5  = q{};
    }
}

@files = sort @files;
foreach my $file_final (@files) {
    my $md5_final = $file2md5{$file_final};
    print "$md5_final  $file_final\n";
}

