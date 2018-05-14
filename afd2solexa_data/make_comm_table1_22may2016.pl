#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "Comparison\tDESeq2_genes\tedgeR_genes\tOverlap_genes\tOverlap_ratio";

my %file2count = ();

while (my $input = <>) {
    chomp $input;
    $input =~ s/\A\s+//;
    $input =~ s/\s+\z//;
    if ( $input =~ /\A (\d+) \s+ (\S+)/xms ) { 
        my $gene_count = $1;
        my $infile     = $2;

        if (! -r $infile) {
            die "Can't read putative infile in: $input\n";
        }
        if ( exists $file2count{$infile} ) {
            die "Redundant gene counts for $infile: $file2count{$infile} vs. $gene_count\n";
        }

        $file2count{$infile} = $gene_count;

        # Slightly redundant, but helps keep data straight to rename it:
        my $edgeR_file  = $infile;
        my $DESeq2_file = q{};

        if ( $edgeR_file =~ /\A edgeR_up_down_genelists\/ (\S+) \.edgeR\.gene_list\.txt \z/xms ) { 
            my $comp_type = $1;

            # Identify correct matching file, verify its existence:
            $DESeq2_file = q{DESeq2_up_down_genelists/} . $comp_type . q{.deseq2.gene_list.txt};
            if (! -r $DESeq2_file ) {
                die "Can't correctly predict readable DESeq2 file ($DESeq2_file) to compare to $edgeR_file\n";
            }

            # Confirm that every single DESeq2 gene count is >= its matching edgeR count.
            if (! exists $file2count{$DESeq2_file} ) {
                die "Failed to count genes in $DESeq2_file\n";
            }
            if ( $file2count{$DESeq2_file} < $file2count{$edgeR_file} ) {
                die "No!  It cannot be!  There are more genes in $edgeR_file ($file2count{$edgeR_file}) than in $DESeq2_file ($file2count{$DESeq2_file})!\n";
            }

            # Get overlap count from system:
            my $overlap = `   comm -12 $edgeR_file $DESeq2_file | wc -l ;`;
            chomp $overlap;
            # Deal with zero-line values:
            if ( $overlap !~ /\S/xms ) { 
                $overlap = 0;
            }

            # Get ratio (when not dividing by zero, in which case it defaults to 'n/a'):
            my $ratio   = 'n/a';
            if ( $file2count{$edgeR_file} >= 1 ) {
                $ratio = ($overlap/$file2count{$edgeR_file});
                $ratio = sprintf "%.3f", $ratio;
                $ratio = "$ratio [$overlap/$file2count{$edgeR_file}]";
            }

            print "$header\n" if $header;
            $header = q{};

            my $DESeq2_count = $file2count{$DESeq2_file};
            $DESeq2_count    = commify($DESeq2_count);
            my $edgeR_count  = $file2count{$edgeR_file};
            $edgeR_count     = commify($edgeR_count);
            $overlap         = commify($overlap);

            print "$comp_type\t$DESeq2_count\t$edgeR_count\t$overlap\t$ratio\n";
        }
    }
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

