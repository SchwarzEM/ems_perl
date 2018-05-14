#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

my $type      = q{};
my $number    = q{};
my $cell      = q{};
my $readcount = q{};

my $raw_reads_file    = $ARGV[0];
my $mapped_reads_file = $ARGV[1];
my $mapped_genes_file = $ARGV[2];

open my $RAW, '<', $raw_reads_file or die "Can't open raw reads file $raw_reads_file: $!";
while ( my $input = <$RAW> ) {
    chomp $input;
    if ( $input =~ / fastq [ ] count [ ] of [ ] \d+_allreads_(\S+)\.LC_no(\d+)_retro_38nt_14may2012.fq /xms ) { 
        $type      = $1;
        $number    = $2;
        $cell      = $type . q{.} . $number;
        $readcount = 0;
    }
    elsif ( $input =~ / Total: [ ]+ \S+ [ ] nt [ ] in [ ] (\S+) [ ] reads\. /xms ) { 
        $readcount = $1;
        $readcount =~ s/,//g;
        $data_ref->{'cell'}->{$cell}->{'raw_readcount'} = $readcount;
    }
    elsif ( $input =~ /\S/xms ) { 
        die "Can't parse input from raw reads file $raw_reads_file: $input\n";
    }
}
close $RAW or die "Can't close filehandle to raw reads file $raw_reads_file: $!";

open my $MAPPED, '<', $mapped_reads_file or die "Can't open mapped reads file $mapped_reads_file: $!";
while ( my $input = <$MAPPED> ) {
    chomp $input;
    if ( $input =~ / Count [ ] of [ ] mapped [ ] reads [ ] in [ ] \d+_(\S+)\.LC_no(\d)_retro_38nt_14may2012\.ce6spl38spk\.bowtie\.txt /xms ) { 
        $type      = $1;
        $number    = $2;
        $cell      = $type . q{.} . $number;
        $readcount = 0;
    }
    elsif ( $input =~ /\A (\d+) \s* \z/xms ) { 
        $readcount = $1;
        $data_ref->{'cell'}->{$cell}->{'mapped_readcount'} = $readcount;
    }
    elsif ( $input =~ /\S/xms ) {
        die "Can't parse input from mapped reads file $mapped_reads_file: $input\n";
    }
}
close $MAPPED or die "Can't close filehandle to raw reads file $mapped_reads_file: $!";

open my $GENES, '<', $mapped_genes_file or die "Can't open mapped genes file $mapped_genes_file: $!";
while ( my $input = <$GENES> ) { 
    chomp $input;
    if ( $input =~ / Count [ ] of [ ] mapped [ ] genes [ ] in [ ] \d+_(\S+)\.LC_no(\d)_retro_38nt_14may2012\.WS220\.rpkm /xms ) {
        $type      = $1;
        $number    = $2;
        $cell      = $type . q{.} . $number;
        $readcount = 0;
    }
    elsif ( $input =~ /\A (\d+) \s* \z/xms ) {
        $readcount = $1;
        $data_ref->{'cell'}->{$cell}->{'mapped_genes'} = $readcount;
    }
    elsif ( $input =~ /\S/xms ) {
        die "Can't parse input from mapped genes file $mapped_genes_file: $input\n";
    }
}
close $GENES or die "Can't close filehandle to raw reads file $mapped_genes_file: $!";

my @types = qw ( L3 L4 nhr );
my @nums = (1..5);

print "\n";

foreach my $type1 (@types) { 
    foreach my $num1 (@nums) { 
        my $cell1    = $type1 . q{.} . $num1;
        my $cell1_pr = $type1 . q{.indiv} . $num1;

        my $raw_reads    = $data_ref->{'cell'}->{$cell1}->{'raw_readcount'};
        my $mapped_reads = $data_ref->{'cell'}->{$cell1}->{'mapped_readcount'};
        my $mapped_genes = $data_ref->{'cell'}->{$cell1}->{'mapped_genes'};

        my $raw_reads_pr    = commify($raw_reads);
        my $mapped_reads_pr = commify($mapped_reads);
        my $mapped_genes_pr = commify($mapped_genes);

        my $percent_mapped = ( ($mapped_reads * 100) / $raw_reads );
        $percent_mapped    = sprintf("%.2f", $percent_mapped);

        my $percent_genes  = ( ($mapped_genes * 100) / 20_252 );
        $percent_genes     = sprintf("%.2f", $percent_genes);

        print $cell1_pr, "\t", $raw_reads_pr, "\t", "$mapped_reads_pr (", $percent_mapped, '%)', "\t", "$mapped_genes_pr (", $percent_genes, '%)', "\n"; 
    }
    print "\n";
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

