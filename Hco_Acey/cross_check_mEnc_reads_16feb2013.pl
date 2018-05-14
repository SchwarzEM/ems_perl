#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $file_list  = q{};
my $data_table = q{};

my $data_ref;

my $verbose;
my $help;

GetOptions ('list=s'  => \$file_list,
            'table=s' => \$data_table,
            'verbose' => \$verbose,
            'help'    => \$help, );

if ( $help or (! $file_list) or (! $data_table) ) { 
    die "cross_check_mEnc_reads_16feb2013.pl\n",
        "    --list|-l [list of actual read files]\n",
        "    --table|-t [data table]\n",
        "    --verbose|-v [loudly announce successes]\n",
        "    --help|-h [print this message]\n",
        ;
}

open my $LIST, '<', $file_list or die "Can't open file list $file_list: $!";
while (my $input = <$LIST>) { 
    chomp $input;
    if ( $input =~ / \A \/sternlab \/redivivus \/data02 \/schwarz \/modENCODE \/raw_SRA_data \/
                        (\w+) \. (SRX\d+) \. (SRR\d+) \/ (\w+) \. (total-RNA [:] \S+ \.fastq \. gz) \z /xms ) { 
        my $type  = $1;
        my $srx   = $2;
        my $srr   = $3;
        my $type2 = $4;
        my $file  = $5;
        if ( $type ne $type2 ) { 
            die "From file list $file_list, inconsistent types in: $input\n";
        }
        if ( exists $data_ref->{'file'}->{$file} ) {
            die "From file list $file_list, redundant data entries for file of: $input\n";
        }
        $data_ref->{'file'}->{$file}->{'type'} = $type;
        $data_ref->{'file'}->{$file}->{'srx'}  = $srx;
        $data_ref->{'file'}->{$file}->{'srr'}  = $srr;
    }
    else {
        die "From file list $file_list, can't parse input: $input\n";
    }
}
close $LIST or die "Can't close filehandle to file list $file_list: $!";

open my $TABLE, '<', $data_table or die "Can't open data table $data_table: $!";
while (my $input = <$TABLE>) {
    chomp $input;
    if ( $input =~ /\A (\w+) \t \S+ \t (SRX\d+) \t (SRR\d+) \t 
                       ftp[:] \/\/ data\.modencode\.org \/ C\.elegans \/ mRNA \/ RNA-seq \/ raw-seqfile_fastq \/ (total-RNA [:] \S+ \.fastq \. gz) \z/xms ) {
        my $type  = $1;
        my $srx   = $2;
        my $srr   = $3;
        my $file  = $4;
        if ( exists $data_ref->{'file'}->{$file} ) { 
            if ( $type ne $data_ref->{'file'}->{$file}->{'type'} ) {
                print "ERROR: Data table $data_table lists discordant types ($type vs. $data_ref->{'file'}->{$file}->{'type'}) for file $file\n";
            }
            if ( $srx ne $data_ref->{'file'}->{$file}->{'srx'} ) { 
                print "ERROR: Data table $data_table lists discordant SRX values ($srx vs. $data_ref->{'file'}->{$file}->{'srx'}) for file $file\n";
            }
            if ( $srr ne $data_ref->{'file'}->{$file}->{'srr'} ) {
                print "ERROR: Data table $data_table lists discordant SRR values ($srr vs. $data_ref->{'file'}->{$file}->{'srr'}) for file $file\n";
            }
            if (     ( $type eq $data_ref->{'file'}->{$file}->{'type'} )
                 and ( $srx  eq $data_ref->{'file'}->{$file}->{'srx'}   )
                 and ( $srr  eq $data_ref->{'file'}->{$file}->{'srr'}   )  ) {
                print "Perfect table/file match for $input!\n" if $verbose;
                delete $data_ref->{'file'}->{$file};
            }
        }
        else { 
            print "ERROR: Data table $data_table lists a file not actually present: $file\n";
        }
    }
    else { 
        die "From data table $data_table, can't parse input: $input\n";
    }
}
close $TABLE or die "Can't close filehandle to data table $data_table: $!";

my @unaccounted_files = sort keys %{ $data_ref->{'file'} };

foreach my $unaccounted_file (@unaccounted_files) { 
    print "Data table $data_table failed to account for actual file: $unaccounted_file\n";
}

