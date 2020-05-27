#!/usr/bin/env perl

# make_wget_reads_26may2020.pl -- Erich Schwarz <ems394@cornell.edu>, 5/26/2020.

use strict;
use warnings;
use autodie;
use Getopt::Long;
use Cwd;

my $infile   = q{};
my $username = q{};
my $passtext = q{};
my $password = q{};
my $id_desc  = q{};
my %id2desc  = ();

my $help;

my $curr_dir = getcwd;

GetOptions ( 'infile=s'    => \$infile,
             'username=s'  => \$username,
             'passtext=s'  => \$passtext,
             'desc=s'      => \$id_desc,
             'help'        => \$help,
);

if ( $help or (! $infile ) or (! $username ) or (! $passtext ) or (! $id_desc ) ) {
    die "Format: wget_reads_26may2020.pl\n",
        "            --infile|-i    [text file of target URLs for wget]\n",
        "            --username|-u  [user name for wget]\n",
        "            --passtext|-p  [text file with password/passphrase for wget]\n",
        "            --desc|-d      [ID to description text TSV]\n",
        "            --help|-h      [print this message]\n",
        ;
}

if ( $username !~ /\A \S+ \z/xms ) {
    die "Cannot parse username: $username\n";
}

open my $PASSTEXT, '<', $passtext;
while (my $input = <$PASSTEXT>) {
    chomp $input;
    if ( (! $password) and ( $input =~ /\A\s*(\S.*\S)\s*\z/xms ) ) {
        $password = $1;
    }
}
close $PASSTEXT;

open my $ID_DESC, '<', $id_desc;
while (my $input = <$ID_DESC>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $id   = $1;
        my $desc = $2;
        if (exists $id2desc{$id}) {
            die "Redundant descriptions of $id: $id2desc{$id} and $desc\n";
        }
        $id2desc{$id} = $desc;
    }
}
close $ID_DESC;

open my $INFILE, '<', $infile;
while (my $url = <$INFILE>) { 
    chomp $url;

    # Sample input:
    # https://jumpgate.caltech.edu/runfolders/volvox02/170823_SN787_0718_AHMVH5BCXY/Unaligned/Project_19227_index1/Sample_19227

    if ( $url =~ /\A https [:] \/\/ .+ \/ ( [^\s\/]+ ) \z/xms ) {
        my $id = $1;

        if (! exists $id2desc{$id} ) {
            die "Cannot parse ID \"$id\" in: $url\n";
        }

        my $desc = $id2desc{$id};

        print "mkdir $curr_dir/$desc ;\n";

        print "cd $curr_dir/$desc ;\n";

        print 'wget -q --user=', $username, ' --password=', q{"}, $password, q{"},
              ' --recursive --level=1 --no-parent --no-directories --no-check-certificate --accept .fastq.gz',
              " $url ;\n",
              ;

        print "\n";
    }
    else { 
       	die "Cannot parse URL: $url\n";
    }
}
close $INFILE or die "Can't close filehandle to infile $infile: $!";

