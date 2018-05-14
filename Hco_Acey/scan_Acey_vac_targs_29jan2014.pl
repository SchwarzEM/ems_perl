#!/usr/bin/env perl

use strict;
use warnings;
use List::MoreUtils qw(uniq);

my %wanted = ();

my @files = @ARGV;

my $site_list = shift @files;

my $data_ref;

open my $SITES, '<', $site_list or die "Can't open wanted site list $site_list: $!";
while (my $input = <$SITES>) { 
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) { 
        $wanted{$input} = 1;
    }
    else { 
        die "From wanted site list $site_list, can't parse: $input\n";
    }
}
close $SITES or die "Can't close filehandle to wanted site list $site_list: $!";

foreach my $infile (@files) {
    my $target = q{};
    open my $INFILE, '<', $infile or die "Can't open target non-site list $infile: $!";
    while (my $input = <$INFILE>) {
        chomp $input;
        if ( $input =~ /\A The [ ] following [ ] sites [ ] are [ ] absent [ ] in [ ] : [ ] (\S+) [ ]+ \-\- /xms ) { 
            $target = $1;
        }
        else {
            $input =~ s/\A\s+//;
            $input =~ s/\s+\z//;
            my @sites = grep { $wanted{$_} } split /\s+/, $input;
            push @{ $data_ref->{'target'}->{$target}->{'sites'} }, @sites;
        }
    }
    close $INFILE or die "Can't close filehandle to target non-site list $infile: $!";
}

my @targets = sort keys %{ $data_ref->{'target'} };
foreach my $target (@targets) { 
    my @sites = @{ $data_ref->{'target'}->{$target}->{'sites'} };
    @sites = sort @sites;
    @sites = uniq @sites;
    my $site_text = join '; ', @sites;
    print "$target\t$site_text\n";
}

