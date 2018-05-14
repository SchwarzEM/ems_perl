#!/usr/bin/perl

while (<>) {
    chomp($input = $_);
    @numbers = split /,/, $input;
    print "$numbers[2]\n";
}
