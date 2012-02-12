#!/usr/bin/env perl

use strict;
use warnings;

# print "Type & Size & Time\n";

undef $/;
my $line = <>;
my $last = "";
# while ($line =~ m{^simulation1\(\) using payload size\s+(\d+)\s+and conversion code from (\w+)\s*\n.*\n^([0-9]+\.?[0-9]+) msec.+$}mg) {
while ($line =~ m{^simulated\w+\(\) using payload size\s+(\d+)\s+and conversion code from (\w+)\s*\n^([0-9]+\.?[0-9]+) msec.+$}mg) {
  if (($2 ne $last) and ($last ne "")) {
    $last = $2;
    printf("\n\n");
  }
   printf("\"%s\" %d %g \n", $2, $1, $3);
  $last = $2;
#   printf("%s %d\n", $2, $1);
};

