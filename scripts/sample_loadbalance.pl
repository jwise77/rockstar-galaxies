#!/usr/bin/perl -w

#Note for non-Perl users: "scalar <>" returns a single line from STDIN.
#This script implements an equal-volume load balancer for Rockstar
my ($num_writers) = split(/ /, scalar <>);
print STDERR "[Loadbalance] Num writers: $num_writers\n";
my (@chunks) = split(/ /, scalar <>);
print STDERR "[Loadbalance] Recommended chunks: @chunks[0..2]\n";
my ($box_size) = split(/ /, scalar <>);
print STDERR "[Loadbalance] Box size: $box_size\n";
my ($scale) = split(/ /, scalar <>);
print STDERR "[Loadbalance] Scale factor: $scale\n";

scalar <>; #Skip format specification

while (<>) {
    last if (/^#/); #Signifies return format specification
    my ($id, $ip, $port) = split;
    print STDERR "[Loadbalance input] $_";
    my $output = "$id $ip $port ";
    my $x_min = ($id % $chunks[0])/($chunks[0]) * $box_size;
    my $x_max = (($id % $chunks[0])+1)/($chunks[0]) * $box_size;
    $id = int($id/$chunks[0]);
    my $y_min = ($id % $chunks[1])/($chunks[1]) * $box_size;
    my $y_max = (($id % $chunks[1])+1)/($chunks[1]) * $box_size;
    $id = int($id/$chunks[1]);
    my $z_min = ($id / $chunks[2]) * $box_size;
    my $z_max = (($id+1) / $chunks[2]) * $box_size;
    $output .= " $x_min $y_min $z_min $x_max $y_max $z_max\n";
    push @output, $output;
    print STDERR "[Loadbalance output] $output";
}

#Note: need to print output *after* all input has been received.
print @output;
