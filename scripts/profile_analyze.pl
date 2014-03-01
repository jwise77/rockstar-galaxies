#!/usr/bin/perl -w

opendir DIR, ".";
my @profiles = grep { /^profile\.\d+$/ } readdir DIR;
closedir DIR;

for (@profiles) {
    my ($num) = /(\d+)/;
    open FILE, "<$_";
    while (<FILE>) {
	my @a = split;
	my ($snap, $chunk) = ($a[1] =~ /S(\d+),C(\d+)/);	
	if (/fof/) {
	    ($foftime{$snap}{$chunk}) = ($a[2] =~ /(\d+)/);
	    ($conftime{$snap}{$chunk}) = ($a[7] =~ /(\d+)/);
	    ($fofs{$snap}{$chunk}) = ($a[3] =~ /(\d+)/);
	    ($parts{$snap}{$chunk}) = ($a[5] =~ /(\d+)/);
	}
	if (/wt/) {
	    ($hp{$snap}{$chunk}) = ($a[2] =~ /(\d+)p/);
	    #($hh{$snap}{$chunk}) = ($a[2] =~ /(\d+)h/);
	    ($hw{$snap}{$chunk}) = ($a[2] =~ /(\d+)w/);
	    ($hpw{$snap}{$chunk}) = ($hw{$snap}{$chunk}) ? (($hp{$snap}{$chunk})  /  ($hw{$snap}{$chunk})) : 0;
	    ($wt{$snap}{$chunk}) = ($a[3] =~ /(\d+)/);
	    ($rcv{$snap}{$chunk},$bp_rcv{$snap}{$chunk}) = ($a[4] =~ /(\d+)s,(\d+)/);
	    ($snd{$snap}{$chunk}) = ($a[5] =~ /(\d+)/);
	    ($wk{$snap}{$chunk}) = ($a[6] =~ /(\d+)/);
	    ($idl{$snap}{$chunk}) = ($a[7] =~ /(\d+)/);
	}
    }
    close FILE;
}

for my $snap (sort {$a <=> $b} keys %wk) {
    print "Snap $snap:\n";
    print "Total Particles: ", total($parts{$snap}), "\n";
    print "Total Fofs: ", total($fofs{$snap}), "\n";
    print "Fof Analysis: ", minmaxavg($foftime{$snap});
    print "Conf. Time:   ", minmaxavg($conftime{$snap});
    print "Total Wrku:   ", total($hw{$snap}), "\n";
    print "Wrku:         ", minmaxavg($hw{$snap}, 1);
    print "Avg. Particles per Wrku:   ", minmaxavg($hpw{$snap}, 1);
    print "Wait Time:     ", minmaxavg($wt{$snap});
    print "Recv Time:     ", minmaxavg($rcv{$snap});
    print "BP Recv Time:  ", minmaxavg($bp_rcv{$snap});
    print "Work Time:     ", minmaxavg($wk{$snap});
    print "Send Time:     ", minmaxavg($snd{$snap});
    print "Idle Time:     ", minmaxavg($idl{$snap});
    print "\n";
}


sub total {
    my $input = shift;
    my $total = 0;
    $total += $_ for (values %$input);
    return $total;
}

sub minmaxavg {
    my $input = shift;
    my $nos = shift;
    my @vals = values %$input;
    my ($min, $max) = ($vals[0], $vals[0]);
    my $total = 0;
    for (@vals) {
	$total += $_;
	$min = $_ if ($min > $_);
	$max = $_ if ($max < $_);
    }
    my $avg = (@vals) ? $total/@vals : 0;
    if ($nos) {
	return sprintf("Avg: %4d  Min: %4d  Max: %4d\n", $total/@vals,
		       $min, $max);
    }
    return sprintf("Avg: %6.1fs  Min: %4ds  Max: %4ds\n", $total/@vals,
		   $min, $max);
 }
