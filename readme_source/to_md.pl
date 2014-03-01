#!/usr/bin/perl -w
our ($verbatim, $infigure, %labels, %refs, %mdref) = (0,0);
our ($sec, $subsec, $subsubsec, $fig) = (0,0,0,0);
my $filename = $ARGV[0];
my $label_fn = $filename;
$label_fn =~ s/\.tex/.labels/;
if (open INPUT, "<$label_fn") {
    while (<INPUT>) {
	chomp;
	my ($a, $b) = split(/=>/, $_);
	$refs{$a} = $b;
    }
    close INPUT;
}

while (<>) {
    if (/tableofcontents/) { parsetoc($filename); next; }
    s/\$<\$/</g;
    s/\$>\$/>/g;
    if (/\\noindent\{\}/) { print "\n" }
    s/``/"/g;
    s/''/"/g;
    s/\\S/Section /g;
    s/\\\\\[[^\]]+\]/\n/;
    s/\\\\//g;
    s/Fig\.\\ /Fig. /g;
    s/\\item/+ /;
    s/\\includegraphics\[[^\]]+\]{([^\.]+)\.\w+}/![$1](http:\/\/slac.stanford.edu\/~behroozi\/Rockstar\/$1.png "$1")/g;
    s/\\footnote\{(.*)\}/ [Footnote: $1]/g;
    s/\\caption\{(.*)\}/"Figure ".($fig).". $1"/eg;
    s/\\(\w+)\{([^}]*)\}/latex($1,$2)/eg;
    s/\\\w+\[[^\]+]+\]\{[^\}]+\}//g;
    s/\$_\\odot\$/sun/g;
    s/\$//g;
    s/\bM_\{?(\w+)\}?/M$1/g;
    s/\bR_\{?(\w)\}?/R$1/g;
    s/^\s*(The Rockstar Halo Finder)/# $1\n/;
    s/^\s+//;
    $_ = "\n" unless length $_;
    my $spaces = 4*((($verbatim>0)?1:0)+(($subsec>0 and !/ ### /)?1:0) + (($subsubsec>0 and !/ #### /)?1:0));
    print " "x$spaces, $_;
}

open OUT, ">$label_fn";
print OUT "$_=>$labels{$_}\n" for (keys %labels);
close OUT;

sub parsetoc {
    my $filename = shift;
    local $_;
    $filename =~ s/\.tex/.toc/;
    open INPUT, "<$filename" or die "Couldn't read table of contents!\n";
    print "## Contents\n";
    while (<INPUT>) {
	my ($num, $title, $ref) = /\{([\d.]+)\}([^\}]+)\}\{\d+\}\{(\w+[\d.]+)\}/;
	my $mdref = "\L$title\E";
	$mdref =~ tr/ /-/;
	$mdref =~ tr!(/)!!d;
	$mdref = "markdown-header-".$mdref;
	$mdref{$ref} = $mdref;
	if (!($num =~ s/\d+\.//)) {
	    print "* [$title](#$mdref)\n";
	}
	elsif (!($num =~ s/\d+\.//)) {
	    print " "x4,"$num. [$title](#$mdref)\n";
	}
	else {
	    $num =~ s/\d+\.//;
	    print " "x8,"$num. [$title](#$mdref)\n";
	}
    }
    close INPUT;
    print "\n";
}

sub latex {
    my ($c, $a) = @_;
    return "" if ($c eq "vspace" or $c eq "usepackage"); 
    if ($c eq "begin") {
	return "" if ($a eq "document" or $a eq "itemize");
	if ($a eq "verbatim") { $verbatim = 1; return ""; }
	if ($a eq "figure") { $infigure=1; $fig++; return "#### Figure $fig\n"; }
    }
    if ($c eq "end") {
	return "" if ($a eq "document" or $a eq "itemize");
	if ($a eq "verbatim") { $verbatim = 0; return ""; }
	if ($a eq "figure") { $infigure=0; return ""; }
    }

    return "_${a}_" if ($c eq "textit");
    return "*${a}*" if ($c eq "textbf");
    if ($c eq "section") {
	$subsec = $subsubsec = 0;
	$sec++;
	return "## $a ##\n";
    }
    if ($c eq "subsection") {
	$subsec++;
	$subsubsec = 0;
	return "$subsec. ### $a ###\n";
    }
    if ($c eq "subsubsection") {
	$subsubsec++;
	return "$subsubsec. #### $a ####\n";
    }
    if ($c eq "label") { 
	if ($infigure) { $labels{$a} = "fig.$fig"; }
	else { 
	    $labels{$a} = "section.$sec";
	    $labels{$a} = "sub". $labels{$a}. ".$subsec" if ($subsec or $subsubsec);
	    $labels{$a} = "sub". $labels{$a}. ".$subsubsec" if ($subsubsec);
	}
	return "";
    }
    if ($c eq "ref") {
	if ($refs{$a} and $mdref{$refs{$a}}) {
	    my $name = $refs{$a};
	    $name =~ s/^\w+\.//;
	    return "[$name](#$mdref{$refs{$a}})";
	} elsif ($refs{$a} and $refs{$a} =~ /fig/) {
	    my $name = $refs{$a};
	    $name =~ s/^\w+\.//;  
	    return "[$name](#markdown-header-figure-$name)";
	}
	else { return "??"; }
    }
    return '(C)' if ($c eq "textcopyright");
    return "<$a>" if ($c eq "url");
    if ($c eq 'ttt' or $c eq 'texttt') {
	$a =~ s/\\_/_/g;
	return "`$a`";
    } 
    return $a;
}
