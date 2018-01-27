#!/usr/bin/perl
#
# deltasvm.pl: A simple script that calculates deltaSVM scores
#
# Copyright (C) 2015 Dongwon Lee
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Last revised 03/30/2015

my $refseqf = $ARGV[0];
my $svmwf = $ARGV[1];

my %rc;
$rc{"A"} = 'T'; $rc{"T"} = 'A'; $rc{"C"} = 'G'; $rc{"G"} = 'C';
sub revcomp {
    my $seq = shift @_;
    my $rcseq = "";
    my $l = length($seq)-1;
    while ($l >= 0) {
        $rcseq = $rcseq . $rc{substr($seq, $l, 1)};
        $l-=1;
    }

    return $rcseq;
}

my $line;

#2. read svmwfile
print STDERR "reading svmwfile $svmwf..\n";
my %svmscore;
open IN, "<$svmwf";
my $bias;
chomp($bias=<IN>);
while (chomp($line=<IN>)) {
    if ($line=~/^#/) {next;}
    if ($line=~/^bias/) {next;}

    my @f = split /\s+/, $line;
    if ($#f == 1) {
        $svmscore{$f[0]} = $f[1];
        $svmscore{revcomp($f[0])} = $f[1];
    } elsif ($#f == 2) {
        $svmscore{$f[0]} = $f[2];
        $svmscore{$f[1]} = $f[2];
    } else {
        print STDERR "[ERROR] unkown format of svm weight file. quit.\n";
        exit(0);
    }
    if ($kmerlen == 0) {
        $kmerlen = length($f[0]);
    }
}
close IN;

#3. scan sequence file
my $seqid_ref;
my $seq_ref;
open IN1, "<$refseqf";
open OUT, ">$outf";
print STDERR "calculating deltaSVM using $refseqf ..\n";
while ($seq_ref = <IN1>) {
	my $score = 0;
	chomp($seq_ref);
	$score = $svmscore{$seq_ref};
	#print $seq_ref, "\n";
	print $score, "\n";
}
close IN1;
