#!/usr/bin/perl -w
use strict;

my $mm=-M "$ENV{srcdir}/.done";
exit 0 if (-e "$ENV{srcdir}/.done" && (-M "$ENV{srcdir}/generate_data.pl" > -M "$ENV{srcdir}/.done"));

my $seed = 15551;
sub rr {
    my ($m)=@_;

    $seed = int(($seed*1103515245+12345)&0xffffffff);
    if (defined($m)) {
return $seed % $m;
    } else {
        return ($seed & 0xffffff) / 0x1000000;
    }
}

# Loads a fasta file into a hash
sub load_fasta {
    my ($fn) = @_;
    my %seqs;

    open(my $fd, "<", $fn) || die "$fn: $!\n";
    my $name = undef;
    my $s = "";
    while (<$fd>) {
	chomp($_);

	if (/^>(\S+)/) {
	    if (defined($name)) {
		$seqs{$name} = $s;
	    }
	    $name = $1;
	    $s = "";
	} else {
	    $s .= $_;
	}
    }
    $seqs{$name} = $s;

    close($fd);

    return \%seqs;
}

#---- Load seq
print "Loading ce.fa\n";

my $seqs = load_fasta("$ENV{srcdir}/data/ce.fa");
my @names = keys(%$seqs);
my %len;
my $n;
my @base = qw/A C G T/;


if (! -w "$ENV{srcdir}/data") {
    chmod(0755, "$ENV{srcdir}/data") || die "$ENV{srcdir}/data: $!";
    chmod(0755, "$ENV{srcdir}")      || die "$ENV{srcdir}: $!";
}

#---- Generate sorted data
print "Generating ce#sorted.sam\n";
open(my $out, ">", "$ENV{srcdir}/data/ce#sorted.sam") ||
    die "$ENV{srcdir}/data/ce#sorted.sam: $!";

# Create @SQ headers
foreach (sort @names) {
    my $len = length($seqs->{$_});
    print $out "\@SQ\tSN:$_\tLN:$len\n";
    $len{$_} = $len;
}

# Sequence lines
$n = 1;
$len{"*"}=$len{$names[0]};
$seqs->{"*"}=$seqs->{$names[0]};
foreach my $chr ((sort @names),"*") {
    my $len = $len{$chr};
    my $s = $seqs->{$chr};
    for (my $i=0; $i<$len-100; $i++) {
	if (rr() < 0.5) { #50x coverage
	    my $dna = substr($s,$i,100);
	    $dna = substr($s, int(rr($len{"*"}-100)), 100) if ($chr eq "*");
	    for (my $j=0; $j<5; $j++) {
		substr($dna, 100*rr(), 1) = $base[4*rr()];
	    }

	    if ($chr eq "*") {
		print $out "Seq",$n++,"\t4\t$chr\t0\t0\t*\t*\t0\t0\t$dna\t*\n";
	    } else {
		print $out "Seq",$n++,"\t0\t$chr\t",$i+1,"\t2\t100M\t*\t0\t0\t$dna\t*\n";
	    }
	}
    }
}
close($out) || die;


#---- Generate unsorted data
print "Generating ce#unsorted.sam\n";
srand 15551;
open($out, ">", "$ENV{srcdir}/data/ce#unsorted.sam") ||
    die "$ENV{srcdir}/data/ce#unsorted.sam: $!";

# Create @SQ headers
foreach (sort @names) {
    my $len = length($seqs->{$_});
    print $out "\@SQ\tSN:$_\tLN:$len\n";
    $len{$_} = $len;
}

# Sequence lines
$n = 1;
for (my $i = 0; $i < 500000; $i++) {
    my $chr = $names[$#names*rr()];
    my $pos = int(rr() * ($len{$chr}-100));
    my $dna = substr($seqs->{$chr},$pos,100);
    for (my $j=0; $j<5; $j++) {
	substr($dna, 100*rr(), 1) = $base[4*rr()];
    }

    print $out "Seq",$n++,"\t0\t$chr\t",$pos+1,"\t2\t100M\t*\t0\t0\t$dna\t*\n";
}
close($out) || die;


#---- Tidy up
open(DONE, ">$ENV{srcdir}/.done")||die;
close(DONE);

