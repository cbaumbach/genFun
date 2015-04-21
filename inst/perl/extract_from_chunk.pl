#!/usr/bin/env perl
# extract_from_chunk.pl
use warnings;
use strict;
use Getopt::Long;
use IO::Uncompress::Gunzip qw($GunzipError);

# Extract snps from gzip-compressed chunk file.

sub usage {
    print STDERR "$0 -s snp_file -c chunk_file\n";
    exit 1;
}

my ($in, $out);                 # multi-purpose filehandles

# ======================================================================
# Parse command-line options.
# ======================================================================
Getopt::Long::Configure(qw(no_auto_abbrev no_ignore_case_always bundling));

my $snp_file;
my $chunk_file;

GetOptions('s=s' => \$snp_file, 'c=s' => \$chunk_file) or usage();

usage() unless -e $snp_file and -e $chunk_file;

# ======================================================================
# Read list of snps.
# ======================================================================
my %selected;
open $in, '<', $snp_file
    or die "Can't open $snp_file for reading: $!";
while (<$in>) {
    chomp;
    ++$selected{$_};
}
close $in;

# ======================================================================
# Search through chunk file.
# ======================================================================
my $z = new IO::Uncompress::Gunzip $chunk_file
    or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
while (<$z>) {
    my (undef, $snp, undef) = split / /, $_, 3;
    print $_ if $selected{$snp};
    delete $selected{$snp};
    last unless keys %selected;
}
