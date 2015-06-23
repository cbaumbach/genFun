#!/usr/bin/env perl
# extract_from_chunk.pl
use warnings;
use strict;
use Getopt::Long;
use IO::Uncompress::Gunzip qw($GunzipError);
use IO::Compress::Gzip qw($GzipError);

my ($in, $out);                 # multi-purpose file handles

# Extract snps from gzip-compressed chunk file.

sub usage {
    print STDERR "$0 -s snp_file -c chunk_file -o output_file [-i column_file]\n";
    exit 1;
}

# ====================================================================
# Parse command-line options.
# ====================================================================
Getopt::Long::Configure(qw(no_auto_abbrev no_ignore_case_always bundling));

my $snp_file;
my $chunk_file;
my $output_file;
my $column_file;

GetOptions('s=s' => \$snp_file,
           'c=s' => \$chunk_file,
           'o=s' => \$output_file,
           'i=s' => \$column_file) or usage();

die "snp file does not exist: \"$snp_file\"\n"
    unless -e $snp_file;
die "chunk file does not exist: \"$chunk_file\"\n"
    unless -e $chunk_file;
die "output file does already exist: \"$output_file\"\n"
    if -e $output_file;
die "column file does not exist: \"$column_file\"\n"
    if defined $column_file && ! -e $column_file;

# ====================================================================
# Read selected columns.
# ====================================================================
my @columns;
if ($column_file) {
    open $in, '<', $column_file
        or die "Can't open $column_file for reading: $!\n";
    @columns = <$in>;
    close $in;

    # Convert to 0-based indexes.
    @columns = map { $_ - 1 } @columns;
}


# ====================================================================
# Read list of snps.
# ====================================================================
my %selected;
open $in, '<', $snp_file
    or die "Can't open $snp_file for reading: $!";
while (<$in>) {
    chomp;
    ++$selected{$_};
}
close $in;

# ====================================================================
# Search through chunk file.
# ====================================================================
my $zin = new IO::Uncompress::Gunzip $chunk_file
    or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
my $zout = new IO::Compress::Gzip $output_file
    or die "IO::Compress::Gzip failed: $GzipError\n";
while (my $line = <$zin>) {
    my (undef, $snp, undef) = split / /, $line, 3;

    if ($selected{$snp}) {
        if (@columns) {
            my @fields = split / /, $line;
            print $zout join(' ', @fields[@columns]);
        }
        else {
            print $zout $line;
        }
        delete $selected{$snp};
        last unless %selected;
    }
}
close $zin;
close $zout;

# Write snps that were not found to stdout.
print $_, "\n" for sort keys %selected;
