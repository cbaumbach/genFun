#!/usr/bin/env perl
# obo2table.pl
use warnings;
use strict;

my %block;                      # collect data related to go term block
my $in_block = 0;               # Are we in a go term block?

$, = "\t";                      # output column separator
$\ = "\n";                      # output record separator

# Columns that go into output table.
my @columns = qw( id name namespace alt_id is_a part_of regulates
                  negatively_regulates positively_regulates def );

# Table dispatch on tag value.
my %deal_with = (
    id        => sub { $block{id} = substr $_[0], 4, 10 },
    name      => sub { $block{name} = substr $_[0], 6 },
    namespace => sub { $block{namespace} = substr $_[0], 11 },
    def       => sub { ($block{def} = $_[0]) =~ s/.*?"(.*?)".*/$1/ },
    alt_id    => sub { push @{$block{alt_id}}, substr $_[0], 8 },
    is_a      => sub { push @{$block{is_a}}, substr $_[0], 6, 10 },
    relationship => sub {
        $_[0] =~ /relationship: (\w+) (GO:\d+)/
            && push @{$block{$1}}, $2
            or die "unknown relationship: $_[0]\n";
    },
    is_obsolete => sub {
        if ($_[0] eq 'is_obsolete: true') {
            %block = ();        # clear block
            $in_block = 0;      # ignore rest of block
        }
        else {
            die "non-true obsolete value: $_[0]\n";
        }
    },
);

print join($,, @columns);       # header of output table

while (<>) {
    chomp;
    if (/^\[Term\]/) {
        $in_block = 1;
        next;
    }
    elsif (/^\s*$/) {
        print_block() if $in_block;
        $in_block = 0;
        next;
    }
    if ($in_block && /^([^:]+):/) {
        $deal_with{$1}->($_) if defined $deal_with{$1};
    }
}

sub print_block {
    my @cells = map {
        if (exists $block{$_}) {
            my $reftype = ref $block{$_};
            if    ($reftype eq '')      { $block{$_} }
            elsif ($reftype eq 'ARRAY') { join(',', @{$block{$_}}) }
            else { die "unknown reftype: $reftype\n" }
        } else {
            'NA';
        }
    } @columns;

    print join($,, @cells);
    %block = ();                # clear block
}
