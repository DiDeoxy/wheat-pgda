package Graph::Cliques::Bron_Kerbosch;

use strict;
use warnings;
use v5.8.3;

use Exporter qw/ import /;

our @EXPORT_OK = qw/ get_cliques /;

my ( %neighbours, @cliques );

sub get_cliques {

    my ( $edges ) = @_;
    %neighbours = ();
    @cliques = ();

    for my $edge ( @$edges ) {
        my ( $n1, $n2 ) = @$edge;
        $neighbours{$n1}{$n2} = 1;
        $neighbours{$n2}{$n1} = 1;
    }
    $_ = [ keys %$_ ] for values %neighbours;

    my ( %r, %p, %x );
    $p{$_} = 1 for map @$_, @$edges;

    _bron_kerbosch( \( %r, %p, %x ) );

    @cliques;
}

sub _bron_kerbosch {
    my ( $r, $p, $x ) = @_;

    unless ( %$p or %$x ) {
        push @cliques, [ keys %$r ];
        return;
    }

    for my $v ( _choose_pivot($p, $x) ) {

        my $nv = $neighbours{$v};

        my %r_ = ( %$r, $v => 1 );
        my %p_ = map { $_ => 1 } _intersect( [ keys %$p ], $nv);
        my %x_ = map { $_ => 1 } _intersect( [ keys %$x ], $nv);

        _bron_kerbosch( \( %r_, %p_, %x_ ) );

        delete $p->{$v};
        $x->{$v} = 1;
    }
}

sub _intersect {
    my ( $aa, $ab ) = @_;
    my %ab = map { $_ => 1 } @$ab;
    grep $ab{$_}, @$aa;
}

# Find an element u of P U X such that as many as possible of its
# neighbours fall in P
#
sub _choose_pivot {
    my ( $p, $x ) = @_;

    my @p = keys %$p;
    my @choice = @p;

    for my $u ( @p, keys %$x ) {
        my $nu = $neighbours{$u};
        my %nu = map +( $_ => 1 ), @$nu;
        my @subset = grep { not $nu{$_} } @p;
        @choice = @subset if @subset < @choice;
    }

    @choice;
}

1;
