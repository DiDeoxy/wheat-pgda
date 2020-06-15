use strict;
use warnings;

use Graph::Cliques::Bron_Kerbosch qw/ get_cliques /;

# Read the data into an array of arrays, converting from the question's R
# output. Each element of @edges contains a pair of nodes of the graph
#
my @edges;
while ( <DATA> ) {
    my @pair = split;
    next unless @pair > 2 and shift( @pair ) =~ /\[/;
    push @edges, \@pair;
}

# Call the utility function to get a list of cliques
#
my @groups = get_cliques( \@edges );

# Extract the hash keys to change the array of hashes into an array of sorted
# arrays, then sort the array first by the size of the clique and then by the
# first value in each group
#
$_ = [ sort { $a <=> $b } @$_ ] for @groups;
@groups = sort { @$a <=> @$b or $a->[0] <=> $b->[0] } @groups;

my $file = "curated_names.txt";
open (my $in, '<', $file) or die "Could not open file '$file' $!";
chomp(my @names = <$in>);
close $in;

print join( ' ', map { sprintf '%s', $names[$_ - 1] } @$_ ), "\n" for @groups;



__DATA__
 [1,] 345  23
 [2,] 352  23
 [3,] 355  23
 [4,] 266  29
 [5,] 111  31
 [6,] 160  36
 [7,] 183  36
 [8,] 239  36
 [9,]  38  37
[10,]  37  38
[11,]  41  40
[12,] 366  40
[13,]  40  41
[14,] 366  41
[15,]  54  53
[16,]  53  54
[17,]  67  66
[18,]  66  67
[19,] 332  84
[20,] 297  89
[21,]  94  93
[22,]  93  94
[23,] 140 105
[24,]  31 111
[25,] 105 140
[26,] 207 141
[27,] 147 146
[28,] 146 147
[29,]  36 160
[30,] 239 160
[31,] 169 168
[32,] 168 169
[33,]  36 183
[34,] 309 197
[35,] 141 207
[36,] 385 218
[37,] 301 221
[38,]  36 239
[39,] 160 239
[40,]  29 266
[41,]  89 297
[42,] 221 301
[43,] 197 309
[44,] 354 314
[45,] 322 321
[46,] 321 322
[47,] 327 326
[48,] 326 327
[49,]  84 332
[50,]  23 345
[51,] 352 345
[52,] 355 345
[53,]  23 352
[54,] 345 352
[55,] 355 352
[56,] 314 354
[57,]  23 355
[58,] 345 355
[59,] 352 355
[60,] 358 357
[61,] 360 357
[62,] 361 357
[63,] 357 358
[64,] 360 358
[65,] 361 358
[66,] 357 360
[67,] 358 360
[68,] 361 360
[69,] 357 361
[70,] 358 361
[71,] 360 361
[72,]  40 366
[73,]  41 366
[74,] 218 385