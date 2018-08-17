#!perl -w
=head1 #===============================================================================
#        USAGE: perl  <input.matrix> > <output.matrix>
#
#  DESCRIPTION: Log the exp table
#
#      OPTIONS: None
# REQUIREMENTS: The input patameter must be
#        NOTES: None
#       AUTHOR: Zhangxiaolong (ST_TC_SINGLECELL), Zhangxiaolong@genomics.cn
# ORGANIZATION: BGI shenzhen
#      VERSION: 1.0
#      CREATED: 2014.10.23
#     REVISION: ---
#===============================================================================
=cut
die `pod2text $0` unless @ARGV == 1;
open IN,shift or die $!;

my $header = <IN>;
$header =~ s/\./_/g;
print "$header";
while (<IN>){
        chomp;
        my @F = split;
        print $F[0];
        for (1..$#F){
                $F[$_] = log($F[$_] + 1)/log(2);
                print "\t$F[$_]";
        }
print "\n";
}


