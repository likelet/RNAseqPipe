#!/usr/bin/perl -w
use strict;
die ("usage: <Matrix file> <Design file> <Compare String> <Output file>") unless @ARGV >= 4;
print "#Matrix file ".$ARGV[0]." with Design file ".$ARGV[1]." with Compare String ".$ARGV[2]." with Output file ".$ARGV[3]."\n";

my @type = split "_vs_",$ARGV[2];
my %type0; my %type1;

open FH,"$ARGV[1]" or die;
while(<FH>){
	chomp;
	my @field=split "\t";
        if ($field[1] eq $type[0]){
                $type0{$field[0]} = $field[1];
        }
	if ($field[1] eq $type[1]){
                $type1{$field[0]} = $field[1];
        }
}

open FH,"$ARGV[0]" or die;
my @pos0; my @pos1;my %h;
while(<FH>){
	chomp;
	my @field=split "\t";
        if ($field[1] eq "GeneName"){
                for (my $i = 3;$i<=$#field;$i++){
                        if (defined($type0{$field[$i]})){
                                push (@pos0,$i);
                        }
                        if (defined($type1{$field[$i]})){
                                push (@pos1,$i);
                        }
                }
        }else{
              	my $test = 0; my $ctr = 0;
                foreach my $k (@pos0){
                        $test += $field[$k] + 1;
                }
                foreach my $k (@pos1){
                        $ctr += $field[$k] + 1;
                }
                my $fc = log($test/$ctr)/log(2);
                $h{$field[1]} = $fc;
        }
}

my @keys = reverse sort { $h{$a} <=> $h{$b} } keys(%h);
my @vals = @h{@keys};
open OUT,"> $ARGV[3]" or die;
foreach my $i (@keys){
        print OUT $i."\t".$h{$i}."\n";
}
