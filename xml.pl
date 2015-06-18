#!/usr/bin/perl 
use strict;
use warnings;

my $semantic = "SemanticTypes_2013AA.txt";
open (FILE, "<$semantic");
my @array;
my %semtypes;
while (my $line = <FILE>){
    @array = split /\|/, $line ;
    $semtypes{$array[0]} = $array[2];
    
}

my $file = "MetaMapOut.xml";

open (FILE, "<$file");
my @final; 
while(my $line = <FILE>){
   if ($line =~ m/<SemType>/){
     $line =~ s/<\w*>//g;
     $line =~ s/<\/\w*>//g; 
     $line =~ s/\s//g;
     $line =~ s/\n//g;
     push @final, $semtypes{$line}
   }
}
close (FILE);

print join ("", @final)
