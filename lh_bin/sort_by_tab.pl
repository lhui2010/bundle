#!/usr/bin/perl -w 
while(<>){
	chomp;
	@e=split /\t/, $_; 
	$hash{$e[2]} = $_; 
	} 

for my $k(sort {$b<=>$a} keys %hash){
	print $hash{$k}, "\n";
	}
