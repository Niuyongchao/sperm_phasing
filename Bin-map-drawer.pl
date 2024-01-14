#!/usr/bin/perl
unless (@ARGV>0) {
	print "perl $0 <*.map> <Chr_length_file> <out_file_prefix>\n";
	exit 0;
}
use GD;
open IN, "$ARGV[0]" || die "$!";
##################--Get Chr_legth info--########################
open CHR, "$ARGV[1]" || die "$!";
my $c=0;
while (<CHR>){
	chomp;
	@_=split/\s+/,$_;
	$acm=$acm+$_[1];
	push @chr_len, $acm;
	if ($c==0){
		$hash{$_[0]}=0;
	}else{
		$hash{$_[0]}=$chr_len[$c-1];
	}
	$c++;
}
close CHR;
###################--Draw bin_map--###############################
open OUT, ">$ARGV[2].png" ||die "$!";
my $im = new GD::Image(6000,5000);
my $white = $im->colorAllocate(255,255,255);
my $black = $im->colorAllocate(0,0,0);       
#my $red = $im->colorAllocate(255,0,0);      
my $red = $im->colorAllocate(255,153,0);
#my $blue = $im->colorAllocate(0,0,255);
my $blue = $im->colorAllocate(0,153,204);
my $gray = $im->colorAllocate(192,192,192);
my $yellow = $im->colorAllocate(255,255,0);
my $head_line=<IN>;
my $first_x=100;
my $x=100;
my $color;
while (my $map_line=<IN>){
	chomp $map_line;
	my @temp=split/\s+/,$map_line;
	$y=100;	
	my $first_y=100;
	foreach $num(2..$#temp){
		$x=($temp[1]+ int($hash{$temp[0]}/100000+0.5))/4+100;
		$y=$y+20;
		if ($temp[$num]eq "A"){
			$color=$red;
		
		}elsif($temp[$num]eq "B"){
        	$color=$blue;
		}
		$im->filledRectangle($first_x,$first_y,$x,$y,$color);
		$first_y=$y;
	}
	$first_x=$x;
	$sample_num=$#temp-1;
}
close IN;
####################--Split CHR by white line--#######################
my $line_y=100;
my $l_y=$sample_num*20+100;
foreach $key(sort keys %hash){
	my $l_x=int ($hash{$key}/100000+0.5)/4+100;
	$im->filledRectangle($l_x-1,$line_y,$l_x+1,$l_y,$white);
	$line_x=$l_x;
}

binmode OUT;
print OUT $im->png;
close OUT;
