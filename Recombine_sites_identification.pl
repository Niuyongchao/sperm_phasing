#!/usr/bin/perl
unless (@ARGV){

	print "perl $0 <*.gt file> <reference chromosome length list (Fromat: chromosome01	45064769)> <Window_size (15) odd_num> <cut_value (13)> <min_bin_size (300000)>\n";
	exit 0;
}
use strict;
###########################ÌáÈ¡È¾É«Ìå±êºÅÐÅÏ¢################
my $head_line_of_chr_info=`head -1 $ARGV[1]`;
my $chr_type=(split/\s+/,$head_line_of_chr_info)[0];
############################################################## 
my $pp1;my $pp2;

my $min_bin_size=$ARGV[4];
#########################ÅÐ¶Ï´°¿ÚµÄ¿¨·½ÁÙ½çÖµ########################################################
my $win_size=$ARGV[2];
my $mid=int($win_size/2);
my $MID=$mid+1;
my $cut_off=0;
foreach my $A_num(1..$win_size-1) {
	my $half=$win_size/2; 
	my $B_num=$win_size-$A_num; 
	my $x2=((abs($A_num-$half)-0.5)**2)/$half + ((abs($B_num-$half)-0.5)**2)/$half; 
	my $cha=$A_num-$B_num; 
	if ($x2>=3.841 and $cha > 0){
		$cut_off=$A_num;  
		last;	
	}
}
my $cut_off=$ARGV[3];
print "Window size is $win_size\nCut value is $cut_off\nMin bin size is $min_bin_size\n";
###############################################################################

my (@hang, @array, @chrom, @bin, @len_detail, @temp);
my ($parent1, $parent2, $head, $head_edge, $max_n);
my ($hetero_12, $hetero_key, $hetero_start, $hetero_end, $line, $line1, $chromosome1, $cstart, $cend, $number1, $chrom_len );
my ($lane_n, $bstart, $bend, $origin, $number, $c, $h, $len, $round, $name);

$parent1=0, $parent2=0, $head=0, $head_edge=0, $max_n=0;
$hetero_12="", $hetero_key="", $hetero_start="", $hetero_end="", $line="", $line1="", $chromosome1="", $cstart="", $cend="", $number1="", $chrom_len="";
$lane_n=0, $bstart="", $bend="", $origin="", $number="", $c=0, $h=0, $len=0, $round=0, $name="";
#####################################################
#Judge the edges of every bin.
#####################################################
if ($ARGV[0]=~/.gz/){open INPUT,"gzip -dc $ARGV[0] |" or die "$!";}else{open INPUT,"$ARGV[0]" or die "$!";}
open OUT1, ">$ARGV[0].win$win_size.edge";
my $n=0;
my $m=1;
$chrom[0]=0;
my $chromosome=$chr_type;

while (<INPUT>) {   # store rlt-file;
	$line = $_;
	chomp($line);
	@hang = (split(/\s+/,$line))[0,1,8,9,-1] unless (/UNKNOWN/ or /Pt/ or /Mt/);
	$lane_n=$#hang;
	for (0..($lane_n)){
		$array[$n][$_]=$hang[$_]; #$n from 0 to $n-1
	}
	$n++; # total_line
	$hang[2]=~s/[a-z]//gi;
	$hang[2]=~s/^0//;
	$chrom[$hang[2]]++;#¸÷È¾É«ÌåreadsÊýÄ¿
}
print "total gt file line: $n\n";
################################´¢´æÈ¾É«Ìå³¤¶ÈºÍÊýÁ¿ÐÅÏ¢#########################
open INPUT2,"<$ARGV[1]" or die "$!";
my $count_total_chr_num=0;
while (<INPUT2>) {    # store chr_length-info
	$line1=$_; 
	chomp($line1);
	($chromosome1,$chrom_len)=split(/\s+/,$line1);
	$number1=$chromosome1;
	$number1=~s/[a-z]//gi;
	$number1=~s/^0//;
	$len_detail[$number1]=$chrom_len;
	$count_total_chr_num++;
}
my $all_chromosome=$count_total_chr_num; #chr num
close INPUT2;
##############################################################

for (1..$all_chromosome){
	$chrom[$_]=$chrom[$_]+$chrom[$_-1]; #È¾É«ÌåreadÊýÀÛ»ý
}

for (1..$all_chromosome){
	$c=$_;
	for (($chrom[$c-1]+$mid)..($chrom[$c]-$MID)) {
		my $win_start=$_;  #´ÓµÚ8ÐÐ¿ªÊ¼£¬Ä©Î²ÉÙ7ÐÐ£¬ÒòÎªÕû¸öÊÇ´Ó0¿ªÊ¼µÄ
		for ((($win_start)-$mid)..($win_start+$mid)){ #¼ÆËã´°¿ÚÄÚµÄ»ùÒòÐÍÊýÄ¿
			if ($array[$_][0] eq "P1"){
				$parent1++;
			}elsif($array[$_][0] eq "P2"){
				$parent2++;
			}
		}
		
		if (($parent1>=$parent2)){
			$array[$win_start][5]="parent1";
		}elsif($parent1<$parent2){
			$array[$win_start][5]="parent2";
		}
		if ($parent1>=$cut_off){
			$array[$win_start][6]="parent1";
			$array[$win_start][7]=$parent1.":".$parent2;
		}elsif ($parent2>=$cut_off){
			$array[$win_start][6]="parent2";
			$array[$win_start][7]=$parent1.":".$parent2;
		}else {
			$array[$win_start][6]="hetero";
			$array[$win_start][7]=$parent1.":".$parent2;
		}
### Get the half-win of every chromosome's fisrt bin.
		if ($win_start eq $chrom[$c-1]+$mid){ #´´°¿Ú´Ó7¿ªÊ¼Ê±
			my $f_bin_p1=0;			
			my $f_bin_p2=0;
			for ($chrom[$c-1]..$win_start-1){#0µ½6
				if ($array[$_][0] eq "P1"){
   	            	$f_bin_p1++;
	   	        }elsif($array[$_][0] eq "P2"){
	   	            $f_bin_p2++;
   	         	}
			}
				if (($f_bin_p1>$f_bin_p2)){
					for ($chrom[$c-1]..$win_start-1){
						$array[$_][5]="parent1";
					}
			    }elsif($f_bin_p1<$f_bin_p2){
		            for ($chrom[$c-1]..$win_start-1){
                        $array[$_][5]="parent2";
                    }
   			    }
		        if ($f_bin_p1>=$cut_off/2){
					for ($chrom[$c-1]..$win_start-1){
		            	$array[$_][6]="parent1";
		           		$array[$_][7]=$f_bin_p1.":".$f_bin_p2;
					}
		        }elsif ($f_bin_p2>=$cut_off/2){
		        	for ($chrom[$c-1]..$win_start-1){
                        $array[$_][6]="parent2";
                        $array[$_][7]=$f_bin_p1.":".$f_bin_p2;
                    }
				}else {
		            for ($chrom[$c-1]..$win_start-1){
						$array[$_][6]="hetero";
			            $array[$_][7]=$f_bin_p1.":".$f_bin_p2;
					}
		        }
        }
		
#######################################################################################################		
### Get the half-win of every chromosome's last bin.
		my $i=0;
		my $j=0;		
		if ($win_start eq $chrom[$c]-$MID){
			my $f_bin_p1=0;
            my $f_bin_p2=0;
            for ($i=$chrom[$c]-1;$i>=$win_start+1;$i--){
                if ($array[$i][0] eq "P1"){
                    $f_bin_p1++;
                }elsif($array[$i][0] eq "P2"){
                    $f_bin_p2++;
                }
            }
            if (($f_bin_p1>$f_bin_p2)){
                for ($i=$chrom[$c]-1;$i>=$win_start+1;$i--){
                    $array[$i][5]="parent1";
                }
            }elsif($f_bin_p1<$f_bin_p2){
                for ($i=$chrom[$c]-1;$i>=$win_start+1;$i--){
                    $array[$i][5]="parent2";
                }
            }
            if ($f_bin_p1>=$cut_off/2){
                for ($i=$chrom[$c]-1;$i>=$win_start+1;$i--){
                    $array[$i][6]="parent1";
                    $array[$i][7]=$f_bin_p1.":".$f_bin_p2;
                }
            }elsif ($f_bin_p2>=$cut_off/2){
                for ($i=$chrom[$c]-1;$i>=$win_start+1;$i--){
                    $array[$i][6]="parent2";
                    $array[$i][7]=$f_bin_p1.":".$f_bin_p2;
                }
            }else {
                for ($i=$chrom[$c]-1;$i>=$win_start+1;$i--){
                    $array[$i][6]="hetero";
                    $array[$i][7]=$f_bin_p1.":".$f_bin_p2;
                }
            }
        }

		$parent1=0;
		$parent2=0;
	}
}

### Adjust the edges.
###################################################
my $key1=1;
my $key2=1;
my $max=0;
for (0..$chrom[$all_chromosome]) {
	my $edge=$_;
	if (($array[$edge][5] ne "") && ($array[$edge+1][5] ne "") && ($array[$edge][5] ne $array[$edge+1][5]) && ($array[$edge][2] eq $array[$edge+1][2]) && ($array[$edge][6] ne "")){ ####
		for(1..30){
			my $temp=$_-15;
			for(1..14){
				if ($array[$edge-$_+$temp][0] eq $array[$edge+$temp][0]){
					$key1++; #´Ó1¿ªÊ¼£¬++ºó±äÎª2
				}elsif($array[$edge-$_+$temp][0] ne $array[$edge+$temp][0]){
					last;
				}
			}
			for(2..15){
				if ($array[$edge+$_+$temp][0] eq $array[$edge+1+$temp][0]){
					$key2++;
				}elsif($array[$edge+$_+$temp][0] ne $array[$edge+1+$temp][0]){
					last;
				}
			}
			$array[$edge+$temp][8]=$key1+$key2;
			if ($array[$edge+$temp][8]>$max){
				$max=$array[$edge+$temp][8];
				$max_n=$edge+$temp;
			}
			$key1=1;
			$key2=1;
		}
		if (($array[$edge][8]>=25) && ($array[$max_n+1][6] ne "")){
			$array[$max_n][9]=$array[$edge][5];
			$array[$max_n+1][9]=$array[$edge+1][5];
			$max=0;
		}
		if (($array[$edge][8]>=25) && ($array[$max_n+1][6] eq "")){
			$array[$edge][9]=$array[$edge][5];
			$array[$edge+1][9]=$array[$edge+1][5];
			$max=0;
		}
		elsif($array[$edge][8]<25){
			$array[$edge][9]=$array[$edge][5];
			$array[$edge+1][9]=$array[$edge+1][5];
			$max=0;
		}
	}
	if (($array[$edge][5] ne "") && ($array[$edge+1][5] ne "") && ($array[$edge][5] ne $array[$edge+1][5]) && ($array[$edge][2] eq $array[$edge+1][2]) && ($array[$edge][6] eq "")){
		$array[$edge][9]=$array[$edge][5];
		$array[$edge+1][9]=$array[$edge+1][5];
	}
}
###################################################
### Judge the heterozygous region.
RRR3: for (0..$chrom[$all_chromosome]) {
	if (($array[$_][6] eq "hetero") && ($array[$_-1][6] ne "hetero")){
		$hetero_start=$_;
		$hetero_key=0;
		$hetero_12=$array[$hetero_start][5];
		next RRR3;
	}
	if (($array[$_][6] eq "hetero") && ($array[$_-1][6] eq "hetero")){
		if ($hetero_12 eq $array[$_][5]){
			next RRR3;
		}elsif ($hetero_12 ne $array[$_][5]){
			$temp[$hetero_key]=$_;
			$hetero_key++;
			
			$hetero_12=$array[$_][5];
			next RRR3;
		}
	}
	if (($array[$_][6] ne "hetero") && ($array[$_-1][6] eq "hetero")){
		$hetero_end=$_-1;
		if ($hetero_key==2 && ($temp[1]-$temp[0]<3)){
			$array[$hetero_start-1][9]=$array[$hetero_start-1][6];
			$array[$hetero_end+1][9]=$array[$hetero_end+1][6];
			$array[$hetero_start][9]="heterozygo";
			$array[$hetero_end][9]="heterozygo";
			for ($hetero_start+1..$hetero_end-1){
				$array[$_][9]="";
			}
		}
		if ($hetero_key>2){
			$array[$hetero_start-1][9]=$array[$hetero_start-1][6];
			$array[$hetero_end+1][9]=$array[$hetero_end+1][6];
			$array[$hetero_start][9]="heterozygo";
			$array[$hetero_end][9]="heterozygo";
			for ($hetero_start+1..$hetero_end-1){
				$array[$_][9]="";
			}
		}
		$hetero_key=0;
		$hetero_start="";
		$hetero_end="";
		next RRR3;
	}
}
##############################	
for (0..$chrom[$all_chromosome]){
	$h=$_;
	for (0..9){
		if ($_==0){
			print OUT1 $array[$h][$_];
		}else{
			print OUT1 "\t".$array[$h][$_];
		}
	}
	print OUT1 "\n";
}			
			
close INPUT;
close OUT1;

#####################################################
#Get every bin based on the edges.#
#####################################################

open IN2,"<$ARGV[0].win$win_size.edge" or die "$!";
my $filename=$ARGV[0].".win$win_size.edge";
my ($prefix,$cenfix,$suffix)=split(/\./,$filename);
open OUT2, ">$ARGV[0].win$win_size.bin.temp";
$n=0;
$m=1;
@array=();
@chrom=();
$chrom[0]=0;
@bin=();
my $ij="";
my $sn=0;
$chromosome=$chr_type;
while (<IN2>) {
	$line = $_;
	chomp($line);
	@hang = split(/\t/,$line);
	$lane_n=9;
	for (0..($lane_n)){
		$array[$n][$_]=$hang[$_];
	}
	$n++;
}
close IN2;
#´ÎÄ¼þ´æÈë¹þÏ£array
for (0..$n) {
	if ($array[$_][2] eq "$chromosome"){
		$chrom[$m]++;
	}else{
		$chromosome=$array[$_][2];
		$chrom[$m+1]++;
		$m++;
	}
}
#Í³¼ÆÃ¿ÌõÈ¾É«ÌåÐÐÊý
for (1..$all_chromosome){
	$chrom[$_]=$chrom[$_]+$chrom[$_-1];
}
#ÀÛ¼Ó
for (1..$all_chromosome){
	$c=$_;
	
	for ($chrom[$c-1]..$chrom[$c]-1) { #·ÖÈ¾É«Ìå´¦Àí
		if ($_==$chrom[$c-1]){
			$bin[$sn][0]=$array[$_][2];
			$bin[$sn][1]=1;
			$ij=$array[$_][5]; #È¾É«ÌåµÚÒ»ÐÐ¸³Öµ
		}
#########################################################################################################################################
#		if (($array[$_][9] ne "") && ($array[$_+1][9] ne "") && ($array[$_][9] ne $array[$_+1][9]) ){ ##BUGBUGBUGBUGBUG
		if (($array[$_][9] ne "") && ($array[$_+1][9] ne "") && ($array[$_][9] ne $array[$_+1][9]) && ($array[$_][2] eq $array[$_+1][2]) ){


####################################################################################################################################
#			$bin[$sn][2]=int(($array[$_][3]+$array[$_+1][3])/2);
			$bin[$sn][2]=int($array[$_][3]+($array[$_+1][3]-$array[$_][3])/2);
#####################################################################################################################################
#			print "$array[$_][2]\t$array[$_][3]\t$array[$_+1][3]\t$bin[$sn][2]\n";
			$bin[$sn][3]=$array[$_][9];
			$bin[$sn][4]=$array[$_][1];
			$bin[$sn][5]=$bin[$sn][2]-$bin[$sn][1]+1; # chromosome01    1       29689708        parent2 FCD1LYHACXX:3:2102:5384:11632#GAATATGG/2        29689708
			$sn++;
			$bin[$sn][0]=$array[$_][2];####È¾É«Ìå
			$bin[$sn][1]=$bin[$sn-1][2]+1;
			$ij=$array[$_+1][9];
		}
		if ($_==$chrom[$c]-1){#È¾É«Ìå×îºóÒ»ÐÐ
			$bin[$sn][2]=$len_detail[$c]; #È¾É«Ìå³¤¶È£
			$bin[$sn][3]=$ij;
			$bin[$sn][4]="chr_end";
			$bin[$sn][5]=$bin[$sn][2]-$bin[$sn][1]+1;
			$sn++;
		}
	}
}

for (0..($sn-1)){
	$h=$_;
	for (0..5){
		if ($_==0){
			print OUT2 $bin[$h][$_];
		}else{
			print OUT2 "\t".$bin[$h][$_];
		}
	}
	print OUT2 "\n";
}
close OUT2;

#####################################################
# Filtering bins which are smaller than XXbp.#
#####################################################
open IN3,"<$ARGV[0].win$win_size.bin.temp" or die "$!";
$n=0;
my $key="";
my @bin=();
while (<IN3>) {
	$line = $_;
	@hang = split(/\s+/,$line);
	$lane_n=4;
	for (0..($lane_n)){
		$bin[$n][$_]=$hang[$_];
	}
	$n++;
}
close IN3;

RRRR:for(1..100){ #######Ñ­»·´ÎÊý
	$round=$_;
	if ($key==1){
		open IN4,"<$ARGV[0].win$win_size.bin" or die "$!";
		$n=0;
		while (<IN4>) {
			$line = $_;
			@hang = split(/\s+/,$line);
			$lane_n=4;
			for (0..($lane_n)){
				$bin[$n][$_]=$hang[$_];
			}
			$n++;
		}
		close IN4;
		$key=0;
	}
	open OUT3, ">$ARGV[0].win$win_size.bin";
	my $skip=0;
	$key=0;
	for (0..($n-1)){
		$h=$_+$skip;
		if ($h<$n-1) {
			if (($bin[$h+1][2]-$bin[$h+1][1]+1)>=$min_bin_size){ #ÏÂÒ»ÐÐ´óÓÚ¹Ì¶¨Æ¬¶Î£¬Ö±½ÓÊä³ö¡£
				for (0..4){
					if ($_==0){
						print OUT3 $bin[$h][$_];
					}else{
						print OUT3 "\t".$bin[$h][$_];
					}
				}
				$len=$bin[$h][2]-$bin[$h][1]+1;
				print OUT3 "\t".$len;
				print OUT3 "\n";
				if ($len<$min_bin_size){
                        $key=1;
                    }
###################################################################################################################
			}elsif((($bin[$h+1][2]-$bin[$h+1][1]+1)<$min_bin_size) && ($bin[$h+2][3] eq $bin[$h][3])){ #ÏÂÒ»ÐÐÐ¡ÓÚ¹Ì¶¨Æ¬¶Î£¬ÏºÍÏÂÁ½ÐÐÀàÐÍÒ»Ñù
				print OUT3 $bin[$h][0]."\t";
				print OUT3 $bin[$h][1]."\t";
				if (($bin[$h+1][4] ne "chr_end") && ($bin[$h+1][1] ne "1")){#²»ÊÇ½áÎ²£¬²»ÊÇÏÂÒ»È¾É«ÌåµÄ¿ªÍ·
					print OUT3 $bin[$h+2][2]."\t";
					print OUT3 $bin[$h+2][3]."\t";
					print OUT3 $bin[$h+2][4]."\t";
					$len=$bin[$h+2][2]-$bin[$h][1]+1;
					print OUT3 $len;
					print OUT3 "\n";
					if ($len<$min_bin_size){
						$key=1;
					}
					$skip=$skip+2;
				}elsif ($bin[$h+1][1] eq "1"){
					print OUT3 $bin[$h][2]."\t";
					print OUT3 $bin[$h][3]."\t";
					print OUT3 $bin[$h][4]."\t";
					$len=$bin[$h][2]-$bin[$h][1]+1;
					print OUT3 $len;
					print OUT3 "\n";
#############################################################
					if ($len<$min_bin_size){
                        $key=1;
                    }
#############################################################

				}elsif ($bin[$h+1][4] eq "chr_end"){
					print OUT3 $bin[$h+1][2]."\t";
					print OUT3 $bin[$h][3]."\t";
					print OUT3 $bin[$h+1][4]."\t";
					$len=$bin[$h+1][2]-$bin[$h][1]+1;
					print OUT3 $len;
					print OUT3 "\n";
					if ($len<$min_bin_size){
                        $key=1;
                    }
					$skip=$skip+1;
				}
#####################################################################################################################
			}elsif((($bin[$h+1][2]-$bin[$h+1][1]+1)<$min_bin_size) && ($bin[$h+2][3] ne $bin[$h][3])){ #ÏÂÒ»ÐÐÐ¡ÓÚ¹Ì¶¨Æ¬¶Î£¬ºÍÏÂÏÂÐÐ²»µÈ¡£
#				for (0..4){
#					if ($_==0){
#						print OUT3 $bin[$h][$_];
#					}else{
#						print OUT3 "\t".$bin[$h][$_];
#					}
#				}
				print OUT3 $bin[$h][0]."\t";
                print OUT3 $bin[$h][1]."\t";
				print OUT3 $bin[$h][2]."\t";
                print OUT3 $bin[$h][3]."\t";
                print OUT3 $bin[$h][4]."\t";
                $len=$bin[$h][2]-$bin[$h][1]+1;
                print OUT3 $len;
                print OUT3 "\n";
				if ($len<$min_bin_size){
                        $key=1;
                    }
                    
                if (($bin[$h+1][4] ne "chr_end") && ($bin[$h+1][1] ne "1")){#²»ÊÇ½áÎ²£¬²»ÊÇÏÂÒ»È¾É«ÌåµÄ¿ªÍ·
                    print OUT3 $bin[$h+1][0]."\t";
                    print OUT3 $bin[$h+1][1]."\t";
					print OUT3 $bin[$h+2][2]."\t";
                    print OUT3 $bin[$h+2][3]."\t";
                    print OUT3 $bin[$h+2][4]."\t";
                    $len=$bin[$h+2][2]-$bin[$h+1][1]+1;
                    print OUT3 $len;
                    print OUT3 "\n";
                    if ($len<$min_bin_size){
                        $key=1;
                    }
                    $skip=$skip+2;
				}
#				$len=$bin[$h][2]-$bin[$h][1]+1;
#				print OUT3 "\t".$len;
#				print OUT3 "\n";
			}
		}
		elsif ($h==$n-1) {#×îºóÒ»ÐÐ
			for (0..4){
				if ($_==0){
					print OUT3 $bin[$h][$_];
				}else{
					print OUT3 "\t".$bin[$h][$_];
				}
			}
			$len=$bin[$h][2]-$bin[$h][1]+1;
			print OUT3 "\t".$len;
			print OUT3 "\n";
			if ($len<$min_bin_size){
				$key=1;
			}
		}
	}
	close OUT3;
	last RRRR if $key==0;
}

#####################################################
#Draw a figure in PNG format indicating bins and SNPs of the RIL.
#####################################################
open OUT4, ">$ARGV[0].png";
use GD; 

my $xtimes=5;
my $im = new GD::Image(10000,4000); 
my $black = $im->colorAllocate(0,0,0); 
my $white = $im->colorAllocate(255,255,255); 
my $red = $im->colorAllocate(255,0,0);
my $blue = $im->colorAllocate(0,0,255);
my $green = $im->colorAllocate(0,255,0);
my $hetero = $im->colorAllocate(255,255,0);
my $background=$im->colorAllocate(200,200,200);
$im->fill(10,10,$white); 
$im->rectangle(0,20,7000,30,$black);
$im->fill(2500,25,$black);
for (1..70){
	$im->line($_*100,5,$_*100,20,$black);
}
for (1..6){
	my $scale_temp=$_*10*$xtimes;
	my $scale=$scale_temp." Mb";
	$im->string(gdGiantFont,$_*1000-25,50,"$scale",$black);
}
open INPUT2,"<$ARGV[1]" or die "$!";
while (<INPUT2>) {
	$line1=$_; 
	chomp($line1);
	($chromosome1,$chrom_len)=split(/\s+/,$line1);
	$number1=$chromosome1;
	$number1=~s/[a-z]//gi;
	$number1=~s/^0//;
	$cend=int(($chrom_len/10000/$xtimes)+0.5);
	$im->rectangle(0,$number1*150,$cend,($number1*150+100),$background);
	$im->fill($cend*0.5,($number1*150+50),$background); 
}
close INPUT2;
open IN5,"<$ARGV[0].win$win_size.bin.2" or die "$!";
while (<IN5>) {
	$line=$_; 
	chomp($line);
	($chromosome,$bstart,$bend,$origin,$name)=split(/\s+/,$line);
	$bstart=int($bstart/10000/$xtimes);
	$bend=int($bend/10000/$xtimes);
	$number=$chromosome;
	$number=~s/[a-z]//gi;
	$number=~s/^0//;
	if ($origin eq "parent1"){
		$im->rectangle($bstart,$number*150+21,$bend,($number*150+40),$blue);
		$im->fill(($bstart+$bend)*0.5,($number*150+30),$blue); 
	}elsif ($origin eq "parent2"){
		$im->rectangle($bstart,$number*150+21,$bend,($number*150+40),$red);
		$im->fill(($bstart+$bend)*0.5,($number*150+30),$red); 
	}
	elsif ($origin eq "heterozygo"){
		$im->rectangle($bstart,$number*150+21,$bend,($number*150+40),$hetero);
		$im->fill(($bstart+$bend)*0.5,($number*150+30),$hetero); 
	}
}
close IN5;
if ($ARGV[0]=~/.gz/){open INPUT4,"gzip -dc $ARGV[0] |" or die "$!";}else{open INPUT4,"$ARGV[0]" or die "$!";}

while (<INPUT4>) {
        next if /Contig/;
        $line=$_;
        chomp($line);
	
        my @temp=split /\s+/,$line;
	
        $bstart=int($temp[9]/10000/$xtimes);
        $number=$temp[8];
        $number=~s/[a-z]//gi;
        $number=~s/^0//;
        if ($temp[0] eq "P1"){
                $im->rectangle($bstart,$number*150+51,$bstart,($number*150+70),$blue);
                $im->fill(($bstart+$bstart)*0.5,($number*150+60),$blue);
		$pp1++;
        }elsif ($temp[0] eq "P2"){
                $im->rectangle($bstart,$number*150+71,$bstart,($number*150+90),$red);
                $im->fill(($bstart+$bstart)*0.5,($number*150+80),$red);
        	$pp2++;
	}
#        elsif ($temp[0] eq "X"){
#                $im->rectangle($bstart,$number*150+61,$bstart,($number*150+80),$hetero);
#                $im->fill(($bstart+$bstart)*0.5,($number*150+70),$hetero);
#        }
}
close INPUT4;
	binmode OUT4; 
	print OUT4 $im->png; 
close OUT4;
print "P1: $pp1\tP2: $pp2\n";

exit; 

