# sperm_phasing
The script provided in this study can identify recombination sites in sperm and construct a binmap
Recombine_sites_identification.pl can be used for identifying the recombine site by using the gt file.
Usage: Recombine_sites_identification.pl <*.gt file> <reference chromosome length list (Fromat: chromosome01        45064769)> <Window_size (15) odd_num> <cut_value (13)> <min_bin_size (300000)>
The format of the gt file:
P2      1_7405  2       3       4       5       6       7       chr01   7405    C
Only columns 1, 9, and 10 are meaningful. P2 means the genotype is same with parent 2. P1 means the genotype is same with parent 1.

Binmap_constructor.pl can generate binmap file by using the input files(*.bin files) generated by Recombine_sites_identification.pl.
Usage: perl Binmap_constructor.pl bin.file.list chr.length.file

Bin-map-drawer.pl is a script used for drawing the binmap.
Usage: perl Bin-map-drawer.pl <*.map> <Chr.length.file> <out_file_prefix>


