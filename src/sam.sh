#!/bin/bash
. $HMHOME/src/root.sh

sam.csize(){
        if [ ! -f $1.bai ];then
                samtools index $1;
        fi
        samtools idxstats $1 | awk -v OFS="\t" '$1 != "*" && $3 > 0 { print $1,$2;}'
}

# reference from https://samtools.github.io/hts-specs/SAMv1.pdf
sam_to_bed(){
	cat $1 | perl -ne 'chomp; my @a=split/\t/,$_;
		if($_=~/^@/){ next;}
		my $flag=$a[1];
		my $chrom=$a[2];
		my $start=$a[3]-1;
		my $mapq=$a[4]; # -10log10 Pr( wrong )
		my $cigar=$a[5];
		my $seq=$a[9];
		my $len=0;
		#my $strand="+"; if ( $flag & 16 ){ $strand="-"; }
		my $strand="+"; if ( $flag & (0x10) ){ $strand="-"; }

		my $gseq=""; # genomic sequence 
		my $offset=0;
		#\*|([0-9]+[MIDNSHPX=])+ 
		while($cigar=~/(\d+)([MIDNSHPX=])/g){ 
			my ($x,$c)=($1,$2);
			if($c=~/[MX=]/){
				$gseq .= substr($seq,$offset,$x);
				$offset += $x; $len += $x;
			}elsif($c=~/[D]/){
				$gseq .= "*"x$x;
				$len += $x;
			}elsif($c=~/[IS]/){
				$offset += $x;
			}elsif($c=~/[N]/){
				$gseq .= "."x$x;
			}else{
				## unknown
			}
			
		}
		my $end=$start+$len;
		print $chrom,"\t",$start,"\t",$end,"\t",$gseq,"\t",$mapq,"\t",$strand,"\n";
	'
}
test__sam_to_bed(){
## example from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002/figure/F1/
echo \
"@HD	VN:1.0	SO:coordinate
@SQ	SN:chr1	LN:249250621
r1	163	chr1	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG *
r2	0	chr1	9	30	3S6M1P1I4M	*	0	0	AAAAGATAAGGATA	*
r3	0	chr1	9	30	5H6M	*	0	0	AGCTAA	*	NM:i:1
r4	0	chr1	16	30	6M14N5M	*	0	0	ATAGCTTCAGC	*
r3	16	chr1	29	30	6H5M	*	0	0	TAGGC	*	NM:i:0
r1	83	chr1	37	30	9M	=	7	-39	CAGCGCCAT	*" \
| sam_to_bed - > obs
echo \
"chr1	6	22	TTAGATAAGATA*CTG	30	+
chr1	8	18	AGATAAGATA	30	+
chr1	8	14	AGCTAA	30	+
chr1	15	26	ATAGCT..............TCAGC	30	+
chr1	28	33	TAGGC	30	-
chr1	36	45	CAGCGCCAT	30	-" > exp
check exp obs
rm -rf exp obs
}

