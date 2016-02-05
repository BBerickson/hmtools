#!/bin/bash
. $HMHOME/src/root.sh
. $HMHOME/src/bed.sh
. $HMHOME/src/stat.sh
asp.view(){
usage="     
USAGE: $FUNCNAME <genes.bed12> <reads.bed12>
INPUT: genes (bed12) and reads (bed12)

OUTPUT:
            @
           / \
          *   #
         / \ / \	
    [    ]-[ ]-[   ]

"
	local tmpd=`mymktempd`;
	mycat $1 > $tmpd/a
	mycat $2 > $tmpd/b
	bed.intron  $tmpd/b | bed.flank - 1 1 0 > $tmpd/j
	run_R '
		jc=read.table("'$tmpd/j'",header=F);
		y=as.numeric(jc[,3]-jc[,2]);
		x=jc[,2:3];
		c=apply(x,1,mean);
		xlim=c(min(x[,1]),max(x[,2]));
		ylim=c(0,max(y));
		png("out.png")
		
		plot(NULL,ylim=ylim,xlim=xlim);
		segments(x[,1],0,c,y);
		segments(c,y,x[,2],0);
		
		dev.off();
	' log		

}
asp.view.test(){
echo \
"chr1	1000	2000	g1	0	+	1000	2000	0,0,0	3	100,100,100	0,500,900"> tmp.gene
echo \
"chr1   1050    1550    r1      1       +       1050    1550    255,0,0 2       50,50   0,450
chr1   1050    1950    r1      1       +       1050    1950    255,0,0 2       50,50   0,850"> tmp.read
asp.view tmp.gene tmp.read
rm tmp.read tmp.gene
}
asp.spi(){
usage=" 
 FUNCT: Calculate SPI
 OUTPU: bed6 a,b,c 2*a/(b+c)
 USAGE: $FUNCNAME <intron.bed6> <read.bed12> <strand> [<count>]
	<strand>: 0: null, 1:count on the same strand, 2: count on the opposite strand
	<count> : 0: use as it is (default), 1: count as 1, 2: phred score
	               (a)	
		_/     sp       \_
	[    5'  ]--------------[   3'      ]
               _(b)_          _(c)_
"
	if [ $# -lt 3 ];then echo "$usage"; return; fi
	local S=${3/1/"-s"}; S=${S/2/"-S"};
	local tmpd=`mymktempd`; 
	#local tmpd=tmpd; rm -rf $tmpd; mkdir -p $tmpd
	mycat $1 | bed.enc - > $tmpd/i
	bed.5p $tmpd/i > $tmpd/i.5
	bed.3p $tmpd/i > $tmpd/i.3
	mycat $2 | bed.score - ${4:-"0"} \
		| awk -v OFS="\t" -v O=$tmpd/r '{ n=1; if($10 > 1){ n=2;} print $0 >> O"."n; }'

	## count spliced
	bed.intron $tmpd/r.2 \
		| intersectBed -a $tmpd/i -b stdin -wa -wb -f 1 -F 1 $S \
		| cut -f4,11 | stat.sum - | sort -k 1,1  > $tmpd/a 

	## count unspliced
	intersectBed -a $tmpd/i.5 -b $tmpd/r.1 -wa -wb $S \
		| cut -f4,11 | stat.sum - | sort -k 1,1 > $tmpd/b

	intersectBed -a $tmpd/i.3 -b $tmpd/r.1 -wa -wb $S \
		| cut -f4,11 | stat.sum - | sort -k 1,1 > $tmpd/c
	
	join -a 1 -a 2 -e 0 -o 0,1.2,2.2 $tmpd/a $tmpd/b \
	| join -a 1 -a 2 -e0 -o 0,1.2,1.3,2.2 - $tmpd/c \
	| awk '{v=$3+$4; if(v==0){ v="inf";}else{ v=2*$2/v;}  print $1,$2,$3,$4,v;}' \
	| tr "@ " "\t"	

	rm -rf $tmpd;
}
asp.spi.test(){
echo \
"c	100	200	intron1	0	+
c	100	199	intron2	0	+" > tmp.intron
echo \
"c	90	210	r1	10	+	90	210	0,0,0	2	10,10	0,110
c	189	199	r2	1	+	190	200	0,0,0	1	10	0
c	190	200	r3	1	+	190	200	0,0,0	1	10	0
c	191	201	r4	1	+	190	200	0,0,0	1	10	0
c	199	209	r5	1	+	190	200	0,0,0	1	10	0
c	200	210	r6	1	+	190	200	0,0,0	1	10	0
c	201	211	r7	1	+	190	200	0,0,0	1	10	0" > tmp.r
head tmp.*
echo "==> result"
	asp.spi tmp.intron tmp.r 1 
	rm tmp.intron tmp.r
}

asp.3ss(){ 
usage=" 
 FUNCT: calculate 3SS statistics for each intron
 OUTPU: intron.bed6 a,b 1 - b/a
 USAGE: $FUNCNAME <intron.bed6> <read.bed12> <window> <strand> [<count>]
	<window>: windowsize for a and b boundaries 
	<strand>: 0, 1(same strand), 2(opposite strand) 
	<count> : 0: use as it is (default), 1: count as 1, 2: phred score


# >>>>> : read
# >>>>>----->>>>> : splicing read
# [     ] : exon
# ------- : intron
# |-----| : predefined window (e.g., 25bp)

# 3SS = 1- #(a)/#(b), [-inf, 1]

                         |--(a)--|--(b)--|
 [  5' exon    ]-----------------[   3' exon     ]
                     >>>>>      >>>>>           (support a)
                                 >>>>>   >>>>>  (support b)
          >>>>>>----------------->>>>>>         (support b)
             >>>>>>-------------->>>>>>         (support b or not)
          >>>>>>--------------------->>>>>>     (support none)

"
	if [ $# -lt 4 ];then echo "$usage"; return; fi
	local W=$3;
	local S=${4/1/"-s"}; S=${S/2/"-S"};

	local tmpd=`mymktempd`;
	mycat $1 | bed.enc - > $tmpd/i 
	bed.3p $tmpd/i | bed.flank - $(( $W -1 )) 0 1  > $tmpd/i.a
	bed.3p $tmpd/i | bed.flank - -1 $W 1 > $tmpd/i.b 
	mycat $2 | bed.score - ${5:-"0"} \
		| awk -v OFS="\t" -v O=$tmpd/r '{ n=1; if($10 > 1){ n=2;} print $0 >> O"."n; }'

	intersectBed -a $tmpd/i.a -b $tmpd/r.1 $S -wa -wb \
		| cut -f4,11 | stat.sum - | sort -k1,1 > $tmpd/a 	
	intersectBed -a $tmpd/i.b -b $tmpd/r.1 $S -wa -wb \
		| awk '$6=="+" && $8>=$2 || $6=="-" && $9<=$3' \
		| cut -f4,11 | stat.sum - | sort -k1,1 > $tmpd/b 	

	bed.intron $tmpd/r.2 \
		| intersectBed -a $tmpd/i -b stdin $S -wa -wb \
		| awk '$6=="+" && $9==$3 || $6=="-" && $8==$2' \
		| cut -f 4,11 | stat.sum - | sort -k1,1 > $tmpd/c

	join -a 1 -a 2 -e 0 -o 0,1.2,2.2 $tmpd/a $tmpd/b \
	| join -a 1 -a 2 -e 0 -o 0,1.2,1.3,2.2 - $tmpd/c \
	| awk '{v=$3+4; if(v==0){ v="-inf";}else{ v=1 - $2/v;}  print $1,$2","$3+$4,v;}' \
	| tr "@ " "\t"
	
	rm -rf $tmpd
}
asp.3ss.test(){ 
echo \
"c	100	200	intron1	0	+
c	100	199	intron2	0	+" > tmp.intron
echo \
"c	90	210	r1	1	+	90	210	0,0,0	2	10,10	0,110
c	80	210	r2	1	+	80	210	0,0,0	2	10,10	0,120
c	189	199	r3	1	+	190	200	0,0,0	1	10	0
c	190	200	r4	1	+	190	200	0,0,0	1	10	0
c	191	201	r5	1	+	190	200	0,0,0	1	10	0
c	199	209	r6	1	+	190	200	0,0,0	1	10	0
c	200	210	r7	1	+	190	200	0,0,0	1	10	0
c	201	211	r8	1	+	190	200	0,0,0	1	10	0" > tmp.r
head tmp.*
echo "==> result"
	asp.3ss tmp.intron tmp.r 10 1 
	rm tmp.intron tmp.r
}
asp.3ss_lessfp(){ 
usage="
FUNNC: calculate 3ss with less false positives caused by:
	1. introns w/o splicing reads 
	2. wrong 5/3 \' junctions
	3. overlapping reads at the 3\' end 

USAGE: $FUNCNAME <intron.bed6> <read.bed12> <window> <strand>
	<window>: windowsize for a and b boundaries 
	<strand>: 0, 1(same strand), 2(opposite strand) 

                     |---a---|---b---|
	[    ]---------------[             ]
	   \                 /       ( wrong 5\' splicing )
	     \                   /   ( wrong 3\' splicing )
	OUTPUT: intron.bed6 + a + b 
"
	if [ $# -lt 4 ];then echo "$usage"; return; fi
	local W=$3;
	local S=${4/1/"-s"}; S=${S/2/"-S"};

	tmpd=`mymktempd`;
	#tmpd=tmpd; mkdir -p $tmpd; 
	mycat $1 | awk -v OFS="\t" '{ $4=$1"@"$2"@"$3"@"$4"@"$5"@"$6;} 1' > $tmpd/a
	mycat $2 > $tmpd/b
	awk '$10==1' $tmpd/b | cut -f1-6 > $tmpd/b.1
	
	## junctio counts
	bed.intron $tmpd/b | intersectBed -a $tmpd/a -b stdin -f 1 -F 1 $S -wa -wb \
	| awk -v OFS="\t" '{ print $4,$11;}' | stat.sum - > $tmpd/c.j

	cat $tmpd/c.j | tr "@" "\t" | awk -v OFS="\t" '{ $4=$1"@"$2"@"$3"@"$4"@"$5"@"$6;} 1' \
	| bed.3p - | bed.flank - $(( $W -1))  $W 1 \
	| intersectBed -a stdin -b $tmpd/b.1 -wa -wb $S \
	| perl -e 'use strict; my $W='$W'; my %res=();
		while(<STDIN>){ chomp; my @a=split/\t/,$_;
			## bed7 + bed6
			my ($s, $e) = ($a[1],$a[2]);
			my ($rs, $re) = ($a[8],$a[9]);

			my $x="x";
			if($s + $W < $rs ){ $x="b";		
			}elsif( $s + $W >= $re ){ $x="a"; }
			if($a[5] eq "-"){
				if($x eq "a"){ $x = "b";}
				elsif($x eq "b" ){ $x = "a";}
			}	
			$res{$a[3]}{"j"} = $a[6];
			$res{$a[3]}{$x} += $a[11];
		}
		foreach my $k (keys %res){
			my $a=defined $res{$k}{"a"} ? $res{$k}{"a"} : 0;
			my $b=defined $res{$k}{"b"} ? $res{$k}{"b"} : 0;
			my $c=defined $res{$k}{"j"} ? $res{$k}{"j"} : 0;
			print join("\t", split /@/,$k),"\t",$a,"\t",$b+$c,"\n";
		}
	'
	rm -rf $tmpd
}
aso.cosi(){
	## BED: 2-3: intron boundary, 5th: exon coordinates
	##      /    a     \         /   b     \ 
	##               __c__     __d__       
	##      /               e              \
	##     ]------------[       ]-----------[
    BED=$1; BAM=$2;
    chroms=( `cut -f1 $BED | sort -u` )
    TMP=`make_temp`
    for CHROM in ${chroms[@]}
    do
        echo " $CHROM .." >&2
        awk -v CHROM=$CHROM '$1==CHROM' $BED > $TMP
        samtools view -bq 255 $BAM $CHROM | bamToBed -bed12 \
		| intersectBed -a $TMP -b stdin -wa -wb \
		| awk -v OFS="\t" '{ split($17,l,",");split($18,s,",");split($5,ex,",");
			is=$2-1;ie=$3; es=ex[1];ee=ex[2]-1; 
			a=0;b=0;c=0;d=0;e=0;
			for(i=1;i<=$16;i++){ ## count unsplicing events
					rs=$8+s[i]; re=rs+l[i]-1;
					if(rs < es && re > es){ c=1;
					}else if(rs < ee && re > ee){ d=1;}
			}
			for(i=2;i<=$16;i++){ ## count splicing events
					rs=$8+s[i-1]+l[i-1]-1; re=$8+s[i];
					if(rs==is && re== es){ a=1;
					}else if(rs==ee && re == ie){ b=1;
					}else if(rs ==is && re == ie){ e=1;}
			}
			if(a+b+c+d+e > 0){
				print $1,$2,$3,$4,$5,$6,a,b,c,d,e;
			}
		}'  | groupBy -g 1,2,3,4,5,6 -c 7,8,9,10,11 -o sum,sum,sum,sum,sum
	done
}
gen_bed(){
	echo "hi" | awk -v OFS="\t" -v s=$1 -v e=$2 -v l=$3 '{
		for(i=s; i< e; i+=1){
			print "c",i,i+l,"p",1,"+";	
		}	
		for(i=s; i< e; i+=1){
			print "c",i,i+l,"n",1,"-";	
		}	
	}'

}
test__3ss(){
echo \
"c	100	200	intron1	0	+
c	100	200	intron1	0	-" > tmpa
echo \
"c	100	200	intron1	0	+	25	25
c	100	200	intron1	0	-	25	25" > exp

gen_bed 0 250 50 | 3ss tmpa - 25 -Sm > obs
check obs exp
rm -rf tmpa exp obs
}
#test__3ss

test__count_us(){
fig="
Count unsplicing for exons:
            100         200
	    [          ]    : exon
       ======      ======   : unsplicing
            ===== ======    : not unsplicing
"
echo \
"chr1	100	200	e	0	+" > a
echo \
"chr1	50	101	r1	1	+
chr1	150	201	r2	10	+
chr1	100	150	r3	100	+
chr1	150	200	r4	1000	+" > b
echo \
"chr1	100	200	e	0	+	11" > exp
count_us a b > obs
check exp obs
rm -f a b exp obs
}
#test__count_us

bed12_to_jc(){
	bed12_to_intron $1 \
	| awk -v OFS="\t" '{ $2=$2-1;$3=$3+1;} 1' \
	| bed_sum - 
}
bede_to_bed6(){
## copy 5th,7th to the end to the name field, set 5th to be 0
	cat $1 | perl -ne 'chomp; my @a=split/\t/,$_;
	print join("\t",@a[0..3]),"@",$a[4],"@",join("@",@a[6..$#a]),"\t0\t",$a[5],"\n"; '
}
bed6_to_bede(){
	cat $1 | perl -ne 'chomp; my @a=split/\t/,$_; my @b=split/@/,$a[3]; $a[3] = shift @b; $a[4] = shift @b; 
	print join("\t",@a),"\t",join("\t",@b),"\n";	'
}

filter_intron_skip(){
	awk '$2-1==$8 && $3+1==$9'
}


count_is(){
	opt=${3:-""};
	intersectBed -a ${1/-/stdin} -b ${2/-/stdin} -wa -wb $opt \
	| awk -v OFS="\t" '$2==$8+1 && $3==$9-1{ print $1,$2,$3,$4"|"$5,$11,$6;}' \
	| bed_sum - \
	| awk -v OFS="\t" '{ split($4,a,"|");print $1,$2,$3,a[1],a[2],$6,$5;}' 
}
test__count_is(){
echo "c	100	200	.	1000	+" > a
echo "c	99	201	.	1	+
c	99	200	.	10	+
c	99	201	.	2	-" >b
	count_is a b -s > obs
	count_is a b -S >> obs
	count_is a b >> obs
echo "c	100	200	.	1000	+	1
c	100	200	.	1000	+	2
c	100	200	.	1000	+	3" > exp
cat obs
	check exp obs
	rm -f obs exp a b

}
#test__count_is;



count_jei(){
usage="
usage: $FUNCNAME <target.bed6> <read.bed12>
function: count exclusion and inclusion of junction reads
<read.bed12> : use modify_bed .. to manage score the column


output: bed6@ count_exc count_inc
       /                  \     exclusion
                 /         \    exclusion
            \          /        inclusion
  -----------[        ]----------
"
if [ $# -lt 2 ]; then echo "$usage"; return; fi
        bed12_to_junction $2 \
        | intersectBed -a $1 -b stdin -wa -wb \
        | awk -v OFS="@" '{
                e=0;i=0; ## exclusion and inclusion
                if( $8 == $3 -1  || $9 == $2 + 1){ i=$11;}
                else{ e=$11;}
                print $1,$2,$3,$4,$5,$6"\t"e"\t"i;
        }' | sort -k1,1 | groupBy -g 1 -c 2,3 -o sum,sum
}


make_temp(){
	mktemp 2>/dev/null || mktemp -t $0
}


fdr_thre(){
	F=$1; T=$3;
	cmd='
		tt = read.table("stdin",header=F);
		x= tt[,C];
		ix = p.adjust(x,method="fdr") <= T;
		cat(paste("res=",max(x[ix]),"\n",sep=""));
	'
	cmd=${cmd/C/$C};
	cmd=${cmd/T/$T};

	tmp=`make_temp`
	echo "$cmd" > $tmp
	cat $F | R --no-save -f $tmp | perl -ne 'chomp;if($_=~/res=([\d|\.]+)/){ print $1,"\n";}'
}

count_a53ss_jc(){
	##          /         a         \        
	##               /    b         \
	##     [   |   ]----------------[       ]
    BED=$1; JBED=$2;
	intersectBed -a $BED -b $JBED -wa -wb -f 1 -r \
	| awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$10; }'
}

count_a53ss(){
	##          /         a         \        
	##               /    b         \
	##     [   |   ]----------------[       ]
    BED=$1; BAM=$2;
    chroms=( `cut -f1 $BED | sort -u` )
    TMP=`make_temp`
    for CHROM in ${chroms[@]}
    do
        echo " $CHROM .." >&2
        awk -v CHROM=$CHROM '$1==CHROM' $BED > $TMP
        samtools view -bq 255 $BAM $CHROM | bamToBed -bed12 \
		| intersectBed -a $TMP -b stdin -wa -wb \
		| awk -v OFS="\t" '{ alen=split($7,a,",");blen=split($8,b,",");split($19,sizes,",");split($20,starts,",");
			for(i=1 ; i<= alen; i++){ ai=a[i]-1;
			for(j=1 ; j<= blen; j++){ bj=b[j]+1;
			for(k=2 ; k<= $18; k++){
				s=$10 + starts[k-1] + sizes[k-1]-1;		
				e=$10 + starts[k] + 1;		
				if(ai == s && bj == e ){
					print $1,$2,$3,$4,ai "," bj,$6,1;
				}
			}}}
		}' | sort -k1 -k2,3n | groupBy -g 1,2,3,4,5,6 -c 7 -o sum | awk -v OFS="\t" '{ print $1,$2,$3,$4,$5,$6,$7;}'
	done
}
#head -n 100 a53ss.bed  > tmp
#count_a53ss tmp ../Tophat/Wt1/accepted_hits.bam -wa -wb


count_bed(){
    BED=$1; BAM=$2;
    chroms=( `cut -f1 $BED | sort -u` )
    TMP=`make_temp`
    for CHROM in ${chroms[@]}
    do
        echo " $CHROM .." >&2
        awk -v CHROM=$CHROM '$1==CHROM' $BED > $TMP
        samtools view -bq 255 $BAM $CHROM | bamToBed -bed12 \
		| intersectBed -a $TMP -b stdin -wa -wb \
		| awk -v OFS="\t" '{ split($17,l,",");split($18,s,","); 
			hit=0;
			if($16==1){ hit=1;
			}else{
				for(i=1;i<=$16;i++){
					if(i==1 && $8+s[i]+l[i] == $3){  ##  [  _]/
						hit=1;
					}else if(i==$16 && $8+s[i] == $2){ ##  \[_  ]
						hit=1;
					}else if($2 == $8+s[i] && $3==$8+s[i]+l[i]){ ## \[__]/
						hit=1;
					}
				}
			}
			if(hit){ print $1,$2,$3,$4,$5,$6,1;}
		}'  | groupBy -g 1,2,3,4,5,6 -c 7 -o sum
	done
}

count_nj(){
    BED=$1; BAM=$2;
    chroms=( `cut -f1 $BED | sort -u` )
    TMP=`make_temp`
    for CHROM in ${chroms[@]}
    do
        echo " $CHROM .." >&2
        awk -v CHROM=$CHROM '$1==CHROM' $BED > $TMP
        samtools view -bq 255 $BAM $CHROM | bamToBed -split \
        | intersectBed -a $TMP -b stdin -wa -wb \
        | awk '{ if($2 > $8 || $3 < $9){ print $1,$2,$3,$4,$5,$6;}}' \
        | uniq -c | awk -v FS=" " -v OFS="\t" '{ print $2,$3,$4,$5,$6,$7,$1;}'
    done
}

intron_excinc(){
    NJ=$1; JC=$2;
    intersectBed -a $NJ -b $JC -wa -wb  \
    | awk -v OFS="\t" '{ if($2-1==$9 && $3+1==$10){ print $1,$2,$3,$4,$5,$6,$11,$7;}}'
}


make_retainable_introns(){
## [    ]—-intron--[   ]
## [                   ]
	intersectBed -a $1 -b $2 -wa -wb \
	| awk -v OFS="\t" '{ split($5,a,",");
		if($8 == a[1] && $9 == a[2]){ ## match boundary
			print $1,$2,$3,$4,0,$6;
		}
	}' | sort -u  
	## unique retained intron candidates
}







## changed from bed12_to_exon
bed12_to_exonEvents(){
	## [s1  ]e1----[s   ]e----[s2   ]e2
	awk -v OFS="\t" '{
		## take introns
		split($11,sizes,",");
		split($12,starts,",");
		for(i=1;i<= $10;i++){
			s = $2 + starts[i];
			e = $2 + starts[i]+sizes[i];	
			s1 = -1; e1 = -1;
			s2 = -1; e2 = -1;
			if( i > 1){ s1 = $2 + starts[i-1]; e1 = s1 + sizes[i-1]; }
			if( i < $10){ s2 = $2 + starts[i+1]; e2 = s2 + sizes[i+1]; }
			split($4,gene,"::");
			## gene base
			print $1,s,e,gene[1],s1 "," e1 "," s2 "," e2,$6;
		}	
	}' | sort -u
}




each_chrom(){
	## lambda function accepts chrom and size parameters
	## lambda(){ 
	##   # .. handle $1 $2 
	## } 
	chrom_size_file=$1; lambda_func=$2;
	tmp=( `cat $chrom_size_file` )
	for (( i=0; i< ${#tmp[@]}; i+=2 ))
	do
		chrom=${tmp[$i]}; size=${tmp[$i+1]};
		$lambda_func $chrom $size
	done
}


bam_to_junctions(){
	## split bam by chrom to reduce memory usage
	bam=$1; chromsize=$2; quality=$3;
	lambda(){
		samtools view -bq $quality $bam $1 | bamToBed -bed12 | awk '$10 > 1' | bed12_to_junction
	}
	each_chrom $chromsize lambda
}
#bam_to_junctions $1 $2 255
exon_excinc_junction(){
    ## input: exon and junction_count
    ## output: exon skipping and junction counts 
	## [   ]2----[5.1  5.2]---3[   ]
    EXON=$1;JUNCTION=$2;
	cat $EXON | awk -v OFS="\t" '{ 
		## take neighbor intron boundary
   		split($5,a,","); ## [a1  ]a2--[2   ]3--[a3  ]a4
		s=$2;e=$3;
		if(a[2] > 0){ s= a[2];} # intron start
		if(a[3] > 0){ e= a[3];} # intron end
		print $1,s, e, $4, $2 "," $3,$6; 
	}' | sort -u | intersectBed -a stdin -b $JUNCTION -wa -wb  | awk -v OFS="\t" '{
    	## [  ]e1----[e2   e3]---e4[   ]
    	split($5,a,","); e1=$2-1; e2=a[1]; e3=a[2]-1; e4=$3; # exon
		js = $8; je = $9-1; jc=$10; # junction 
        x=0;y=0;z=0;
        if(js==e1 && je==e2 ){ ## left splicing  
            x=jc;
        }else if(js==e3 && je==e4){ ## right splicing 
            y=jc;
        }else if(js==e1 && je==e4){ ## [  ]/  [  ] \[  ]skipping
            z=jc;
        }
    	if( x+y+z > 0){ 
			## output: bed6: [   2]--[5.1   5.2]---[3   ], counts:  7,8,9
			#print $1,$2,$3,$4,$5,$6,x,z,y; ## left splicing, skipping, right splicing
			print $1,$2,$3,$4,$5,$6,z,x+y;  ## exclusion, inclusion
		}
    }' | sort -k1,1 -k2,3n -k4,6 | groupBy -g 1,2,3,4,5,6 -c 7,8 -o sum,sum 
# 	| awk -v OFS="\t" '{
#    	split($5,a,","); # [   2]---[a1  a2]---[3  ]
#		if($8>0 && $7> 0){
#			print $1,$2,$3,$2,a[1],$6,$8,$7;
#		}
#		if($8>0 && $9>0){
#			print $1,$2,$3,a[2],$3,$6,$8,$9;
#		}
#	}' | sort -u
}

count_nonjunction_events(){
	#   --  --  ---  : left unspliced, within, right unspliced
    #   -----------  : within = 0, left unspliced == right unspliced 
	#    [       ]
	bed6=$1; bam=$2; chromsize=$3; quality=$4
	lambda(){
		samtools view -bq $quality $bam $1| bamToBed -bed12 \
		| intersectBed -a $bed6 -b stdin -wa -wb \
		| awk -v OFS="\t" '{ 
			L = 0; C = 0; R = 0;
			if($8 < $2 && $9-1 > $2){ L += 1;}  # left unspliced
			if($8 < $3-1 && $9 > $3){ R += 1;}  # right unspliced 
			if($8 > $2 && $9 < $3){ C += 1;}    # within 
			if( L + C + R > 0){
				print $1,$2,$3,$4,$5,$6,L,C,R;
			}
		}' | sort -k1,1 -k2,3n -k4,6 | groupBy -g 1,2,3,4,5,6 -c 7,8,9 -o sum,sum,sum  
	}
	each_chrom $chromsize lambda
}
#count_nonjunction_events $1 $2 $3 $4

count_intron_junction_events(){
	intron=$1; junction=$2;
	intersectBed -a $intron -b $junction -wa -wb  | awk -v OFS="\t" '{
		##     /       \
		## [   ]-------[    ]
		if($2-1 == $8 && $3+1 == $9){ 
			print $1,$2,$3,$4,$5,$6,$10;
		}
	}'
}
		
count_intron_events(){
	intron=$1;bam=$2;junction=$3;chromsize=$4;quality=$5;
	tmp1=`make_temp`
	tmp2=`make_temp`
	count_intron_junction_events  $intron $junction  > $tmp1
	count_nonjunction_events $intron $bam $chromsize $quality > $tmp2
	intersectBed -a $tmp2 -b $tmp1 -wa -wb -f 1 -r -s | awk -v OFS="\t" '{
		if($4 == $13){ ## shared by different genes
			## exclusion and inclusion
			print $1,$2,$3,$4,$5,$6,$16,$7+$9;	
		}
	}'
}
#count_intron_events Events/rintrons.bed ../Tophat/Wt1/accepted_hits.bam Events/Wt1/jc.bed Data/chrom.size 255

quote(){ 
	perl -ne 'chomp; my @a = map{ "\"$_\"" } split /,/,$_; print join ",",@a; '; 
}

testI(){
	## INPUT: comma separated control and treatment EI file ( bed6 + exclusion + inclusion)
	## OUTPUT: bed6 + logFC + pvalue 
	a=`echo $1 | quote`
	b=`echo $2 | quote`
	rcmd='
	fa=c(FILEA)
	fb=c(FILEB)
	#fa=c("Events/Wt1/exons_jc.bed","Events/Wt2/exons_jc.bed")	
	#fb=c("Events/C41/exons_jc.bed","Events/C42/exons_jc.bed")	
	out="OUT"
	group=factor(c(rep(1,length(fa)),rep(2,length(fb))));
	D=NULL;
	i=1;
	for( f in c(fa,fb)){
		tt=read.table(f,header=F);
		colnames(tt)=c("chr","start","end","name","score","strand",paste(i,c("exc","inc"),sep="."))
		if(is.null(D)){ D=tt;
		}else{ D=merge(D,tt,by=1:6,all=T); }
		i=i+1;
	}
	D[is.na(D)]=0;
	ix = apply(D[,grep("inc",colnames(D))],1, min) > 2;
	d=D[ix,1:6];
	m=D[ix,grep("inc",colnames(D))];
	o=order(d$name);
	## make exclusion event from the gene sum
	d=d[o,]; m=m[o,]; s=apply(m,2,function(x){ ave(x,d$name,FUN=sum)});
	s = s-m; colnames(s)=paste(1:ncol(s),"exc");
	m = cbind(s, m);


	## test
    library(edgeR)
    #j=ncol(m)/2; #y=DGEList(counts=m[,1:j]+m[,(j+1):ncol(m)],group=group)
    y=DGEList(counts=m,group=rep(group,2))
    y=calcNormFactors(y);

    event.this=factor(rep(1:2,each=length(group)));
    group.this=factor(rep(group,2));
    H1 <- model.matrix(~ event.this + group.this + event.this:group.this, data=y$samples )
    H0 <- model.matrix(~ event.this + group.this )
    coef <- (ncol(H0)+1):ncol(H1)
    #y=estimateCommonDisp(y)
    #y=estimateTagwiseDisp(y, trend="movingave")
    y = estimateGLMCommonDisp(y,H1);
    y = estimateGLMTrendedDisp(y,H1);
    y = estimateGLMTagwiseDisp(y,H1);

    fit=glmFit(m, H1, y$tagwise.dispersion,offset=0,prior.count=0)
    llh=glmLRT(fit,coef=coef)

    ex.h0=apply( m[,group.this == 1 & event.this == 1], 1, sum);
    in.h0=apply( m[,group.this == 1 & event.this == 2], 1, sum);

    res=data.frame(d[,1:6], logIR=log( in.h0/ ex.h0), logFC=llh$table$logFC, pval=llh$table$PValue)
    ## chrom start end logFC pval
    write.table(res, out, col.names=T,row.names=F,sep="\t",quote=F);
	'
	tmp=`mktemp`
	rcmd=${rcmd/FILEA/$a}
	rcmd=${rcmd/FILEB/$b}
	rcmd=${rcmd/OUT/$tmp}
	echo "$rcmd" | R --no-save >&2
	cat $tmp
}

test_exin_fisher(){
	## input : bed6 + exclusion + inclusion counts
	intersectBed -a $1 -b $2 -wa -wb -f 1 -r -s \
	| awk -v OFS="@" '{ print $1,$2,$3,$4,$5,$6"\t"$7"\t"$8"\t"$15"\t"$16;}' \
	| fisher_test -
}
test_exin_edger(){
	## INPUT: comma separated control and treatment EI file ( bed6 + exclusion + inclusion)
	## OUTPUT: bed6 + logFC + pvalue 
	a=`echo $1 | quote`
	b=`echo $2 | quote`

	rcmd='
	fa=c(FILEA)
	fb=c(FILEB)
	#fa=c("Events/Wt1/exons_jc.bed","Events/Wt2/exons_jc.bed")	
	#fb=c("Events/C41/exons_jc.bed","Events/C42/exons_jc.bed")	
	out="OUT"

	group=factor(c(rep(1,length(fa)),rep(2,length(fb))));
	D=NULL;
	i=1;
	for( f in c(fa,fb)){
		tt=read.table(f,header=F);
		colnames(tt)=c("chr","start","end","name","score","strand",paste(i,c("exc","inc"),sep="."))
		if(is.null(D)){ D=tt;
		}else{ D=merge(D,tt,by=1:6,all=T); }
		i=i+1;
	}
	D[is.na(D)]=0;
	ix=apply(D[,7:ncol(D)], 1, min) > 0 & apply(D[,7:ncol(D)],1,max) > 10
	d=D[ix,]
	#d[d==0]=0.5
	m=cbind(d[,grep("exc",colnames(d))],d[,grep("inc",colnames(d))])

	library(edgeR)
	#j=ncol(m)/2; #y=DGEList(counts=m[,1:j]+m[,(j+1):ncol(m)],group=group)
	y=DGEList(counts=m,group=rep(group,2))
	y=calcNormFactors(y);

	event.this=factor(rep(1:2,each=length(group)));
	group.this=factor(rep(group,2));
	H1 <- model.matrix(~ event.this + group.this + event.this:group.this, data=y$samples )
	H0 <- model.matrix(~ event.this + group.this )
	coef <- (ncol(H0)+1):ncol(H1)
	#y=estimateCommonDisp(y)
	#y=estimateTagwiseDisp(y, trend="movingave")
	y = estimateGLMCommonDisp(y,H1);
	y = estimateGLMTrendedDisp(y,H1);
	y = estimateGLMTagwiseDisp(y,H1);

	fit=glmFit(m, H1, y$tagwise.dispersion,offset=0,prior.count=0)
	llh=glmLRT(fit,coef=coef)

	ex.h0=apply( m[,group.this == 1 & event.this == 1], 1, sum);
	in.h0=apply( m[,group.this == 1 & event.this == 2], 1, sum);

	res=data.frame(d[,1:6], logIR=log( in.h0/ ex.h0), logFC=llh$table$logFC, pval=llh$table$PValue)
	## chrom start end logFC pval
	write.table(res, out, col.names=T,row.names=F,sep="\t",quote=F);
	'

	tmp=`mktemp`
	rcmd=${rcmd/FILEA/$a}
	rcmd=${rcmd/FILEB/$b}
	rcmd=${rcmd/OUT/$tmp}
	echo "$rcmd" | R --no-save >&2
	cat $tmp
}

#data='chr7	99647389	99662661	ENSG00000166529::ENST00000543588	0	+	99647389	99662661	0,0,0	6	75,133,435,193,192,956	0,1711,7204,7931,14021,14316
#chr7	99647396	99662661	ENSG00000166529::ENST00000456748	0	+	99647396	99662661	0,0,0	5	68,495,193,192,956	0,7137,7924,14014,14309
#chr7	99654523	99662660	ENSG00000166529::ENST00000379635	0	+	99654523	99662660	0,0,0	4	349,81,193,1250	0,424,797,6887
#chr7	99647416	99655464	ENSG00000166529::ENST00000438937	0	+	99647416	99655464	0,0,0	4	48,133,495,144	0,1684,7117,7904
#chr7	99654533	99661878	ENSG00000166529::ENST00000477297	0	+	99654533	99661878	0,0,0	3	495,193,173	0,787,7172
#chr7	99647396	99662659	ENSG00000166529::ENST00000292450	0	+	99647396	99662659	0,0,0	4	68,495,193,1249	0,7137,7924,14014';

bed12_to_exon(){
	awk -v OFS="\t" '{
		## take introns
		split($11,sizes,","); split($12,starts,",");
		for(i=1;i<= $10;i++){ ## 1-base index not perl
			ord=i; if($5 eq "-"){ ord = $10-ord+1; }
			s = $2 + starts[i]; e = $2 + starts[i]+sizes[i];	
			print $1,s,e,$4,0,$6;
		}	
	}' $1
}

bed12_to_a53ss(){
	## event id is composed of start,end,gene
	## 
	##     /  splicing event    \
	## [a   b]--------------------[ c   ]d : numbers represent columns
	## event_id:  gene:a-d:c or gene:a-d:b
	awk -v OFS="\t" '{ split($11,sizes,",");split($12,starts,",");split($4,gt,"::");
		for(i=2;i<=$10;i++){
			a=$2+starts[i-1]; b=a+sizes[i-1]-1;
			c=$2+starts[i]; d=c+sizes[i];
			print $1 "," gt[1] ":" a "-" d ":" b "," $6,b,c+1;
			print $1 "," gt[1] ":" a "-" d ":" c "," $6,b,c+1;
		}
	}'  | sort -uk 1 | groupBy -g 1 -c 2,3 -o collapse,collapse \
	| perl -ne ' chomp; my @a = split /\t/,$_;
		my @starts = split /,/,$a[1];
		next if scalar @starts < 2;
		my ($chrom,$eid,$strand) = split /,/,$a[0];
		my ($event_gene, $event_range, $event_pivot) = split /:/,$eid;
		my @ends = split /,/,$a[2];
		if($starts[0] == $event_pivot){ ## [   ]-----[  ][  ]
			my $type = $strand eq "+" ? "A3" : "A5";
			my @ix = sort { $ends[$a] cmp $ends[$b] } 0 .. $#ends;
			for my $i (0 .. $#ix){
				print join "\t",($chrom,$starts[$ix[$i]],$ends[$ix[$i]],$eid,$type.":".$i,$strand),"\n";	
			}
		}else{ ## [   ][   ]-------[   ]
			my $type = $strand eq "+" ? "A5" : "A3";
			my @ix = sort { $starts[$b] cmp $starts[$a] } 0 .. $#starts;
			for my $i (0 .. $#ix){
				print join "\t",($chrom,$starts[$ix[$i]],$ends[$ix[$i]],$eid,$type.":".$i,$strand),"\n";	
			}
		} 
	'
}


bed12_to_rintrons(){
	local tmp=`make_temp`	
	local tmp1=`make_temp`
	local tmp2=`make_temp`
	cat > $tmp
	cat $tmp | bed12_to_exon > $tmp1
	cat $tmp | bed12_to_dumbbells > $tmp2
	intersectBed -a $tmp2 -b $tmp1 -r -f 1  -s -wa 
}

bed12_to_a53ss1(){
# 5' splicing      [$2                 a1]----[a2 $3]
#                  [    ]                     [     ]
# 3' splicing      [$2 a1]---------[a2           $3 ]
	bed12_to_2exons | groupBy -g 1,2,3,4,6 -c 5,5 -o collapse,count | awk '$7 > 1' #\
#	 | intersectBed -a stdin -b $tmp1 -wa -f 1 -r -s -v \
#	 | perl -ne ' chomp; my ($c,$s,$e,$n,$st,$tmp) = split /\t/,$_;
#		my %h5=();
#		my %h3=();
#		my %ri_filter=();
#		foreach my $e (split /,/,$tmp){
#			my ($es, $ee) = split /:/,$e;
#			$h5{$ee}{$es} = 5; # 5 prim splicing [   |   ]------[   ]
#			$h3{$es}{$ee} = 3; # 3 prim splicing [   ]--------[  |  ]
#		}	
#		foreach $k (keys %h5){
#			my @a=keys %{$h5{$k}};
#			if(scalar @a > 1){
#				my $st1 = ($st eq "+" ? 5 : 3);
#				print "$c\t$s\t$e\t$n\t$st1\t$st\t";
#				print join( ",",@a),"\t",$k,"\n";
#			}
#		}
#		foreach $k (keys %h3){
#			my @a=keys %{$h3{$k}};
#			if(scalar @a > 1){
#				my $st1 = ($st eq "+" ?  3 : 5);
#				print "$c\t$s\t$e\t$n\t$st1\t$st\t";
#				print $k,"\t",join( ",",@a),"\n";
#			}
#		}
#	' 
} 
## test
#echo -e "$data" | bed12_to_a53ss


