#!/bin/bash


mycat(){
        if [[ ${1##*.} = "gz" ]];then
                gunzip -dc $1;
        elif [[ ${1##*.} = "bz2" ]];then
                bunzip2 -dc $1;
        elif [[ ${1##*.} = "zip" ]];then
                unzip -p $1;
        else
                cat $1;
        fi
}

check(){ ## obtained from bamtools test
	if diff $1 $2; then
		echo ok
		return 1
	else
		echo fail
		return 0
	fi
}

rmempty(){
	for f in $@;do
		(( $(stat -c %s "$f") )) || rm "$f"
	done
}

make_tempdir(){
	mktemp -d 
}

run_R(){
usage="
usage: $FUNCNAME <Rscript> 
"
		
        cmd=$1;
        cmd=${cmd/stdout/tmp.out}
        echo "$cmd" > tmp.cmd
        R --no-save -q -f tmp.cmd &> tmp.log;
	if [ -f tmp.out ];then
		cat tmp.out
	else
		cat tmp.log
	fi
	rm -f tmp.cmd tmp.out tmp.log
}


###########################################################
# test 
###########################################################
test(){
	echo "too clean to test"
}


