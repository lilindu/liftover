#!/bin/bash
# 目的是基于 reference，构建 assembly 的 chain 文件，便于后续进行坐标转换等分析

# 参考：https://iamphioxus.org/tag/genome-assembly/

# History
#   First_release	Wen	2022-05-04 citun@163.com
#	    以 reference genome 为基准，分析 I II III MT
#	    示例见 21-11-20 实验一 及 22-05-04 实验一
#   Second_release    Wen 2023-04-04
#       set the outdir automatically
#   Third_release   Wen 2023-04-17
#       Use all the reference contigs
#   Fourth_release   Wen 2024-08-17
#       Could specify a reference

assembly=$1
ID=$2
ref_dir=$3 # /data/a/liwen/pombase/version_210306/spombe_liftover
reference=$4 # /data/a/liwen/pombase/version_210306/Spombe_ref_genome_210306.fa
# set a parameter ref_id, default is Ref
ref_id=${5:-Ref}

outdir=${ref_dir}/${ID}To${ref_id}_chain

if [ ! -d $outdir ]
	then mkdir -p $outdir
fi

log=${outdir}/${ID}To${ref_id}_chain.log

# blat
mkdir ${outdir}/ref_psl

for i in `ls ${ref_dir}/*_split/*.fa`
do
    prefix=`basename $i`
    prefix=${prefix%.split.fa}
    blat $assembly $i -t=dna -q=dna \
        -tileSize=12 -fastMap -minIdentity=95 -noHead -minScore=100 \
        ${outdir}/ref_psl/${prefix}.psl >> $log 2>&1
done

# liftUp
mkdir ${outdir}/ref_liftup

for i in `ls ${ref_dir}/*_lift/*.lft`
do
    prefix=`basename $i`
    prefix=${prefix%.lft}
    liftUp -pslQ ${outdir}/ref_liftup/${prefix}.liftup.psl $i \
        warn ${outdir}/ref_psl/${prefix}.psl >> $log 2>&1
done

# axtChain
mkdir ${outdir}/ref_chain_raw

for i in `ls ${outdir}/ref_liftup/*.liftup.psl`
do
    prefix=`basename $i`
    prefix=${prefix%.liftup.psl}
    axtChain -linearGap=medium -faQ -faT -psl $i $assembly $reference \
        ${outdir}/ref_chain_raw/${prefix}.chain >> $log 2>&1
done

# combine
chainMergeSort ${outdir}/ref_chain_raw/*.chain | \
	chainSplit ${outdir}/ref_chain_split stdin

faSize $assembly -detailed > ${outdir}/${ID}.chr_length.txt

# alignment nets
mkdir ${outdir}/ref_net

for i in `ls ${outdir}/ref_chain_split/*.chain`
do
    prefix=`basename $i`
    prefix=${prefix%.chain}
    chainNet $i ${outdir}/${ID}.chr_length.txt ${ref_dir}/ref.chr_length.txt \
        ${outdir}/ref_net/${ID}_${prefix}.net ${outdir}/tmp_${prefix}.null >> $log 2>&1
done

# liftOver chain file
mkdir ${outdir}/ref_over

for i in `ls ${outdir}/ref_chain_split/*.chain`
do
    prefix=`basename $i`
    prefix=${prefix%.chain}
    netChainSubset ${outdir}/ref_net/${ID}_${prefix}.net \
        $i ${outdir}/ref_over/${ID}_${prefix}.chain >> $log 2>&1
done

cat ${outdir}/ref_over/*chain \
	> ${outdir}/${ID}To${ref_id}.over.chain