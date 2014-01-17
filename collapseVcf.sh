#!/usr/bin/env bash
file1_path=`dirname $1`
file1_name=`basename $1`
file2_path=`dirname $2`
file2_name=`basename $2`

cat $1 | vcf-sort | bgzip -c > $1.gz ; tabix -p vcf $1.gz
cat $2 | vcf-sort | bgzip -c > $2.gz ; tabix -p vcf $2.gz

vcf-isec -c $1.gz $2.gz > $file1_path/collapseVcf.$file1_name.and.$file2_name
vcf-isec -c $2.gz $1.gz >> $file1_path/collapseVcf.$file1_name.and.$file2_name
vcf-isec -n =2 $1.gz $2.gz >> $file1_path/collapseVcf.$file1_name.and.$file2_name

rm $1.gz $1.gz.tbi $2.gz $2.gz.tbi
