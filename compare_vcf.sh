#!/usr/bin/env bash
file1_path=`dirname $1`
file1_name=`basename $1`
file2_path=`dirname $2`
file2_name=`basename $2`

cat $1 | vcf-sort | bgzip -c > $1.gz ; tabix -p vcf $1.gz
cat $2 | vcf-sort | bgzip -c > $2.gz ; tabix -p vcf $2.gz

vcf-isec -c $1.gz $2.gz > $file1_path/compare_vcf.uniqueTo.$file1_name.versus.$file2_name
vcf-isec -c $2.gz $1.gz > $file2_path/compare_vcf.uniqueTo.$file2_name.versus.$file1_name
vcf-isec -n =2 $1.gz $2.gz > $file1_path/compare_vcf.commonTo.$file1_name.and.$file2_name

rm $1.gz $1.gz.tbi $2.gz $2.gz.tbi
