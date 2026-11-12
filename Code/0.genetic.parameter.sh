#heterozygosity

cd ${bam_path}
for ind in `cat ${group_list}`
do
   ${samtools}/samtools mpileup -t AD,ADF,ADR,DP,SP -Bugf ${ref} ${ind}.RG.bam | ${bcftools}/bcftools call -c - \
    | ${bcftools}/vcfutils.pl vcf2fq -d 10 -D 100 > ${output_path}${ind}.fq
done 

for ind in `cat ${group_list}`
do
   seqtk seq -a -q20 -n N ${output_path}${ind}.fq > ${output_path}${ind}.f
done 


for ind in `cat ${group_list}`
do
  sed -i '/>/d' ${output_path} ${ind}.f 
done

for ind in `cat ${group_list}`
do
   cat ${output_path}${ind}.f | tr a-z A-Z > ${output_path}${ind}.fa
done

for ind in `cat ${group_list}`
do
   awk -v RS='A' 'END {print "A "--NR}' ${output_path}${ind}.fa >> ${output_path}${ind}.a
   awk -v RS='T' 'END {print "T "--NR}' ${output_path}${ind}.fa >> ${output_path}${ind}.a
   awk -v RS='C' 'END {print "C "--NR}' ${output_path}${ind}.fa >> ${output_path}${ind}.a
   awk -v RS='G' 'END {print "G "--NR}' ${output_path}${ind}.fa >> ${output_path}${ind}.a
   awk -v RS='K' 'END {print "K "--NR}' ${output_path}${ind}.fa >> ${output_path}${ind}.a
   awk -v RS='M' 'END {print "M "--NR}' ${output_path}${ind}.fa >> ${output_path}${ind}.a
   awk -v RS='R' 'END {print "R "--NR}' ${output_path}${ind}.fa >> ${output_path}${ind}.a
   awk -v RS='S' 'END {print "S "--NR}' ${output_path}${ind}.fa >> ${output_path}${ind}.a
   awk -v RS='W' 'END {print "W "--NR}' ${output_path}${ind}.fa >> ${output_path}${ind}.a
   awk -v RS='Y' 'END {print "Y "--NR}' ${output_path}${ind}.fa >> ${output_path}${ind}.a
   awk -v RS='N' 'END {print "N "--NR}' ${output_path}${ind}.fa >> ${output_path}${ind}.a 
done

for ind in `cat ${group_list}`
do
   cat ${output_path}${ind}.a | awk '{print $2}' | head -n 10 | awk '{ sum += $1; } END { print  sum}' >> ${output_path}${ind}.B
   cat ${output_path}${ind}.a | awk '{print $2}' | sed -n '5,10p' | awk '{ sum += $1; } END { print  sum}' >> ${output_path}${ind}.B
   awk '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS $i:$i}END{for(i=0;i++<NF;)print a[i]}' ${output_path}${ind}.B > ${output_path}${ind}.B2
   awk '(sum=$2/$1){print sum}' ${output_path}${ind}.B2 >  ${output_path}${ind}.B3
   rm ${output_path}${ind}.B
   rm ${output_path}${ind}.B2
done

#Genetic diversity

angsd -bam pop.bam.list -anc ref.fasta \
	-out pop -doSaf 1 -gl 1 -minInd `wc -l pop.bam.list` -minMapQ 30 -minQ 20 -nThreads 30
realSFS pop.saf.idx -maxIter 100 -P 8 > pop.sfs
realSFS saf2theta pop.saf.idx -outname pop -sfs pop.sfs
thetaStat do_stat pop.thetas.idx -win 50000 -step 10000 -outnames pop.window.theta

#gone2
gone2 -g 0 -t 20 -r 0.94 -s 100000 -o EU.pop pop.vcf

#tajimaD
vcftools --vcf pop.vcf --TajimaD 50000 --out pop
