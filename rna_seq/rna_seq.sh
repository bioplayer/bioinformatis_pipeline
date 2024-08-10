files=("Aspirin" "CTR" "NaAc")
data="/mnt/data/whn/wongLab/wux/rna_20240629/01.RawData"
ws="/mnt/data/whn/wongLab/wux/rna_20240629"
ref="/mnt/data/whn/reference/Mus_musculus"
thread=12

if [ 1 = 2 ];then
  mkdir -p ${ws}/1.fastqc
  cd ${ws}/1.fastqc || exit
  for i in {0..2}
  do
    for j in {1..3}
    do
      f1=${data}/${files[i]}${j}_1.fq.gz
      f2=${data}/${files[i]}${j}_2.fq.gz
      trim_galore --fastqc --paired -j ${thread} --basename ${files[i]}_${j} ${f1} ${f2}
    done
  done
fi;

if [ 1 = 1 ];then
  mkdir -p ${ws}/2.STAR
  cd ${ws}/2.STAR || exit
  for i in {0..2}
  do
    for j in {1..3}
    do
      STAR --genomeDir ${ref}/StarIndex/ --runThreadN ${thread} --readFilesIn ${ws}/1.fastqc/${files[i]}_${j}_val_1.fq.gz ${ws}/1.fastqc/${files[i]}_${j}_val_2.fq.gz --readFilesCommand "zcat" --outFileNamePrefix ${files[i]}_${j}_ --outSAMtype BAM SortedByCoordinate
      samtools view -bq 30 -o ${files[i]}_${j}.bam -@ ${thread} ${files[i]}_${j}_Aligned.sortedByCoord.out.bam
    done
  done
fi;

if [ 1 = 1 ];then
  mkdir -p ${ws}/3.DEG
  cd ${ws}/3.DEG
  for i in {0..2}
  do
    for j in {1..3}
    do
      featureCounts -p -a ${ref}/basic.gtf -o ${files[i]}_${j}_fc.txt -T ${thread} ${ws}/2.STAR/${files[i]}_${j}_Aligned.sortedByCoord.out.bam
      n_read=`awk '{n=n+$2}END{print n}' ${files[i]}_${j}_fc.txt.summary`
      echo ${n_read}
    done
  done
fi;
