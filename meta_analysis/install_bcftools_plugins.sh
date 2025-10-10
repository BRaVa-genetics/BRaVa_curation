mkdir -p $HOME/bin /well/lindgren/dpalmer/GRCh3{7,8} && cd /tmp

wget http://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
tar xjvf bcftools-1.20.tar.bz2
cd bcftools-1.20
./configure --disable-bz2 --prefix=/well/lindgren/dpalmer
make
make install

/bin/rm -f plugins/{score.{c,h},{munge,liftover,metal,blup}.c,pgs.{c,mk}}
wget -P plugins http://raw.githubusercontent.com/freeseek/score/master/{score.{c,h},{munge,liftover,metal,blup}.c,pgs.{c,mk}}
/bin/rm plugins/pgs.{c,mk}
make
/bin/cp bcftools plugins/{munge,liftover,score,metal,pgs,blup}.so $HOME/bin/
wget -P $HOME/bin http://raw.githubusercontent.com/freeseek/score/master/assoc_plot.R
chmod a+x $HOME/bin/assoc_plot.R

export PATH="$HOME/bin:$PATH"
export BCFTOOLS_PLUGINS="$HOME/bin"

module load BWA/0.7.17-GCCcore-12.3.0
module load SAMtools/1.18-GCC-12.3.0

wget -O- ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz | \
  gzip -d > /well/lindgren/dpalmer/GRCh37/human_g1k_v37.fasta
samtools faidx /well/lindgren/dpalmer/GRCh37/human_g1k_v37.fasta
bwa index /well/lindgren/dpalmer/GRCh37/human_g1k_v37.fasta
wget -P /well/lindgren/dpalmer/GRCh37 http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
wget -P /well/lindgren/dpalmer/GRCh37 http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain.gz
ref="/well/lindgren/dpalmer/GRCh37/human_g1k_v37.fasta"

wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | \
  gzip -d > /well/lindgren/dpalmer/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx /well/lindgren/dpalmer/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
bwa index /well/lindgren/dpalmer/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
wget -P /well/lindgren/dpalmer/GRCh38 http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
wget -P /well/lindgren/dpalmer/GRCh38 http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz
wget -P /well/lindgren/dpalmer/GRCh38 http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
ref="/well/lindgren/dpalmer/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

module load BEDTools
wget -O - "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz" |\
gunzip -c | grep 'transcript_type "protein_coding"' |\
awk '($3=="exon") {printf("%s\t%s\t%s\n",$1,int($4)-1,$5);}' |\
sort -T . -V -t $'\t' -k1,1 -k2,2n | bedtools merge > /well/lindgren/dpalmer/protein_coding_regions_hg38_no_padding_no_UTR_v39.bed
