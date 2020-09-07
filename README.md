# lesv

`lesv` is a high accuracy structural variantions caller for Nanopore noisy long reads.
For a 45x human raw read dataset, `lesv` achieves 93%, 91% and 92% of precesion, recall and f1 metrics, respectively.
For a 40x human ultra-long raw read >100k dataset, `lesv` achieves 94%, 97%, and 95% of precision, recall and f1 metrics, respectively.
Details can be found in the [Benchmarks](#Benchmarks) below.

`lesv` currently only supports Deletions and Insertions Structural Variants. More types of SVs will be implemented subsequently.


## Installation

Installation and usage of `lesv` can be found [here](https://github.com/xiaochuanle/lesv/blob/master/install_lesv.md).

## Benchmarks

### Experimental Setup

Our experiments are all conducted on a computer with 48 2.50GHz CPUs and 512GB memory. The operating system is Ubuntu 16.04.

### Datasets

We carry out benchmarks on two datasets.

* The first dataset is generated using our in-house sequencing and will be released soon. We call it `HG002-normal-45x`.

* The second dataset consists of 40x 100kb+ ultra-long raw reads. It is extracted from the publically released data [HG002_ucsc_ONT_lt100kb.fastq.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/HG002/hpp_HG002_NA24385_son_v1/nanopore/downsampled/greater_than_100kb/HG002_ucsc_ONT_lt100kb.fastq.gz) We call it`HG002-100k-45x`.

#### Statistical Information of `HG002-normal-45x`

```shell
sequences:              10,546,295
residues:               134,084,883,764
max:                    69,999
min:                    1
avg:                    12,713
median:                 10,663

(N10, L10):             42,081    266,833   
(N20, L20):             33,018    629,889   
(N25, L25):             29,853    843,677   
(N30, L30):             27,150    1,079,309 
(N40, L40):             22,679    1,620,717 
(N50, L50):             18,994    2,267,411 
(N60, L60):             15,861    3,040,812 
(N70, L70):             13,216    3,967,872 
(N75, L75):             12,055    4,499,182 
(N80, L80):             10,982    5,081,991 
(N90, L90):             8,276     6,454,807 
```

#### Statistical Information of `HG002-100k-40x`

```shell
sequences:              744,007
residues:               120,002,274,309
max:                    2,432,464
min:                    120,142
avg:                    161,291
median:                 140,392

(N10, L10):             299,074    29,408    
(N20, L20):             214,579    77,867    
(N25, L25):             193,630    107,388   
(N30, L30):             179,214    139,657   
(N40, L40):             160,609    210,697   
(N50, L50):             148,861    288,466   
(N60, L60):             140,430    371,568   
(N70, L70):             133,932    459,135   
(N75, L75):             131,128    504,417   
(N80, L80):             128,576    550,634   
(N90, L90):             124,057    645,684  
```

#### Sequencing Error Distribution


| Error Rate    | `HG002-normal-45x`  | `HG002-100k-40x` |
| ------------: |---------------:     | ---------------: |
| 1%            |   0.23%             |    0.59%         |
| 2%            |   0.21%             |    0.39%         |
| 3%            |   0.26%             |    0.73%         |
| 4%            |   0.52%             |    1.74%         |
| 5%            |   1.10%             |    3.43%         |
| 6%            |   2.31%             |    5.03          |
| 7%            |   3.51%             |    6.34%         |
| 8%            |   4.85%             |    6.17%         |
| 9%            |   7.00%             |    5.28%         |
| 10%           |   8.52%             |    4.39%         |
| 10-15%        |  32.28%             |    13.35%        |
| 15-20%        |  14.33%             |    8.02%         |
| 20-25%%       |   8.66%             |    9.39%         |
| 25-30%%       |   6.85%             |    14.37%        |
| >30%          |   9.36%             |    20.76%        |

#### Download Reference Genome

Download hg19 reference with decoys and map non-ACGT characters to N following the same way as [PacBio CCS 15kb dataset](https://github.com/PacificBiosciences/sv-benchmark).

```
curl -s ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz > ref/human_hs37d5.fasta.gz
gunzip ref/human_hs37d5.fasta.gz
sed -i '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' ref/human_hs37d5.fasta
```
We further extract the 24 chrosomes from `human_hs37d5.fasta` and store them in `hg19_sv_ref.fasta`.

#### Download the Repeat Annotations

```shell
curl -s https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_hs37d5.trf.bed > ref/human_hs37d5.trf.bed
```

#### Download the giab SV results

```shell
FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/
curl -s ${FTPDIR}/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed > giab/HG002_SVs_Tier1_v0.6.bed
curl -s ${FTPDIR}/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz > giab/HG002_SVs_Tier1_v0.6.vcf.gz
```

### Running Configurations

#### `HG002-normal-45x` configurations

```shell
PROJECT=hg002_short
RAW_READS=/data1/cy/ontsv/HG002_NA24385_Son/read_list.txt
REFERENCE=/data1/cy/hg19/hg19_sv_ref.fasta
TRF_FILE=/data1/cy/run_ccs/ref/human_hs37d5.trf.bed
THREADS=48

# split long read into short subreads
MAX_SUBSEQ_SIZE=50000
SUBSEQ_OVLP_SIZE=0
MIN_LAST_SUBSEQ_SIZE=20000

# reference mapping options
MAP_OPTIONS="-max_query_vol_res 8g -query_batch_size 2g -max_target_seqs 5 -max_hsps 10"

# sv read options
SVR_MIN_SEQ_SIZE=3000
SVR_MIN_SVE_PERC_IDENTITY=70.0
SVR_MAX_OVERHANG=300

# sv signature options
SVSIG_MIN_INDEL_SIZE=40
```

#### `HG002-100k-40x` configuations

```shell
PROJECT=hg002_100k
RAW_READS=/data1/cy/data/human-long-fltd.fasta.gz
REFERENCE=/data1/cy/hg19/hg19_sv_ref.fasta
TRF_FILE=/data1/cy/run_ccs/ref/human_hs37d5.trf.bed
THREADS=48

# split long read into short subreads
MAX_SUBSEQ_SIZE=50000
SUBSEQ_OVLP_SIZE=0
MIN_LAST_SUBSEQ_SIZE=20000

# reference mapping options
MAP_OPTIONS="-kmer_size 19 -kmer_window 20 -max_query_vol_res 8g -query_batch_size 2g -max_target_seqs 5 -max_hsps 10"

# sv read options
SVR_MIN_SEQ_SIZE=3000
SVR_MIN_SVE_PERC_IDENTITY=80.0
SVR_MAX_OVERHANG=300

# sv signature options
SVSIG_MIN_INDEL_SIZE=40
```

### Experimental Results

#### `HG002-normal-45x` results

The following shell script is used for evaluating the SV results on `HG002-normal-45x`:
```shell
$ cat evaluate_hg002_short.sh
pbsv_vcf="/data1/cy/run_short_ont/hg002_short/7-sv-call/pbsv/hg002_short.pbsv.vcf.gz"
/data1/cy/run_ccs/truvari/truvari.py -f /data1/cy/hg19/hg19_sv_ref.fasta -b giab/HG002_SVs_Tier1_v0.6.vcf --includebed giab/HG002_SVs_Tier1_v0.6.bed -o bench-pbsv --passonly --giabreport -r 1000 -p 0.00 -c ${pbsv_vcf}
```

```shell
    "TP-base": 8815,
    "TP-call": 8815,
    "FP": 669,
    "FN": 826,
    "base cnt": 9641,
    "call cnt": 9484,
    "base size filtered": 6309,
    "call size filtered": 2387,
    "base gt filtered": 0,
    "call gt filtered": 0,
    "precision": 0.9294601433994095,
    "recall": 0.9143242402240431,
    "f1": 0.9218300653594771
```

#### `HG002-100k-40x` results

The following shell script is used for evaluating the SV results on `HG002-100k-40x`:
```shell
$ cat evaluate_hg002_long.sh
pbsv_vcf="/data1/cy/run_long_ont/hg002_100k/7-sv-call/pbsv/hg002_100k.pbsv.vcf.gz"
/data1/cy/run_ccs/truvari/truvari.py -f /data1/cy/hg19/hg19_sv_ref.fasta -b giab/HG002_SVs_Tier1_v0.6.vcf --includebed giab/HG002_SVs_Tier1_v0.6.bed -o bench-pbsv --passonly --giabreport -r 1000 -p 0.00 -c ${pbsv_vcf}
```

```shell
    "TP-base": 9340,
    "TP-call": 9340,
    "FP": 591,
    "FN": 301,
    "base cnt": 9641,
    "call cnt": 9931,
    "base size filtered": 6309,
    "call size filtered": 1962,
    "base gt filtered": 0,
    "call gt filtered": 0,
    "precision": 0.9404893766992246,
    "recall": 0.9687791722850326,
    "f1": 0.9544246883302677
```

#### `PacBio CCS 15kb dataset` results

We reproducing the benchmarks using the same way as in [PacBio CCS 15kb dataset](https://github.com/PacificBiosciences/sv-benchmark). We obtain a similar but slightly different metrics.

```shell
    "TP-base": 9441,
    "TP-call": 9441,
    "FP": 685,
    "FN": 200,
    "base cnt": 9641,
    "call cnt": 10126,
    "base size filtered": 6309,
    "call size filtered": 4274,
    "base gt filtered": 0,
    "call gt filtered": 0,
    "precision": 0.932352360260715,
    "recall": 0.9792552639767659,
    "f1": 0.9552284109880103
```











