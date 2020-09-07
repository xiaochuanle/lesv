# Install `lesv`

In this section we introduce the installation and usage of `lesv` by calling structural variants on the human chromosome 21 with real (i.e. not simulated) Nanopore noisly long reads.

The experiment is conducted on a computer running Ubuntu 16.04:
```shell
$ uname -a
Linux zyserver 4.15.0-112-generic #113~16.04.1-Ubuntu SMP Fri Jul 10 04:37:08 UTC 2020 x86_64 x86_64 x86_64 GNU/Linux
```

Note that in this test we get very nice results. The precision, recall, and f1 metrics are 95%, 97% and 96% respectively!


[TOC]



## 1. Working Directory

We will work in the directory `/home/zyserver/chenying/my_sv`. For easy reference we assign it to the environment variable `SV_HOME`:
```shell
$ mkdir -p /home/zyserver/chenying/my_sv
$ cd /home/zyserver/chenying/my_sv
$ pwd
/home/zyserver/chenying/my_sv
$ export SV_HOME=$(pwd)
$ echo ${SV_HOME}
/home/zyserver/chenying/my_sv
```

## 2. Download Dataset

Download the dataset from [here](https://github.com/chenying2016/hg002_chr21_sv). It is recommanded that you download it from the web browser by cliking **Code -> Download ZIP**. Put the downloaded `zip` file into `${SV_HOME}`:
```shell
$ pwd
/home/zyserver/chenying/my_sv
$ ls
hg002_chr21_sv-master.zip
```

Unpack the `zip` data:
```shell
$ unzip hg002_chr21_sv-master.zip
$ ls
hg002_chr21_sv-master  hg002_chr21_sv-master.zip
$ cd hg002_chr21_sv-master/
$ cat hg002_chr21_sv.bz2.a* | tar xjv
$ rm hg002_chr21_sv.bz2.a*
$ cd hg002_chr21_sv/
$ ls
chr21_giab  chr21_raw_reads  chr21_ref
```

Data information:
1. `chr21_giab/chr21.bed` and `chr21_giab/chr21.vcf`: Structural Variations on human chromosome 21 extracted from `HG002_SVs_Tier1_v0.6.bed` and `HG002_SVs_Tier1_v0.6.vcf` respectively.
2. `chr21_raw_reads/chr00000020_queries.fasta`: Real (i.e. not simulated) Nanopore noisy long reads sampled from human chromosome 21. These reads are extracted from the dataset `HG002-normal-45x`.
3. `chr21_ref/21.fasta`: The human chromosome 21 sequence extracted from `hg19_sv_ref.fasta`.
4. `chr21_ref/chr21.trf.bed`: Tandem repeat annotations on `chr21_ref/21.fasta`. It is extracted from `human_hs37d5.trf.bed`.

Deails of datasets mentioned above can all be found in the **Benchmarks** Section.

We can have a quick look at the statistical information of the raw reads:
```shell
$ qx2viewdb hg002_chr21_sv-master/hg002_chr21_sv/chr21_raw_reads/chr00000020_queries.fasta
sequences:              76,617
residues:               1,229,047,537
max:                    152,906
min:                    3,000
avg:                    16,041
median:                 13,111

(N10, L10):             41,859    2,391     
(N20, L20):             33,230    5,711     
(N25, L25):             30,187    7,653     
(N30, L30):             27,459    9,791     
(N40, L40):             23,001    14,694    
(N50, L50):             19,362    20,523    
(N60, L60):             16,353    27,434    
(N70, L70):             13,812    35,620    
(N75, L75):             12,645    40,273    
(N80, L80):             11,582    45,350    
(N90, L90):             9,433     57,025   
```

## 3. Installing the `lesv` Pipeline

#### Step 1: Install Anaconda and bioconda
There are two different methods. Choose the one according to one of the following two cases.

*Case One:* **If you are not in Chaina mainland**, install `Anaconda` and add the `bioconda` channel according the the `Get tools` Section from [here](https://github.com/PacificBiosciences/sv-benchmark), and then proceed to Step 2.

*Case Two:* **If you are in Chaina mainland (如果你的位置位于中国内地)**, install `Anaconda` and add `bioconda` as follows:

Install `Anaconda`. If you have alreadly installed `Anaconda` like this:
```shell
$ which conda
/home/zyserver/anaconda3/bin/conda
```
Then you can skip the following commands. Otherwise please run:
```shell
$ cd ${SV_HOME}
$ pwd
/home/zyserver/chenying/my_sv
$ wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/archive/Anaconda3-5.3.1-Linux-x86_64.sh
$ chmod +x Anaconda3-5.3.1-Linux-x86_64.sh
$ ./Anaconda3-5.3.1-Linux-x86_64.sh
```

After installing `Anaconda`, add the following channels to `conda`:
```shell
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/pro
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
$ conda config --set show_channel_urls yes
```

#### Step 2: Creating `sv` Environment

Create the `sv` environment:
```shell
$ conda create -n sv python=3.6
```

You can remove the `sv` environment by (**please don't do it now!!**)
```shell
$ conda remove -n sv --all
```

Activate the `sv` environment:
```shell
$ conda activate sv
```

To deactivate or exit the `sv` environment, run (**please don't do it now!!**)
```shell
$ conda deactivate
```

#### Step 3: Install Auxiliary Tools

**Make sure you have activated the `sv` environemnt**.

Install `pbsv`:
```shell
$ conda install pbsv
$ pbsv --version
pbsv 2.3.0 (commit v2.3.0)
```

Install `pbmm2`:
```shell
$ conda install pbmm2
$ pbmm2 --version
pbmm2 1.3.0 (commit 1.3.0)
```

Install `samtools`. **`lesv` only works for `samtools` of version `0.1.19`.**
```shell
$ conda install samtools=0.1.19
```

#### Step 4: Install `truvari`

I fail to install [`truvari`](https://github.com/spiralgenetics/truvari) in the `sv` environment above. It tells me that the version of `Python` does not match. So I have to install it in a separate environment. To do so, I fist deactivate the `sv` environment:
```shell
$ conda deactivate
```
And then create a new one:
```shell
$ conda create -n truvari
```
And install [`truvari`](https://github.com/spiralgenetics/truvari) in that environment:
```
$ conda activate truvari
$ conda install truvari
$ conda deactivate
```


#### Step 5: Install `lesv`
```shell
$ pwd
/home/zyserver/chenying/my_sv
$ git clone https://github.com/xiaochuanle/lesv.git
$ cd lesv/src/
$ make -j
$ cd ../Linux-amd64/bin
$ export PATH=$PATH:$(pwd)
```
The last command above is used for adding `${SV_HOME}/lesv/Linux-amd64/bin` to the system `PATH`. If you encounter problems like
```shell
$ lesv.sh
lesv.sh: command not found
```
It is likely because you have not added `${SV_HOME}/lesv/Linux-amd64/bin` to the system `PATH`.

## 4. Run `lesv`

#### Step 1: Activate the `sv` environment if you have not:
```shell
$ conda activate sv
```

#### Step 2: Create a config file
```shell
$ lesv.sh config cfg
```
The command above creates a config file `cfg`:
```shell
$ echo cfg
PROJECT=
RAW_READS=
REFERENCE=
TRF_FILE=
THREADS=4

# split long read into short subreads
MAX_SUBSEQ_SIZE=50000
SUBSEQ_OVLP_SIZE=0
MIN_LAST_SUBSEQ_SIZE=20000

# reference mapping options
MAP_OPTIONS=

# sv read options
SVR_MIN_SEQ_SIZE=3000
SVR_MIN_SVE_PERC_IDENTITY=70.0
SVR_MAX_OVERHANG=300

# sv signature options
SVSIG_MIN_INDEL_SIZE=40
```

After setting the options in `cfg`, we have
```shell
$ cat cfg
PROJECT=sv_chr21
RAW_READS=/home/zyserver/chenying/my_sv/hg002_chr21_sv-master/hg002_chr21_sv/chr21_raw_reads/chr00000020_queries.fasta
REFERENCE=/home/zyserver/chenying/my_sv/hg002_chr21_sv-master/hg002_chr21_sv/chr21_ref/21.fasta
TRF_FILE=/home/zyserver/chenying/my_sv/hg002_chr21_sv-master/hg002_chr21_sv/chr21_ref/chr21.trf.bed
THREADS=16

# split long read into short subreads
MAX_SUBSEQ_SIZE=50000
SUBSEQ_OVLP_SIZE=0
MIN_LAST_SUBSEQ_SIZE=20000

# reference mapping options
MAP_OPTIONS="-max_target_seqs 5 -max_hsps 10"

# sv read options
SVR_MIN_SEQ_SIZE=3000
SVR_MIN_SVE_PERC_IDENTITY=70.0
SVR_MAX_OVERHANG=300

# sv signature options
SVSIG_MIN_INDEL_SIZE=40
```

#### Step 3: Run `lesv`
```shell
$ lesv.sh run cfg
```
The sv results can be found in `${SV_HOME}/sv_chr21/7-sv-call/pbsv/sv_chr21.pbsv.vcf.gz`.

#### Step 4: Assess the SV quality

We first deactivate and exit the `sv` environment and activate the `truvari` environment:
```shell
$ conda deactivate sv
$ conda activate truvari
```
The following shell script is used:
```shell
$ cat run_benchmarks.sh
#!/bin/bash

reference=/home/zyserver/chenying/my_sv/hg002_chr21_sv-master/hg002_chr21_sv/chr21_ref/21.fasta
giab_vcf=/home/zyserver/chenying/my_sv/hg002_chr21_sv-master/hg002_chr21_sv/chr21_giab/chr21.vcf
giab_vcf_bed=/home/zyserver/chenying/my_sv/hg002_chr21_sv-master/hg002_chr21_sv/chr21_giab/chr21.bed
pbsv_vcf=/home/zyserver/chenying/my_sv/sv_chr21/7-sv-call/pbsv/sv_chr21.pbsv.vcf.gz

truvari -f ${reference} -b ${giab_vcf} --includebed ${giab_vcf_bed} -o bench-pbsv --passonly --giabreport -r 1000 -p 0.00 -c ${pbsv_vcf}
```

We run the following command to assess the quality:
```shell
$ ./run_benchmarks.sh
```

If you encounter the following error
```shell
Traceback (most recent call last):
  File "/home/zyserver/anaconda3/envs/truvari/bin/truvari", line 990, in <module>
    run(sys.argv[1:])
  File "/home/zyserver/anaconda3/envs/truvari/bin/truvari", line 739, in run
    for entry in regions.iterate(vcf_base):
  File "/home/zyserver/anaconda3/envs/truvari/bin/truvari", line 555, in iterate
    if self.include(entry):
  File "/home/zyserver/anaconda3/envs/truvari/bin/truvari", line 567, in include
    return overlaps and len(self.tree[entry.CHROM].search(astart, aend)) == 1
AttributeError: 'IntervalTree' object has no attribute 'search'
```
Run the following command to fix it (remember we are still activated in the `truvari` environment). You may fail to install it due to TimeOut error. Please try a few more times until it is installed successfually!
```
$ pip install intervaltree==2.1.0
$ rm -r bench_pbsv
$ ./run_benchmarks.sh
```

If you encounter the following error
```shell
Traceback (most recent call last):
  File "/home/zyserver/anaconda3/envs/truvari/bin/truvari", line 990, in <module>
    run(sys.argv[1:])
  File "/home/zyserver/anaconda3/envs/truvari/bin/truvari", line 811, in run
    for comp_entry in vcf_comp.fetch(base_entry.CHROM, max(0, fetch_start - 1), fetch_end + 1):
  File "/home/zyserver/anaconda3/envs/truvari/lib/python2.7/site-packages/vcf/parser.py", line 620, in fetch
    raise Exception('pysam not available, try "pip install pysam"?')
Exception: pysam not available, try "pip install pysam"?
```
Run the following command to fix it:
```shell
$ conda install pysam
$ rm -r bench_pbsv
$ ./run_benchmarks.sh
```

And we finally get the benchmark results:
```shell
2020-09-06 15:21:56,289 [INFO] Stats: {
    "FP": 9, 
    "f1": 0.9608938547486033, 
    "base gt filtered": 0, 
    "precision": 0.9502762430939227, 
    "TP-call": 172, 
    "call cnt": 181, 
    "base size filtered": 107, 
    "FN": 5, 
    "TP-base": 172, 
    "base cnt": 177, 
    "recall": 0.9717514124293786, 
    "call gt filtered": 0, 
    "call size filtered": 36
}
```
And of course we are still encountering new error (Python is a marvelous language!):
```shell
2020-09-06 15:21:56,289 [INFO] Creating GIAB report
Traceback (most recent call last):
  File "/home/zyserver/anaconda3/envs/truvari/bin/truvari", line 990, in <module>
    run(sys.argv[1:])
  File "/home/zyserver/anaconda3/envs/truvari/bin/truvari", line 985, in run
    make_giabreport(args, stats_box)
  File "/home/zyserver/anaconda3/envs/truvari/bin/truvari", line 454, in make_giabreport
    count_by(size_keys, tp_base, sum_out)
  File "/home/zyserver/anaconda3/envs/truvari/bin/truvari", line 372, in count_by
    cnt[x[main_key]] += 1
TypeError: unhashable type: 'list'
```
I don't know how to fix it anymore. However, this error will not happen when we assess the quality on the whose human genomes (see the **Benchmarks** Section).
