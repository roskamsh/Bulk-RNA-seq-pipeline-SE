# Co-culture Bulk RNA-sequencing pipeline SE

Pipeline to run basic RNA-seq analysis on single-end data, which was performed in a xenograft mouse model, and incorporates filtering out of contaminating mouse reads.

This package implements the Snakemake management workflow system and is currently implemented to work with 
the cluster management and job scheduling system SLURM. This snakemake workflow utilizes conda installations to download and use packages for further analysis, so please ensure that you have installed miniconda prior to use.

Questions/issues
======================

Please add an issue to the Co-culture-Bulk-RNA-seq-pipeline repository. We would appreciate if your issue included sample code/files (as appropriate) so that we can reproduce your bug/issue. 


Contributing
======================

We welcome contributors! For your pull requests, please include the following:

* Sample code/file that reproducibly causes the bug/issue
* Documented code providing fix
* Unit tests evaluating added/modified methods. 

Use
======================

Locate raw files:
* After sequencing, your raw fastq files are placed in `/path/to/sequencing/files`.

```
$ cd /path/to/raw/data
$ ls -alh
```

Check md5sum.

```
$ md5sum –c md5sum.txt > md5sum_out.txt
```

Move your files into the archive to be stored.

```
$ mv /path/to/raw/data /path/to/archive
```

Check md5sum again to ensure your sequencing files are not corrupted.

```
$ md5sum –c md5sum.txt > md5sum_out.txt
```

Clone this Pipeline into your working directory.

```
$ git clone https://github.com/ohsu-cedar-comp-hub/Bulk-RNA-seq-pipeline-SE.git
```

Create a `samples/raw` directory, and a `logs` directory in your `wdir()`.

```
$ mkdir logs
$ mkdir samples
$ cd samples
$ mkdir raw
```

Symbollically link the fastq files of your samples to the `wdir/samples/raw` directory using a bash script loop in your terminal.

```
ls -1 /path/to/data/LIB*gz |     while read gz; do
  R=$( basename $gz | cut -d '_' -f 3 | awk '{print $1".fastq.gz"}' )
  echo $R
  ln -s ${gz} ./${R}
done
```

Upload your metadata file to the `data` directory, with the correct formatting:
* Columns should read:
```StudyID  Column2   Column3   ...```
* Each row should be a sample, with subsequent desired information provided (RNA extraction date, etc.)
* Edit omic_config.yaml to include only columns included in this metadata file:
  * This includes `meta_columns_to_plot` and `pca labels`
* All values in this file should be tab-separated

Edit the `omic_config.yaml` in your `wdir()`:
* Change the `project_id` to a unique project identifier
* Add appropriate contrasts based on your samples under the `[diffexp][contrasts]` section
* Add the path to your metadata file for the `omic_meta_data` and `samples` parameters
* Change `base_dir` to your current working directory
* Ensure you have the correct `assembly` specified
    * Current options for this are: hg19, hg38.89 (ensembl v89) and hg38.90 (ensembl v90)

Do a dry-run of snakemake to ensure proper execution before submitting it to the cluster (in your wdir).

```
$ snakemake -np --verbose
```

Once your files are symbolically linked, you can submit the job to exacloud via your terminal window.

```
$ sbatch submit_snakemake.sh
```

To see how the job is running, look at your queue.

```
$ squeue -u your_username
```


