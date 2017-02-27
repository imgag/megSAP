# *megSAP - RNA pipeline installation instructions*

**Note: These installation instructions are in alpha status. Do not try this at home!** 

### Step 1) Install DNA pipeline

The RNA pipeline uses some tools of the DNA pipeline.  
Thus it is assumed that you have [installed the DNA pipeline as documented](../Readme.md).

### Step 2) Install dependencies

For STAR, several non-standard Perl modules are required:

	> sudo apt install cpanminus libdb-dev
	> sudo -E cpanm Set::IntervalTree
	> sudo -E cpanm URI::Escape
	> sudo -E cpanm DB_File

### Step 3) Initial setup

First, we need to download the tools the pipeline relies on:

	> cd megSAP/data/
	> chmod 755 download_*.sh
	> ./download_tools_rna.sh

Finally, we need to download and some open-source databases for annotations and created genome indices for spliced mapping:

	> ./download_dbs_rna.sh

### Step 4) Test

Now you can test the installation by running the test:

	> cd megSAP/test/data_rna/
	> make










