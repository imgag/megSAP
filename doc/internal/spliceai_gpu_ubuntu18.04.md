# Installation of SpliceAI for GPU on Ubuntu 18.04


## Install python

Install pyhton 3.6.9. Newer versions lead to a compile error in pyfaidx dependency.

	> mkdir -p Python3.6.9
	> wget https://www.python.org/ftp/python/3.6.9/Python-3.6.9.tgz
	> tar -zxvf Python-3.6.9.tgz
	> cd Python-3.6.9 && ./configure --prefix=/mnt/storage3/users/ahsturm1/Sandbox/2021_12_01_update_splicing_scores/Python3.6.9
	> cd Python-3.6.9 && make
	> cd Python-3.6.9 && make install
	> rm -R Python-3.6.9
	> rm Python-3.6.9.tgz

## Install SpliceAI

To use GPUs, you have to install `tensorflow-gpu` instead of the normal `tensorflow` package.  
Note: `keras` has to be installed in a lower version (2.3.1) instead of (2.7.0) to avoid an error finding the 'get_config' function.

	> mkdir -p spliceai
	> cd spliceai && ../Python3.6.9/bin/python3 -m venv splice_env
	> source spliceai/splice_env/bin/activate
	> pip install keras==2.3.1
	> pip install spliceai 
	> pip install tensorflow-gpu
	> deactivate

## Install CUDA 10.0

Each tensorflow library needs a specific cuda toolkit installation.  
Here we need version 10.0:

	> sudo apt-get install cuda-toolkit-10.0

## Install cudnn library

Finally `libcudnn7` is needed, which is not available in the package system.  
Thus we need to download and install it manually:

	> wget https://developer.nvidia.com/compute/machine-learning/cudnn/secure/7.6.5.32/Production/10.0_20191031/Ubuntu18_04-x64/libcudnn7_7.6.5.32-1%2Bcuda10.0_amd64.deb
	> sudo dpkg -i /mnt/share/software/libcudnn7_7.6.5.32-1+cuda10.0_amd64.deb

## SpliceAI reforged

Code and installation instructions can be found here:  
https://github.com/skoblov-lab/spliceai-reforged
Installation was done using conda.


# Benchmarks

Benchmark of 1000 variants:

<table border=1>
  <tr><th>server</th><th># threads</th><th>time</th></tr>
  <tr><td>SRV011 - no GPU</td><td>1</td><td>24:19</td></tr>
  <tr><td>SRV011 - no GPU</td><td>5</td><td>25:11</td></tr>
  <tr><td>SRV011 - no GPU</td><td>10</td><td>25:30</td></tr>
  <tr><td>SRV019 - 2x Tesla V100</td><td>1</td><td>4:09</td></tr>
  <tr><td>SRV019 - 2x Tesla V100</td><td>5</td><td>4:28</td></tr>
  <tr><td>SRV019 - 2x Tesla V100</td><td>10</td><td>4:35</td></tr>
  <tr><td>SRV019 - 2x Tesla V100 - SpliceAI reforged</td><td>1</td><td>9:48</td></tr>
  <tr><td>SRV019 - 2x Tesla V100 - SpliceAI reforged</td><td>5</td><td>10:18</td></tr>
  <tr><td>SRV019 - 2x Tesla V100 - SpliceAI reforged</td><td>10</td><td>10:41</td></tr>
</table>

# Conclusion

Here the conclusions from the benchmarks:

- The GPU version is about 5 times faster than the CPU version.
- SpliceAI reforged was not faster than the normal GPU version of SpliceAI.
- Neither the CPU version nor the GPU version of SpliceAI scale with the number of threads > start several instances of SpliceAI if you want to parallelize.
- Putting the genome on a ramdisk was also tested and gave no imporovements.


