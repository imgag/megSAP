# Installation of SpliceAI for GPU on Ubuntu 18.04


## Install python

Install pyhton 3.6.9. Newer versions lead to a compile error in pyfaidx dependency.

	> mkdir -p Python3.6.9
	> wget https://www.python.org/ftp/python/3.6.9/Python-3.6.9.tgz
	> tar -zxvf Python-3.6.9.tgz
	> cd Python-3.6.9 && ./configure --prefix=/mnt/users/ahsturm1/Sandbox/2021_12_01_update_splicing_scores/Python3.6.9
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



# Benchmarks

<table border=1>
  <tr><th>server</th><th># variants</th><th># threads</th><th>time</th></tr>
  <tr><td>SRV018 - no GPU</td><td>100</td><td>10</td><td>4:16</td></tr>
  <tr><td>SRV018 - no GPU</td><td>100</td><td>5</td><td>5:14</td></tr>
  <tr><td>SRV018 - no GPU</td><td>100</td><td>1</td><td>13:50</td></tr>
  <tr><td>SRV019 - 2x Tesla V100</td><td>100</td><td>10</td><td>1:56</td></tr>
  <tr><td>SRV019 - 2x Tesla V100</td><td>100</td><td>5</td><td>1:18</td></tr>
  <tr><td>SRV019 - 2x Tesla V100</td><td>100</td><td>1</td><td>1:16</td></tr>
  <tr><td>SRV019 - 2x Tesla V100</td><td>1000</td><td>10</td><td>4:35</td></tr>
  <tr><td>SRV019 - 2x Tesla V100</td><td>1000</td><td>5</td><td>4:28</td></tr>
  <tr><td>SRV019 - 2x Tesla V100</td><td>1000</td><td>1</td><td>4:09</td></tr>
  <tr><td>SRV019 - 2x Tesla V100 + genome on ramdisk</td><td>1000</td><td>10</td><td>4:26</td></tr>
  <tr><td>SRV019 - 2x Tesla V100 + genome on ramdisk</td><td>1000</td><td>5</td><td>3:52</td></tr>
  <tr><td>SRV019 - 2x Tesla V100 + genome on ramdisk</td><td>1000</td><td>1</td><td>4:09</td></tr>
  <tr><td>SRV019 - 2x Tesla V100</td><td>10000</td><td>10</td><td>29:05</td></tr>
  <tr><td>SRV019 - 2x Tesla V100</td><td>10000</td><td>5</td><td>29:01</td></tr>
  <tr><td>SRV019 - 2x Tesla V100</td><td>10000</td><td>1</td><td>24:26</td></tr>
  <tr><td>SRV019 - SpliceAI reforged via conda + genome on ramdisk</td><td>1000</td><td>10</td><td>10:41</td></tr>
  <tr><td>SRV019 - SpliceAI reforged via conda + genome on ramdisk</td><td>1000</td><td>5</td><td>10:18</td></tr>
  <tr><td>SRV019 - SpliceAI reforged via conda + genome on ramdisk</td><td>1000</td><td>1</td><td>9:48</td></tr>
  <tr><td>SRV019 - SpliceAI via conda + genome on ramdisk</td><td>1000</td><td>10</td><td></td></tr>
  <tr><td>SRV019 - SpliceAI via conda + genome on ramdisk</td><td>1000</td><td>5</td><td></td></tr>
  <tr><td>SRV019 - SpliceAI via conda + genome on ramdisk</td><td>1000</td><td>1</td><td></td></tr>
</table>
