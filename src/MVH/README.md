# Modellvorhaben Genomsequenzierung

## Installation of GRZ QC workflow (not used right now because of the Nextflow problems

### Installation nexflow

	> cd /mnt/storage2/MVH/tools/nextflow
	> curl -s https://get.nextflow.io | bash

### Installation of nf-core tools

	> pipx install nf-core

### Installation of GRZ QC pipeline

	> export PATH=$PATH:/mnt/storage2/MVH/tools/nextflow/
	> nf-core pipelines download BfArM-MVH/GRZ_QC_Workflow --container-system singularity
	> nextflow plugin install nf-schema@2.1.1


## Installation of tools (for running this script without GRZ QC workflow)

### Installation of mosdepth

	> wget https://github.com/brentp/mosdepth/releases/download/v0.3.11/mosdepth --no-check-certificate

### Installation of fastp

	> wget http://opengene.org/fastp/fastp
	> chmod a+x ./fastp

### Installation of fastplong

	> wget http://opengene.org/fastplong/fastplong
	> chmod a+x ./fastplong

## Installation of GRZ-CLI

see <https://github.com/BfArM-MVH/grz-tools/blob/main/packages/grz-cli/README.md> for details

- Install miniforge at /mnt/storage2/megSAP/tools/miniforge3/
	> curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
	> bash Miniforge3-$(uname)-$(uname -m).sh
- Install GRZ-CLI		
	> /mnt/storage2/MVH/tools/miniforge3/bin/conda create -n grz-tools -c conda-forge -c bioconda "grz-cli"
	> /mnt/storage2/MVH/tools/miniforge3/bin/conda activate grz-tools
- Updates with:
	> /mnt/storage2/MVH/tools/miniforge3/bin/conda update -n base -c conda-forge conda -c bioconda
	> /mnt/storage2/MVH/tools/miniforge3/bin/conda update -n grz-tools -c conda-forge -c bioconda grz-cli
	> cd /mnt/storage2/MVH/tools/GRZ_QC_Workflow && git pull
- List all package versions:
	> /mnt/storage2/MVH/tools/miniforge3/bin/conda list -n grz-tools 

## Installation of python3

	> python3 -m venv /mnt/storage2/MVH/tools/python
	> /mnt/storage2/MVH/tools/python3/bin/pip install grz-pydantic-models pandas argparse importlib


## Installation of consent mapper

Install gradle:

	> mkdir /mnt/storage2/MVH/tools/gradle
	> wget https://services.gradle.org/distributions/gradle-8.10.2-bin.zip -P /tmp
	> unzip -d /mnt/storage2/MVH/tools/gradle /tmp/gradle-8.10.2-bin.zip
	> echo 'export PATH=/mnt/storage2/MVH/tools/gradle/latest/bin:$PATH' >> ~/.bashrc

Create `~/.gradle/gradle.properties` and add the following lines to set up the proxy for gradle:

	# === HTTP Proxy Settings ===
	systemProp.http.proxyHost=httpproxy.zit.med.uni-tuebingen.de
	systemProp.http.proxyPort=88
	systemProp.http.proxyUser=AH1inges
	systemProp.http.proxyPassword=iD5PvnGy+Vm@

	# === HTTPS Proxy Settings ===
	systemProp.https.proxyHost=httpproxy.zit.med.uni-tuebingen.de
	systemProp.https.proxyPort=88
	systemProp.https.proxyUser=AH1inges
	systemProp.https.proxyPassword=iD5PvnGy+Vm@

	# === Exclusions (don't proxy local/internal addresses) ===
	systemProp.http.nonProxyHosts=localhost|127.0.0.1|*.example.local

	# === Optional: Increase timeout (useful behind slow proxies) ===
	systemProp.http.connectionTimeout=60000
	systemProp.http.readTimeout=60000

Install mapper:

	> git clone https://github.com/KohlbacherLab/mii_broad_consent_mapper.git
	> gradle wrapper
	> ./gradlew build 
	
## TODOs

- add tests when first final version is done
  - KDK-SE: WGS, lrGS, no_seq
  - GRZ: SE WGS, SE lrGS, SE WGS trio, T/N
	
