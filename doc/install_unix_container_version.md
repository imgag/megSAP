# Building megSAP Container Version

## Dependencies

These instructions are for **Ubuntu 20.04**, but as long as you can install **Apptainer**, the container should work on other systems as well.

## Install Apptainer (Ubuntu 20.04)

```sh
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
```

## Downloading

Download the megSAP Apptainer container:

```sh
wget --no-check-certificate -O megSAP_master.sif https://megsap.de/download/container/megSAP_master.sif
```

## Initial Setup

Although megSAP is encapsulated within an Apptainer container, some required tools, genomes, and databases must be downloaded separately.

### Create a Data Folder

First, create a folder on your **host system** to store the downloaded data:

```sh
mkdir -p <path-to-host-data-folder>
```

### Download Tools and Containers

Next, download a few tools and Apptainer containers for the rest of the tools:

```sh
apptainer exec --pwd /megSAP/data -B <path-to-host-data-folder>:/megSAP/data/data_folder/ megSAP_[version].sif ./download_tools.sh
apptainer exec --pwd /megSAP/data -B <path-to-host-data-folder>:/megSAP/data/data_folder/ megSAP_[version].sif ./download_container.sh
```

### Download and Index the Reference Genome

```sh
apptainer exec --pwd /megSAP/data -B <path-to-host-data-folder>:/megSAP/data/data_folder/ megSAP_[version].sif ./download_GRCh38.sh
```

### Download and Convert Databases

```sh
apptainer exec --pwd /megSAP/data -B <path-to-host-data-folder>:/megSAP/data/data_folder/ megSAP_[version].sif ./download_dbs.sh
apptainer exec -B <path-to-host-data-folder>:/megSAP/data/data_folder/ megSAP_[version].sif php /megSAP/src/Install/db_download.php -data_folder /megSAP/data/data_folder/
```

**Note:** OMIM, HGMD, and COSMIC databases are not downloaded automatically due to licensing restrictions. If you have the required licenses, follow the instructions in the [download_dbs.sh](../data/download_dbs.sh) script to download and convert them manually. You can use the ngs-bits container downloaded earlier for the ngs-bits command and the megSAP container to execute the php scripts.

## Settings

Changing the settings is **optional**, as most entries have default values. If needed, copy the default settings file from within the container to your host system:

```sh
apptainer exec megSAP_[version].sif cp /megSAP/settings.ini.default ./settings.ini
```

Now, edit `settings.ini` as needed. **Do not change any paths**—they reference directories **inside** the megSAP container. Once modified, bind your settings file when invoking the container:

```sh
-B <path-to-new-settings.ini>:/megSAP/settings.ini
```

For a detailed description of settings, refer to [settings.md](settings.md).

## NGSD Initialization

If you plan to use **NGSD** and it is not initialized yet, follow these steps:

1. Install a MariaDB server.

2. Create a database and user in the SQL database.

3. Add NGSD connection details to `settings.ini`.

4. Create tables using:

   ```sh
   apptainer exec -B <path-to-new-settings.ini>:/megSAP/settings.ini,<path-to-host-data-folder>:/megSAP/data/data_folder/ megSAP_[version].sif php /megSAP/src/Install/db_init.php
   ```

   During this step, a `ngsbits_settings.ini` is created from the megSAP settings.ini and saved under `<path-to-host-data-folder>/tools/ngsbits_settings.ini`. This `ngsbits_settings.ini` file is required to run the ngs-bits container independently of the megSAP container.

5. Import base genomic data using `ngs-bits` tools:

   ```sh
   NGSDImportQC --help  
   NGSDImportHGNC --help  
   NGSDImportEnsembl --help  
   NGSDImportHPO --help  
   NGSDImportGeneInfo --help  
   NGSDImportOMIM --help  
   NGSDImportORPHA --help  
   ```

**Note:** To call `ngs-bits` tools, use: `apptainer exec ngs-bits_[version].sif [tool] [parameters] --settings <path-to-host-data-folder>/tools/ngsbits_settings.ini`. The `ngs-bits` container is downloaded to `<path-to-host-data-folder>/tools/apptainer_container`.

**Note:** To annotate variants with NGSD in-house counts, classifications, etc., NGSD data has to be exported regularly. To do so, adapt the file `data\dbs\NGSD\Makefile` and execute `make export` once a week using a cronjob. To do so in the megSAP container version you have to copy the Makefile to your host system, change it and mount it when executing the `make export`.

## Execution

Now, all required tools and data are ready. You can list available pipelines:

```sh
apptainer exec megSAP_[version].sif ls /megSAP/src/Pipelines
```

Refer to the [documentation](../README.md) for more details.

### Running a Pipeline

Execute a pipeline with the following command:

```sh
apptainer exec -B \
    <path-to-host-data-folder>:/megSAP/data/data_folder/, \
    <path-to-new-settings.ini>:/megSAP/settings.ini, \
    </tmp/folder/>:/tmp/local_ngs_data/, \
    </path/to/input_data/>:/working_dir/ \
    megSAP_[version].sif \
    php /megSAP/src/Pipelines/analyze.php \
        -folder /working_dir \
        -name [processed-sample-name] \
        -steps ma,vc,cn,sv,re \
        -threads [number-of-threads]
```

### Explanation of Bind Mounts:

- ``<path-to-host-data-folder> → /megSAP/data/data_folder/``: Mounts the data folder containing the downloaded data.
- ``<path-to-new-settings.ini> → /megSAP/settings.ini``: Mounts your modified settings file.
- ``</tmp/folder/> → /tmp/``: Mounts your tmp folder for intermediate and locally stored files.
- ``</path/to/input_data/> → /working_dir/``: Mounts the folder with input data.

**Note:** If your host system's temporary directory is not ``/tmp`` change tmp folder mounting to ``</tmp/folder/> → </tmp/folder/>`` and change ``local_data = /tmp/local_ngs_data/`` to ``local_data = </tmp/folder/>/local_ngs_data`` in your `settings.ini`

**Note:** This example shows how to run `analyze.php`. To check other pipelines, use: `apptainer exec megSAP_[version].sif php /megSAP/src/Pipelines/[pipeline_name].php --help`
