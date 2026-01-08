//TODO

## Tools encapsulated in Apptainer container

The following tools are encapsulated in Apptainer container and used throughout the megSAP pipelines. All containers can be found in the `container_folder` specified in the settings.ini.
By default they are downloaded to `data/tools/apptainer_container/`. They can also be downloaded individually from `https://megsap.de/download/container/tool_version.sif`.

| container name                        | tool                      | version           | command                                                                                           |
|---------------------------------------|---------------------------|-------------------|---------------------------------------------------------------------------------------------------|
| abra2_v2.23.sif	                    | abra2                     | 2.23	            | java -Xmx16G -jar /opt/abra2.jar                                                                  |
| arriba_v2.5.1-20250904.sif            | arriba	                | 2.5.1	            | arriba                                                                                            |
| blastn_v2.9.0+.sif                    | blastn                    | 2.9.0+	        | blastn                                                                                            |
| bwa_v0.7.18.sif	                    | bwa                       | 0.7.18	        | bwa                                                                                               |
| bwa-mem2_v2.2.1.sif                   | bwa-mem2	                | 2.2.1             | bwa-mem2                                                                                          |
| circos_0.69.9.sif	                    | circos	                | 0.69.9	        | circos                                                                                            |
| clair3_v1.2.0-20251119.sif            | clair3	                | 1.2.0 	        | run_clair3.sh                                                                                     |
| clair3-trio_v0.7.sif                  | clair3-trio               | 0.7	            | /opt/bin/run_clair3_trio.sh                                                                       |
| ClinCNV_v1.18.3.sif	                | ClinCNV                   | 1.18.3	        | clinCNV.R                                                                                         |
| deepsomatic_1.9.0.sif	                | DeepSomatic               | 1.9.0 	        | run_deepsomatic                                                                                   |
| deepvariant_1.9.0.sif	                | DeepVariant               | 1.9.0 	        | run_deepvariant                                                                                   |
| expansionhunter_v5.0.0.sif            | expansionhunter	        | 5.0.0	            | ExpansionHunter                                                                                   |
| freebayes_v1.3.6.sif	                | freebayes                 | 1.3.6	            | freebayes                                                                                         |
| gatk_4.6.0.0.sif	                    | gatk                      | 4.6.0.0	        | gatk                                                                                              |
| glnexus_v1.4.1.sif                    | GLnexus                   | 1.4.1 	        | glnexus-cli                                                                                       |
| happy_v0.3.14.sif	                    | happy                     | 0.3.14	        | hap.py                                                                                            |
| hla-genotyper_2025-04.sif	            | hla-genotyper             | 2025-04	        | genotyper.py                                                                                      |
| htslib_1.16.sif                       | bgzip/tabix               | 1.16              | `bgzip` or `tabix`                                                                                |
| kraken2_v2.1.3.sif	                | kraken2	                | 2.1.3         	| kraken2                                                                                           |
| longphase_v1.7.3.sif	                | longphase                 | 1.7.3	            | longphase                                                                                         |
| manta_v1.6.0.sif	                    | manta                     | 1.6.0	            | python2 /opt/manta/bin/configManta.py                                                             |
| methylartist_v1.5.2-20251105.sif      | Methylartist              | 1.5.2             | methylartist                                                                                      |
| minimap2_v2.30-20251119.sif           | minimap2	                | 2.30	            | minimap2                                                                                          |
| modkit_v0.5.0-20251105.sif            | modkit	                | 0.5.0	            | modkit                                                                                            |
| msisensor-pro_v1.2.0.sif              | msisensor-pro             | 1.2.0	            | msisensor-pro                                                                                     |
| ngs-bits_2025-12-20251210.sif	        | ngs-bits	                | 2025-12	        | "tool_name" (e.g. BedAdd)                                                                         |
| orad_v2.6.1.sif	                    | orad	                    | 2.6.1	            | orad                                                                                              |
| paraphase_v3.3.1-20251125.sif         | paraphase                 | 3.3.1	            | paraphase                                                                                         |
| python_v3.10.9-20250729.sif           | python	                | 3.10.9	        | python3                                                                                           |
| REViewer_v0.2.7.sif                   | REViewer	                | 0.2.7	            | REViewer                                                                                          |
| samblaster_v0.1.26.sif	            | samblaster	            | 0.1.26	        | samblaster                                                                                        |
| samtools_1.20-20250812.sif            | samtools	                | 1.20	            | samtools                                                                                          |
| scarHRD_v1.sif	                    | scarHRD	                | 1	                | cli_scarHRD.R                                                                                     |
| SigProfilerExtractor_v1.1.24.sif      | SigProfilerExtractor      | 1.1.24	        | python3 -c 'from SigProfilerExtractor import sigpro as sig; sig.sigProfilerExtractor("parameter")'|
| sniffles_v2.7.1-20251112.sif          | sniffles	                | 2.7.1	            | sniffles                                                                                          |
| spliceai_v1.3.1.sif	                | spliceai	                | 1.3.1	            | spliceai                                                                                          |
| STAR_v2.7.11b.sif                     | STAR	                    | 2.7.11b	        | STAR                                                                                              |
| straglr_v1.5.5-20251119.sif	        | straglr	                | 1.5.5	            | straglr.py                                                                                        |
| straglrOn_v0.2.4-20250730.sif	        | straglrOn                 | 0.2.4	            | straglron.py                                                                                      |
| strelka2_v2.9.10.sif	                | strelka2                  | 2.9.10	        | python2 /opt/strelka2/bin/"script.py" (e.g.:runWorkflow.py)                                       |
| subread_v2.0.6.sif                    | subread                   | 2.0.6	            | featureCounts                                                                                     |
| umi-tools_v1.1.5.sif                  | umi-tools                 | 1.1.5	            | umi_tools                                                                                         |
| umiVar_2025-12-20251212.sif           | umiVar                    | 2025-12           | "script.py" (e.g. umiVar.py)                                                                      |
| varscan2_v2.4.6.sif	                | varscan2                  | 2.4.6	            | java -jar /opt/VarScan.jar                                                                        |
| vcflib_v1.0.3.sif	                    | vcflib	                | 1.0.3	            | "tool_name" (e.g. vcfallelicprimitives)                                                           |
| vep_release-112.0.sif                 | vep	                    | release-112.0	    | vep                                                                                               |
| whatshap_v2.8-20251119.sif            | whatshap                  | 2.8               | whatshap                                                                                          |

With Apptainer installed each container can be invoked with the following command:

    > singularity exec -B bind/paths tool_version.sif command parameters

## Building a new Apptainer container

### Writing the definition file

To build a new Apptainer container for a specific tool you have to first write a definition file. Below you can see a template definition file:

    Bootstrap: (e.g. docker)
    From: (e.g. ubuntu:20.04)

    %files
        <source>
        <source> <destination>

    %post
        LANG=C.UTF-8
        LC_ALL=C.UTF-8

        # Update package list and install build dependencies
        apt-get update --fix-missing 
        apt-get install -y \
            (e.g. wget) \
            (e.g. build-essential)
        
        # Download and build tool
        (install instructions for your tool)
        
        # Cleanup build dependencies
        cd /
        apt-get remove -y \
            (e.g. wget) \
            (e.g. build-essential)
        apt-get autoremove -y
        apt-get clean
        rm -rf /var/lib/apt/lists/*

    %environment
        export LANG=C.UTF-8
        export LC_ALL=C.UTF-8
        export PATH=/path/to/executable:$PATH

`HEADER`

In the header you define the bootstrap agent that will be used to create the base operating system you want to use.

`%files`

In this section you can copy files from the host system into the container. If you don't define a <destination> path it will be assumed to be the same as <source>.

`%post`

The %post section is where you actually download and build your software and all it's dependencies. 
To keep the container small and efficient, you should remove build-only dependencies after installing your tool.

`%environment`

In the %environment section you can define environment variables that will be set at runtime. These variables are only available in the container at runtime, 
but not at build time. Environment variables needed while building need to be defined in the %post section. 
It is recommended to add the files you want to execute to PATH here to simplify the invocation of the container.

More optional sections can be added to your definition file. For example a %runscript and/or %startscript section to define a default behavior when starting an instance of the container or when running it.
To learn about all possible sections and their use cases, visit https://apptainer.org/docs/user/latest/definition_files.html.

### Building the Apptainer container

To build your container from the created definition file execute:

    > apptainer build tool_version.sif definition_file.def

When creating new containers for megSAP they should be named `(tool)_(version).sif`, where `tool` is the name of the main tool installed inside the container and `version` its version number.
A new container must be added to the [apptainer-container] section of the settings files as `container_tool = version`.
The new container should be uploaded to https://megsap.de/download/container/ and the definition file should be added to the container definition file folder in the megSAP repository (data/tools/container_recipes).

### Executing tools inside your new container

To implement the invocation of a containerized tool in megSAP you can use the `execApptainer()` function:

    execApptainer(
                $container,                         tool
                $command,                           command to be executed inside the container (e.g. run_clair3.sh)
                $parameters,                        parameters for the given command (e.g. --bam_fn={$bam} --ref_fn=$genome --threads={$threads} ...)
                $in_files = array(),                host system directories/files used as input (needed for read file access from inside the container)
                $out_files = array(),               host system directories used as output (needed for read and write access from inside the container)
                $command_only=false,                if `true` only the singularity exec command is returned without executing it (needed when the tool execution is part of a pipeline)
                $log_output=true,                   Flag (true/false) to turn on/off logging of stdout, stderr and execution time of the containerised tool
                $abort_on_error=true,               Flag (true/false) whether to throw an error when execution fails
                $warn_on_error=true)                Flag (true/false) whether to throw a warning when execution fails
