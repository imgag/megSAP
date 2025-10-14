# megSAP src folder

The megSAP `src` folder contains the following sub-folders:

## Pipelines

This folder contains the main anaylsis pipelines.

## Common

This folder contains common functionaity used by many script, e.g. for file and database handling.
To include all functionality use 'all.php', it includes all other files.

## Install

This folder contains scripts that are used during installation only.

## Tools

This folder contains scripts that are used in pipelines or other scripts, e.g. wrappers for BWA, DeepVariant, etc.

## Auxilary

This folder contains auxilary scripts that are not part of the pipeline.  
They are mainly used by the developers.

## IMGAG

This folder contains scripts that are used by the developers in Tübingen only.

## MVH

This folder contains scripts used for the "Modellvorhaben Genomsequenzierung".

## Migration

This folder contains scrips that are used to migrate between megSAP versions.  
All scripts are prefixed with the date, so that it is easy to find out which scripts are relevant for a migration from version X to Y.

## Deprecated

This folder contains deprecated scrips that will be removed in the next release unless we get feedback that they are still needed.
