# Modellvorhaben Genomsequenzierung


## Projekt�bersicht

### Grundlagen

Das Institut f�r Medizinische Genetik und Angewandte Genomik (IMGAG) f�hrt den Datenexport f�r das Modellvorhaben durch 

- f�r KDK-Daten des SE Netzwerks (KDK-Export f�r OE wird vom ZPM gemacht)
- f�r GRZ-Daten aller Netzwerke

Voraussetzungen/Zust�ndigkeiten:

- Modellvorhaben-Fall abgeschlossen/abgebrochen
- Documentation abgeschlossen in:
  - RedCap `Fallverwaltung Modellvorhaben Genomsequenzierung` f�r SE (ZSE) und OE (ZPM)
  - RedCap `Netzwerk Seltene Erkankungen` f�r SE (ZSE)
  - NGSD (NGS-Datenbank IMGAG)

### Schritt 1 � Datenaggregation + Pr�fung

Um den Datenexport durchzuf�hren wurde die Applikation MVHub entwickelt:

![MVHub](MVHub.png)

MVHub erm�glicht einen Gesamt�bersicht aller Daten zu den F�llen im Modellvorhaben.
Dazu werden f�r jeden Fall die Daten aus folgenden Quellen in einem Datenbank-Backend von MVHub gesammelt:

- Daten aus RedCap `Fallverwaltung Modellvorhaben Genomsequenzierung` via XML-Export
- Daten aus RedCap `Netzwerk Seltene Erkankungen` via XML-Export
- Proben-IDs der Genetik aus NGSD (NGS-Datenbank IMGAG)
- Einwilligung Forschung (REST API des meDIC f�r Broad Consent)

Die prim�re Datenquelle ist dabei das RedCap der Fallverwaltung.
Alle anderen Daten werden �ber die Patienten-ID von SAP mit den Daten der Fallverwaltung zusammengef�hrt.

Vor dem Datenexport muss gepr�ft werden ob die Daten bereit sind f�r den Upload:

- Fall abgeschlossen oder abgebrochen
- Grobe Pr�fung der Daten (alles vorhanden, keine Wiederspr�che)

### Schritt 2 - Upload

Nach positiver Pr�fung eines Datensatzes kann aus MVHub der Upload der Daten angesto�en werden.  
Dazu wird der Auftrag zum Export an KDK/GRZ in der MVHub-Datenbank hinterlegt.

Auf dem IMGAG Applikationsserver l�uft ein Background-Prozess der den eigentlichen Export durchf�hrt. Das Ergebnis jedes Exports (erfolgreich/abgebrochen und Ausgabe auf stdout/stderr) werden in der Datenbank gespeichert (siehe Screenshot). Falls ein Upload nicht erfolgreich ist, kann er nach Korrektur der Metadaten wiederholt werden.

F�r jeden Upload wird eine TAN �ber das Trustcenter des meDIC generiert.

### Schritt 3 - Meldebest�tigung

Nach erfolgreichem Upload, wird im RedCap � Fallverwaltung� die Information zum Upload hinterlegt (Datum, TAN, �):

![Pruefbericht](Pruefbericht.png)

Nach Pr�fung der hochgeladenen Daten schickt das KDK/GRZ einen Pr�fbericht ans BfArM.
Das BfArM schickt dann eine Meldebest�tigung mit der TAN per Email an ein Funktionspostfach des UKT.

ZSE/ZPM dokumentieren die Meldebest�tigung dann im RedCap und SAP:


![Meldebestaetigung](Meldebestaetigung.png)

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

Install mapper (outside of UKT network):

	> git clone https://github.com/KohlbacherLab/mii_broad_consent_mapper.git
	> ./gradlew build 
