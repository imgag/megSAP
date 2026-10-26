# Modellvorhaben Genomsequenzierung


## Projektübersicht

### Grundlagen

Das Institut für Medizinische Genetik und Angewandte Genomik (IMGAG) führt den Datenexport für das Modellvorhaben durch 

- für KDK-Daten des SE Netzwerks (KDK-Export für OE wird vom ZPM gemacht)
- für GRZ-Daten aller Netzwerke

Voraussetzungen/Zuständigkeiten:

- Modellvorhaben-Fall abgeschlossen/abgebrochen
- Documentation abgeschlossen in:
  - RedCap `Fallverwaltung Modellvorhaben Genomsequenzierung` für SE (ZSE) und OE (ZPM)
  - RedCap `Netzwerk Seltene Erkankungen` für SE (ZSE)
  - NGSD (NGS-Datenbank IMGAG)

### Schritt 1 – Datenaggregation + Prüfung

Um den Datenexport durchzuführen wurde die Applikation MVHub entwickelt:

![MVHub](MVHub.png)

MVHub ermöglicht einen Gesamtübersicht aller Daten zu den Fällen im Modellvorhaben.
Dazu werden für jeden Fall die Daten aus folgenden Quellen in einem Datenbank-Backend von MVHub gesammelt:

- Daten aus RedCap `Fallverwaltung Modellvorhaben Genomsequenzierung` via XML-Export
- Daten aus RedCap `Netzwerk Seltene Erkankungen` via XML-Export
- Proben-IDs der Genetik aus NGSD (NGS-Datenbank IMGAG)
- Einwilligung Forschung (REST API des meDIC für Broad Consent)

Die primäre Datenquelle ist dabei das RedCap der Fallverwaltung.
Alle anderen Daten werden über die Patienten-ID von SAP mit den Daten der Fallverwaltung zusammengeführt.

Vor dem Datenexport muss geprüft werden ob die Daten bereit sind für den Upload:

- Fall abgeschlossen oder abgebrochen
- Grobe Prüfung der Daten (alles vorhanden, keine Wiedersprüche)

### Schritt 2 - Upload

Nach positiver Prüfung eines Datensatzes kann aus MVHub der Upload der Daten angestoßen werden.  
Dazu wird der Auftrag zum Export an KDK/GRZ in der MVHub-Datenbank hinterlegt.

Auf dem IMGAG Applikationsserver läuft ein Background-Prozess der den eigentlichen Export durchführt. Das Ergebnis jedes Exports (erfolgreich/abgebrochen und Ausgabe auf stdout/stderr) werden in der Datenbank gespeichert (siehe Screenshot). Falls ein Upload nicht erfolgreich ist, kann er nach Korrektur der Metadaten wiederholt werden.

Für jeden Upload wird eine TAN über das Trustcenter des meDIC generiert.

### Schritt 3 - Meldebestätigung

Nach erfolgreichem Upload, wird im RedCap „ Fallverwaltung” die Information zum Upload hinterlegt (Datum, TAN, …):

![Pruefbericht](Pruefbericht.png)

Nach Prüfung der hochgeladenen Daten schickt das KDK/GRZ einen Prüfbericht ans BfArM.
Das BfArM schickt dann eine Meldebestätigung mit der TAN per Email an ein Funktionspostfach des UKT.

ZSE/ZPM dokumentieren die Meldebestätigung dann im RedCap und SAP:


![Meldebestaetigung](Meldebestaetigung.png)

## Installation of GRZ QC workflow (for running this script without GRZ QC workflow, which is not used because of the NextFlow proxy issues)

See container definition at <https://github.com/imgag/megSAP/blob/master/data/tools/container_recipes/grz_qc_workflow_master.def>

## Installation of GRZ-CLI

see <https://github.com/BfArM-MVH/grz-tools/blob/main/packages/grz-cli/README.md> for details

Install miniforge at /mnt/storage2/MVH/tools/miniforge3/

	> curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
	> bash Miniforge3-$(uname)-$(uname -m).sh

Install GRZ-CLI

	> /mnt/storage2/MVH/tools/miniforge3/bin/conda create -n grz-tools -c conda-forge -c bioconda "grz-cli"
	> /mnt/storage2/MVH/tools/miniforge3/bin/conda activate grz-tools

Updates with:

	> /mnt/storage2/MVH/tools/miniforge3/bin/conda update -n base -c conda-forge conda -c bioconda
	> /mnt/storage2/MVH/tools/miniforge3/bin/conda update -n grz-tools -c conda-forge -c bioconda grz-cli

List all package versions:

	> /mnt/storage2/MVH/tools/miniforge3/bin/conda list -n grz-tools

## Installation of consent mapper

Install mapper (outside of UKT network):

	> git clone https://github.com/KohlbacherLab/mii_broad_consent_mapper.git
	> ./gradlew build 

Note: now a JAR is available as well: https://github.com/orgs/KohlbacherLab/packages?repo_name=mii_broad_consent_mapper


# Updates

## 26.10.2026 - general update

- grz-cli 1.7.2 > 3.0.1

## 22.06.2026 - downgrade to versions used in official GRZ workflow

- fastp: 1.3.3 > 1.0.1
- fastplong: 0.2.2 > 0.2.2

## 05.06.2026 - general update

- fastp: 0.24 > 1.3.3
- fastplong: 0.2.2 > 0.4.1
- grz-cli 1.5.1 > 1.7.2
- GRZ_QC_Workflow v1.3.0 > v1.3.2
