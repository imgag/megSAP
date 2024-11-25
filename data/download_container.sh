#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`dirname $(pwd)`
CONTAINER_FOLDER="$root/data/tools/apptainer_container/"

# Path to your settings file
SETTINGS_FILE=$root/settings.ini

if [ ! -f "$SETTINGS_FILE" ]; then
    SETTINGS_FILE="$root/settings.ini.default"
fi

#Ignore this - used for local installation
#CONTAINER_FOLDER=/mnt/storage2/megSAP/tools/apptainer_container

BASE_URL="https://megsap.de/download/container/"

# make sure container folder exists
mkdir -p $CONTAINER_FOLDER
cd $CONTAINER_FOLDER

# Scan the settings file for lines that contain container definitions, ignoring container_folder
grep -E "^container_" "$SETTINGS_FILE" | grep -v "container_folder" | while IFS=' = ' read -r container version; do
    # Extract the toolname from container_<toolname>
    toolname=$(echo "$container" | sed 's/container_//')
    url="${BASE_URL}/${toolname}_${version}.sif"

    # Construct the apptainer pull command
    echo "Downloading ${toolname} version ${version} from $url"
    wget --no-check-certificate -O "${toolname}_${version}.sif" "$url"
done

# Print a message indicating the file has been written
echo "All containers found in megSAPs settings.ini have been downloaded."