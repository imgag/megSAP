#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`

# Path to your settings file
SETTINGS_FILE=$root/settings.ini

#Ignore this - used for local installation
#CONTAINER_FOLDER=/mnt/storage2/megSAP/tools/apptainer_container

# Directory containing the containers (change in settings.ini if you want to pull containers elsewhere)
CONTAINER_FOLDER=$(grep "^container_folder" "$SETTINGS_FILE" | awk -F' = ' '{print $2}')
BASE_URL="https://download.imgag.de/public/megSAP_container"

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
    wget --no-check-certificate --no-proxy -O "${toolname}_${version}.sif" "$url"
done

# Print a message indicating the file has been written
echo "All containers found in megSAPs settings.ini have been downloaded."