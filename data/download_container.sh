#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`dirname $(pwd)`

# Path to your settings file
SETTINGS_FILE=$root/settings.ini

if [ ! -f "$SETTINGS_FILE" ]; then
    SETTINGS_FILE="$root/settings.ini.default"
fi

# Extract container_folder from the settings file
CONTAINER_FOLDER=$(grep -E "^container_folder" "$SETTINGS_FILE" | awk -F ' = ' '{print $2}' | sed "s|\[path\]|$root|")

# Verify the container folder was found
if [ -z "$CONTAINER_FOLDER" ]; then
    echo "Error: container_folder not found in $SETTINGS_FILE"
    exit 1
fi

# Ensure the container folder exists
mkdir -p "$CONTAINER_FOLDER"
cd "$CONTAINER_FOLDER"

BASE_URL="https://megsap.de/download/container/"

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