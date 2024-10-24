#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`

# Path to your settings file
SETTINGS_FILE=$root/settings.ini

#Ignore this - used for local installation
#CONTAINER_FOLDER=/mnt/storage2/megSAP/tools/singularity_container

# Directory containing the containers (change in settings.ini if you want to pull containers elsewhere)
CONTAINER_FOLDER=$(grep "^container_folder" "$SETTINGS_FILE" | awk -F' = ' '{print $2}')

# make sure container folder exists
mkdir -p $CONTAINER_FOLDER
cd $CONTAINER_FOLDER

# File to store the pull commands //TODO remove when actual pull commands are executed
PULL_COMMANDS_FILE="$CONTAINER_FOLDER/pull_commands.txt"

# Clear the file if it already exists //TODO remove when actual pull commands are executed
> "$PULL_COMMANDS_FILE"

# Scan the settings file for lines that contain container definitions, ignoring container_folder
grep -E "^container_" "$SETTINGS_FILE" | grep -v "container_folder" | while IFS=' = ' read -r container version; do
    # Extract the toolname from container_<toolname>
    toolname=$(echo "$container" | sed 's/container_//')

    # Construct the apptainer pull command
    echo "apptainer pull docker://imgag/${toolname}:${version}" >> "$PULL_COMMANDS_FILE" #TODO modify when actual pull commands are executed
done

# Print a message indicating the file has been written #TODO change to message indicating all containers have been pulled
echo "Pull commands have been written to $PULL_COMMANDS_FILE"