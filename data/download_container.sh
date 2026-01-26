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
mkdir -p checksums

BASE_URL="https://megsap.de/download/container/"

# Scan the settings file for lines that contain container definitions, ignoring container_folder
grep -E "^container_" "$SETTINGS_FILE" | grep -v "container_folder" | while IFS=' = ' read -r container version; do
    # Extract the toolname from container_<toolname>
    toolname=$(echo "$container" | sed 's/container_//')
    filename="${toolname}_${version}.sif"
    url="${BASE_URL}/${filename}"
    remote_md5_url="${BASE_URL}/checksums/${filename}.md5"
    local_md5_file="checksums/${filename}.md5"

    needs_download=true

    # Always try to get the remote checksum
    if wget --no-check-certificate -q -O "${filename}.md5" "$remote_md5_url"; then
        # If local container exists
        if [ -f "$filename" ]; then
            echo "Found existing container: $filename"

            # Create local checksum if missing
            if [ ! -f "$local_md5_file" ]; then
                echo "Local checksum missing. Creating $local_md5_file"
                md5sum -b "$filename" | cut -d ' ' -f1 > "$local_md5_file"
            fi

            # Compare checksums
            if diff "${filename}.md5" "$local_md5_file" > /dev/null; then
                echo "MD5 sums match: ${filename} is up to date. Skipping download."
                needs_download=false
                rm -f "${filename}.md5"
            else
                echo "MD5 sums differ: Redownloading ${filename} (updated version on server)"
                rm -f "$filename" "$local_md5_file"
            fi
        fi
    else
        echo "WARNING: MD5 sum file for container ${filename} missing on megsap.de. Could not validate container."
        needs_download=true
    fi

    # Download if needed
    if [ "$needs_download" = true ]; then
        echo "Downloading ${toolname} version ${version} from $url"
        wget --no-check-certificate -O "$filename" "$url"
        chmod 775 "$filename" || true

        # Recreate checksum
        md5sum -b "$filename" | cut -d ' ' -f1 > "$local_md5_file"

        if [ -f "${filename}.md5" ]; then
            if diff "${filename}.md5" "$local_md5_file" > /dev/null; then
                echo "MD5 sums match after download: ${filename} successfully downloaded."
            else
                echo "WARNING: MD5 sums differ after download -> possible corruption in ${filename}!"
            fi
            rm -f "${filename}.md5"
        fi
    fi
done

# Print a message indicating the file has been written
echo "All containers found in megSAPs settings.ini have been downloaded."