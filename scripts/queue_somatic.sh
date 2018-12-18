#!/bin/bash

set -e -u -o pipefail

# resolve script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# usage
usage() {
    cat <<EOT
Usage: $0 [-map] <tumor sample ID> <control sample ID>

If -map is specified, single sample analysis for tumor and control are queued in addition to somatic analysis.
EOT
}

# arg check
if [[ $# -lt 2 ]]; then
    usage
    exit 2
fi

if [[ $1 == "-map" ]]; then
    shift
    echo "php ${SCRIPT_DIR}/../src/NGS/db_queue_analysis.php -type 'single sample' -samples $1 -args '-steps ma,db -somatic'"
    echo "php ${SCRIPT_DIR}/../src/NGS/db_queue_analysis.php -type 'single sample' -samples $2 -args '-steps ma,db -somatic'"
fi

echo "php ${SCRIPT_DIR}/../src/NGS/db_queue_analysis.php -type somatic -samples $1 $2 -info tumor normal"