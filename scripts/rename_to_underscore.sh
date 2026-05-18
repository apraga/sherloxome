#!/bin/bash
set -euo pipefail

DIR="${1:-data/exp_pro}"

find "$DIR" -depth | while IFS= read -r path; do
    parent=$(dirname "$path")
    base=$(basename "$path")
    new_base="${base//-/_}"
    if [[ "$base" != "$new_base" ]]; then
        mv -v "$parent/$base" "$parent/$new_base"
    fi
done
