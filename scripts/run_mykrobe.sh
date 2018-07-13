#!/usr/bin/env sh
barcode="$1"
species="$2"
input=$(python3 -c "import os.path; print(os.path.abspath('$3'))")
output=$(python3 -c "import os.path; print(os.path.abspath('$4'))")
container=$(python3 -c "import os.path; print(os.path.abspath('$5'))")

mkdir -p data/mykrobe/"$barcode"
cd data/mykrobe/"$barcode" || return 1
mykrobe predict "$barcode" "$species" \
  --ont --seq "$input" --output "$output"
