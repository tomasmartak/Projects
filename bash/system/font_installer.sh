#!/bin/bash

## init variables
# colours
RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
ENDCOLOR="\e[0m"
# vars
input1="$1"
input2=${2:-"truetype"} # defaults to truetype
mapfile -t matched_output < <(ls /usr/share/fonts/ | grep -i "$input2")
unique_ouputs=${#matched_output[@]}

## Main
if [ "$EUID" -ne 0 ]; then
    echo -e "${RED}Please run as root to run the script.${ENDCOLOR}"
    exit 1
else
    # Check if the script was run with an argument
    if [ $# -ne 2 ]; then
        echo -e "Error: no input was provided.\nUsage: $0 <origin dir> <font type>.\n"
        echo -e "Choose one of the following available directories:"
        for dir in "${matched_output[@]}"; do
            echo -e "${YELLOW}$(basename "$dir")${ENDCOLOR}"
        done
        exit 1
    fi

    # Check whether there is a non-ambiguous outdir
    if [ $unique_ouputs -gt 1 ]; then
        echo -e "Error: multiple output directories matched (${matched_output[@]}), try and be more precise."
    else
        outdir="/usr/share/fonts/${matched_output}"
    fi

    # move font folder to outdir
    if [ -d "${input1}" ] && [ -d "${outdir}" ]; then
        sudo mv $input1 $outdir
        echo -e "Installed font: ${GREEN}$(basename "$input1")${ENDCOLOR}\nUpdating cache...\n"
        sudo fc-cache -f -v > /dev/null 2>&1
        echo "Done!"
    else
        # If no matching remote is found, throw an error
        echo -e "Error: Incorrect input (${RED}${input1}${ENDCOLOR}) or target (${RED}${outdir}${ENDCOLOR}) directory."
        exit 1
    fi
fi