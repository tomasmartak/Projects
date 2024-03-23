#!/bin/bash

## init variables
# colours
RED="\e[31m"
GREEN="\e[32m"
ENDCOLOR="\e[0m"
# vars
rclone_remotes=$(rclone listremotes)
input="$1"
matched_remote=$(echo "$rclone_remotes" | grep -i "$input")
mount_directory="$HOME/Documents/GDrive/${matched_remote%:}"

# Check if the script was run with an argument
if [ $# -ne 1 ]; then
    echo -e "Error: no input was provided.\nUsage: $0 <rclone remote_name>"
    echo -e "Choose one of the following available remotes:"
    for remote in "${rclone_remotes[@]}"; do
        echo -e "${RED}${remote}${ENDCOLOR}"
    done
    echo -e ""
    exit 1
fi

# Use rclone to list existing remotes and grep for the target remote
if [ -n "$matched_remote" ]; then
    # Check if the directory already exists; if not, create
    if [ ! -d "$mount_directory" ]; then
        mkdir -p "$mount_directory"
    fi
    
    # Mount the matching remote if it exists
    rclone mount --daemon --vfs-cache-mode minimal "$matched_remote" "${mount_directory}"
    echo -e "Mounted remote: ${GREEN}${matched_remote}${ENDCOLOR}"
else
    # If no matching remote is found, throw an error
    echo -e "Error: Remote ${RED}${input}${ENDCOLOR} not found."
    exit 1
fi