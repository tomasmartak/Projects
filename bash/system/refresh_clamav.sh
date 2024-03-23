#!/bin/bash

# init bash colour variables
RED="\e[31m"
GREEN="\e[32m"
ENDCOLOR="\e[0m"

# code
if [ "$EUID" -ne 0 ]; then
    echo -e "${RED}Please run as root to refresh clamav.${ENDCOLOR}"
    exit 1
else
    echo -e "${GREEN}Refreshing clamav...${ENDCOLOR}"
    systemctl stop clamav-freshclam
    freshclam
    systemctl start clamav-freshclam
    echo -e "${RED}Done!${ENDCOLOR}"
fi



