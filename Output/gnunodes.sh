#!/bin/bash

cd "$(dirname "$0")"
if [ ! -f "$1" ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]
then 
    echo  -e "\e[1mUsage\e[0m: gnunodes.sh \e[3m[INPUT FILE] [OUTPUT FILE]\e[0m"
else
    if [ "$#" -lt "2" ]
    then
        OUT="output"
    else
        OUT="$2"
    fi
    echo Writing to "$OUT.png".
    gnuplot -e "filename='$1';"  -c ./gnunodes.plt
    #xdg-open $OUT
fi

