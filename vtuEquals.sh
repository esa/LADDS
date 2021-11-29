#!/bin/bash

# needs exactly two arguments
if [[ $# -ne 2 ]]
then
    echo "Invalid number of arguments!"
    echo "Usage: $0 fileA fileB"
    exit -1
fi

fileA=$1
fileB=$2

MAXERROR="10^-8"

# Simple check: Compare number of lines
fileALines=$(wc -l < "${fileA}")
fileBLines=$(wc -l < "${fileB}")
if [[ "${fileALines}" -ne  "${fileBLines}" ]];
then
    echo "Number of lines differ:"
    printf "  %-40s: %5s\n" ${fileA} ${fileALines}
    printf "  %-40s: %5s\n" ${fileB} ${fileBLines}
    exit 1
fi

# Compare all lines
diffLines=$(diff --side-by-side --suppress-common-lines --expand-tabs "${fileA}" "${fileB}")

# if the diff is empty the files are exactly equal
if [[ -z "$diffLines" ]]
then
    exit 0
fi

# check all differing lines and compare them within some margin of error
while IFS= read -r line
do
    FLOAT="[0-9.-]\+"
    differenecs=$(sed --quiet -e "s/ *\($FLOAT\) *\($FLOAT\) *\($FLOAT\) *| *\($FLOAT\) *\($FLOAT\) *\($FLOAT\).*/\1 - \4\n\2 - \5\n\3 - \6/p" <<< $line)
    while IFS= read -r equation
    do
        result=$(bc -l <<< "
        define abs(x) {if (x<0) {return -x}; return x;}
        abs($equation) > $MAXERROR
        ")
        if [[ $result -eq 1 ]]
        then 
            echo "Lines don't match!"
            echo $line
            echo "|$equation| > $MAXERROR"
            exit 2
        fi
    done <<< "${differenecs}"
done <<< "${diffLines}"
