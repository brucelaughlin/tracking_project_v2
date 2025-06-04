#!/bin/bash

# Accidentally modified this, not sure if will still work

symLinkDir="$2"

wildCardArg="$3"

for fileSourceDir in "${1}/"*
do
    date=$( grep -o "[0-9]\{8\}$" <<< $fileSourceDir )
    inputFileName=( "${fileSourceDir}/out/"*"${wildCardArg}"* )
    outputFile="${symLinkDir}/${date}_$(basename $inputFileName)"
    ln -s $(realpath "$inputFileName") $outputFile
done



