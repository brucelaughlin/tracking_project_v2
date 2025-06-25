#!/bin/bash


# This is for Clark's ROMS output, meaning it works with his naming conventions and storage style (see <inputFileNam>, which is a custom path)

symLinkDir="$2"

wildCardArg="$3"

for fileSourceDir in "${1}/"*
do
    date=$( grep -o "[0-9]\{8\}$" <<< $fileSourceDir )
    inputFileName=( "${fileSourceDir}/out/"*"${wildCardArg}"* )
    outputFile="${symLinkDir}/${date}_$(basename $inputFileName)"
    ln -s $(realpath "$inputFileName") $outputFile
done



