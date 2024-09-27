#!/bin/bash

for a in $@ ; do
        echo "$a"
if [[ "$a" == *".xyz" ]] ; then
        echo -e "\tsplitting multicomponent xyz file"
        split -l $(echo "$(head -n 1 "$a") + 2" | bc) --numeric-suffixes --suffix-length=4 --additional-suffix=".xyz" "$a" ${a%%.*}_
        echo -e "\tcomplete" ; else
        echo -e "\tnot an xyz file - run on multiple structure .xyz file \n\tskipping file"
        continue
fi
done
