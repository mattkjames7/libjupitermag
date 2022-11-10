#!/bin/bash
CWD=$(pwd)
declare -a libgits=("https://github.com/mattkjames7/libinternalfield.git" \
                    "https://github.com/mattkjames7/libcon2020.git" \
                    "https://github.com/mattkjames7/libspline.git")
declare -a libdirs=("lib/libinternalfield" "lib/libcon2020" "lib/libspline" )
for i in "${!libdirs[@]}"
do
    echo $i ${libdirs[$i]} ${libgits[$i]}
    cd ${libdirs[$i]}
    git stash
    git pull origin main
    cd ${CWD}
done