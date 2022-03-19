#!/bin/zsh

#
# Tests a given LADDS exe using a given yaml file with a list of different AutoPas configs.
#

set -euo pipefail

if [[ "$#" -ne 2 ]]
then
    echo "usage: $0 LADDS CFG"
    exit 1
fi

# make sure arguments are absolute paths
EXE_ABS=$(realpath "$1")
CFG_ABS=$(realpath "$2")
CFG=$(basename "${CFG_ABS}")
CSV_OUT="$PWD/testConfigs.csv"

# configure here what to test. Format= CONTAINER:TRAVERSAL
CONFIGS=(
    "LinkedCells:lc_c08"
    "LinkedCells:lc_c04_HCP"
    "VarVerletListsAsBuild:vvl_as_built"
    "VerletListsCells:vlc_sliced_balanced"
    "VerletClusterLists:vcl_c06"
)

echo "Total[ns],Collision[ns],Integrator[ns],Update[ns],Container,Traversal" #> "${CSV_OUT}"

# for all configs create a dir, cfg, and output
for config in "${CONFIGS[@]}"
do
    CONTAINER=${config%%:*}
    TRAVERSAL=${config#*:}
    mkdir -p ${CONTAINER}-${TRAVERSAL} && cd ${CONTAINER}-${TRAVERSAL}

    # copy and adapt config
    cp "${CFG_ABS}" .
    sed --in-place --regexp-extended "s/(\s+Container\s*:).*/\1 ${CONTAINER}/" ${CFG}
    sed --in-place --regexp-extended "s/(\s+Traversal\s*:).*/\1 ${TRAVERSAL}/" ${CFG}

    # run the sim
    "${EXE_ABS}" "${CFG}" > output.out

    # collect interesting info in the csv
    echo \
$(sed --quiet --regexp-extended "s/\s*Total\s*:\s*([0-9]+) ns.*/\1/p" output.out),\
$(sed --quiet --regexp-extended "s/\s*Collision detection\s*:\s*([0-9]+) ns.*/\1/p" output.out),\
$(sed --quiet --regexp-extended "s/\s*Integrator\s*:\s*([0-9]+) ns.*/\1/p" output.out),\
$(sed --quiet --regexp-extended "s/\s*Container update\s*:\s*([0-9]+) ns.*/\1/p" output.out),\
${CONTAINER},${TRAVERSAL}

    cd ..
done

