# This script parses an arbitrary number of files each containing the CLI output of LADDS,
# extracts the timer information and merges it into one CSV for further processing.

if [[ $# -lt 1 ]]
then
    echo "Usage: " $0 " outFiles..."
    exit 1
fi

outFiles="$@"

getFilenames () {
    echo "Filename"
    for f in ${outFiles}
    do
        echo $(basename $f)
    done

}

# get all lines about number of migrants (sed), calculate the average (awk) and print it rounded.
getMigrants () {
    echo "Avg Migrants"
    for f in ${outFiles}
    do
        sed --quiet 's/.*Migrated particles: \([0-9]\+\) .*/\1/p' ${f} \
            | awk '{ total += $1 } END {printf "%.0f\n", total/NR}'
    done
}

# Create a column of the tag and the corresponding values from all out files
get () {
    local tag="$@"
    echo $tag
    for f in ${outFiles}
    do
        # after the tag only spaces are allowed
        local result=$(sed --quiet "s/.*${tag}\s*:.*( *\([0-9.]\+\)s).*/\1/p" "${f}")
        # alternative regex for values not in ( )
        if [[ -z "$result" ]]
        then
            local result=$(sed --quiet "s/.*${tag}\s*: *\([^ ]\+\).*/\1/p" "${f}")
        fi
        # error output
        if [[ -z "$result" ]]
        then
            echo "Couldn't find $tag in $f"
        else
            echo $result
        fi
    done
}

# combine all info into one csv
paste -d ','                                \
    <(getFilenames)                         \
    <(get MPI Ranks)                        \
    <(get OpenMP Threads per Rank)          \
    <(get Container)                        \
    <(get desiredCellsPerDimension)         \
    <(get 'Total (ranks accumulated)')      \
    <(get Initialization)                   \
    <(get Simulation)                       \
    <(get Integrator)                       \
    <(get Resolving Burn ups)               \
    <(get Constellation insertion)          \
    <(get Communication)                    \
    <(get Collision detection)              \
    <(get Collision detection immigrants)   \
    <(get Collision detection emmigrants)   \
    <(get Collision writer)                 \
    <(get Evasion writer)                   \
    <(get Container update)                 \
    <(get Output)                           \
    <(get 'Total (wall-time)')              \
    <(get One iteration)                    \
    <(getMigrants)                          \

