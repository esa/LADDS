if [[ $# < 1 ]]
then
    echo "Usage: " $0 " outFiles..."
    exit -1
fi

outFiles="$@"

getFilenames () {
    echo "Filename"
    for f in ${outFiles}
    do
        echo $f
    done

}

# Create a column of the tag and the corresponding values from all out files
get () {
    local tag="$@"
    echo $tag
    for f in ${outFiles}
    do
        # after the tag only spaces are allowed
        sed --quiet "s/.*${tag}\s*:.*( *\([0-9.]\+\)s).*/\1/p" "${f}"
    done
}

# combine all info into one csv
paste -d ','                                \
    <(getFilenames)                         \
    <(get 'Total (ranks accumulated)')      \
    <(get Initialization)                   \
    <(get Simulation)                       \
    <(get Integrator)                       \
    <(get Resolving Burn ups)               \
    <(get Constellation insertion)          \
    <(get Communication)                    \
    <(get Collision detection)              \
    <(get Collision detection immigrants)   \
    <(get Collision writer)                 \
    <(get Container update)                 \
    <(get Output)                           \
    <(get 'Total (wall-time)')              \
    <(get One iteration)                    \

