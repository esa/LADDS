# Detailed guide
Satellite constellations (e.g. Starlink, OneWeb) are usually described by a list of orbital shells.
An orbital shell is described by a 4-tuple with information about `altitude`, `inclination`, `number of
planes`, and `number of satellites` per plane. The notebook 
`ConstellationGeneration.ipynb` can be used to generate position and 
velocity vectors for each satellite as well as metadata in a .yaml file. 
The file contains further information needed by LADDS
and it logs the parameters that have been used to create the vectors.

---

## Detailed guide:

### 0. Execute the first cell to import libraries and initialize variables

### 1. Initialize Constellation and reset planet_list and other lists
Provide the constellation name, the time of launch as well as
the time it takes to fully deploy the constellation and execute the cell.
### 2.1 Generate shell = ( altitude, inclination, nPlanes, nSats )
The execution of this cell creates a single shell ( planet_list_tmp ) that is not yet a part of the constellation. It is possible
to visualize the generated satellites by executing the cell `Plot current shell`
before adding it to the constellation. Append the shell to the constellation
by executing cell `2.2`. Each time the code of `2.1` is run, the temporary data is
overridden.

Create a shell by specifying the parameters of an orbital shell which are the 
following:
* `altitude` the altitude of each satellite in kilometers, relative to the
earths surface
* `inclination` the inclination of each orbital plane in degrees
* `nPlanes` the number of orbital planes
* `nSats` the number of satellites in each plane

In some cases it is necessary to specify more parameters: `offsetM`, `argPeriapsis`,
`startingW`, and `W_area`. It is best to modify them only if necessary and to use
the default values otherwise. Refer to the end of the file to see their usage
information (additional parameters).

### 2.2 Append to planet_list and other
Append the temporary set of satellites/planets (planet_list_tmp) to the
constellation (planet_list). The same happens to lists with the shells
meta information.

### 3. Create position and velocity vectors
Execute the code in order to create position and velocity vectors that will
be stored in the 'objects' variable.

### 4. Plot and store results
Prepare .csv and .yaml output and simultaneously plot the constellation by
executing this cell.

### 5. Save Constellation as files in directory
Create a directory with the same name that was specified in `1` in the variable
constellation_name in the `data` directory of the project. It contains three
files: Two .csv files with positions and velocities and a .yaml file with the
metadata of the constellation.

---

### Additional parameters
### OffsetM
`offsetM` offset for the mean anomaly parameter of an orbital element set that is
accumulating (relative phasing).
The mean anomaly of the first satellite of any plane is always bigger by
offsetM than the first satellite of the previous plane. In a Walker Constellation, 
which is a typical model for constellations, the parameter follows the formula with F 
being a natural number in {0,...,nPlanes-1}:

`offsetM = F * 360 / (nPlanes*nSats)` 

Note: F is a natural number: F is element of {0, ... , nPlanes - 1}

Note: the property of a set phasing difference
between neighboring planes gets lost if the shell is added over time



### argPeriapsis
`argPeriapsis` argument of periapsis is another parameter of the orbital element
set. The base to which the Mean Anomaly is added. Choosing an irrational number
prevents satellites from two symmetrical planes to have the same position.

### startingW
A set of planes that span 360° are created by adding the stepsize `360/nPlanes`
the orbital element W (Longitude of ascending node). `startingW` determines
the base this offset is added to. For shells with the same
altitude and the same or similar inclination, a way to prevent
the overlap exists:

The first of the overlapping shells holds the value 0 for startingW.

The other shells should follow the formula:

`(360/G)*(i/N)`

Note: G stands for the smallest common denominator of the nPlanes parameters of
each shell

Note: N is the number of planes overlapping

Note: i is an index in {1,...,N-1}. each shell should have a distinct i

### W_area
this parameter determines which degree the planes of the shell span.

