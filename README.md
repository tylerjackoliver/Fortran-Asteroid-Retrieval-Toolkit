# Asteroid Retrieval

Global repository for all asteroid retrieval-related codes. The following sub-folders exist:

---

* **asteroid_optimiser**: Codes required to optimise the condition of one asteroid; main codes used to optimise these asteroids.
* cost_function: Not currently used. Programs to identify the lowest-velocity transfers for a given asteroid (reduces the amount of Lamberts to compute.) No longer used.
* optimizer_matlab: MATLAB version of the asteroid_optimiser routines for testing in fmincon. Untested.
* orbit_construction: MATLAB routines used to generate orbits in the CR3BP. Not refactored.
* Pre-filter: Contains all files necessary to perform the pre-filtering routines for *one* asteroid in the database. Some versions of the code have been parallelized, either using shared memory paradigms or by GPU computing. These versions are as yet untested.
* utility_functions: Utility functions: various other routines that are helpful.

