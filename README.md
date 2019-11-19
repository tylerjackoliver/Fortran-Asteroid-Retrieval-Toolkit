# Asteroid Retrieval

Global repository for all asteroid retrieval-related codes. The following sub-folders exist:

---

* asteroid_prefilter: Contains all files necessary to perform the pre-filtering routines for *one* asteroid in the database. Some versions of the code have been parallelized, either using shared memory paradigms or by GPU computing. These versions are as yet untested.
* orbit_construction: Codes required to generate the initial manifold conditions
* asteroid_optimiser: Codes required to optimise the condition of one asteroid
