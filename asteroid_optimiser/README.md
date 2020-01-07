## Asteroid optimisers

Main routines are global/ and local/.

* global/ implements a Derivative-Free Boender-Timmer-Rinnoy Kan algorithm to attempt to find the global minimum. The main program used.
* local/ implements the BOBYQA derivative-free local optimiser. Not often used; rarely makes improvements on the global solver.
* ephemeris/ implements the global solver but by using ephemerides at all points of the optimisation. Not refactored.
* global_test/ is a `playground' for testing various things in the optimisation.
