# Defines

- `BOLTZI_MEASURE_TIMINGS` - Show timings of expensive functions.
- `BOLTZI_DASHBOARD` - Use dashboard with graphs for temperature and amount of
                       phonons during nonlinear DSMC simulation.
- `BOLTZI_ALTERNATIVE_SCATTER_PROBABILITY` - Use small time step approximation
                                            for the selection of phonons to
                                            scatter in order to eliminate
                                            some numerical errors for big time steps.

- `BOLTZI_ALWAYS_RANDOMIZE_THREE_PHONON_SCATTERING` - Randomize the direction
                                of phonons scattered in three-phonon-processes
                                no matter if they are Umklapp-processes or not.

- `BOLTZI_BRANCH_TRANSITIONS` - Allow transitions between branches when scattering.
- `BOLTZI_PHONON_SPAWN_TIMES` - Save spawntimes of phonons when using linear algorithm
                                within file phonon_spawn_times.list
- `BOLTZI_BOUNDARY_DIRS` - Log phonons transitioning boundaries
- `BOLTZI_COLLISION_POINTS` - Log phonon collision points with sample surface if 
                              inelastically scattered
