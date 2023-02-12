# boltzi

## Description

Boltzi is an implementation of the DSMC-Method for the linearized,
non linearized and differential kind of the algorithm.
It can be used for the simulation of the Boltzmann-Transport-Equation of phonons.
A comprehensive description of the algorithm can be found in ^[1].
It was designed for my master thesis, which is not yet published.
This specific implementation aims to support arbitrarily complex sample geometries.

Samples are given in wavefront format (.obj) and can be generated with Blender.
A detailed description and hints for generating samples
with Blender can be found in doc/sample-invariants.md and doc/object_naming_language.md.
Other sample properties are given through a YAML-config-file, described in doc/example_config.yml.

Physical quantities can be measured through "artifacts",
also described in doc/example_config.yml when defined by config
or given as a commandline option to the simulation.
All artifacts are position dependent and recorded in a grid with dimensions
given in the config file or given through the -r NX NY NZ commandline option.
There is also the possibility of recording the heat current through
certain surfaces with the -b FILENAME commandline option which comes
with supplemental information like where the phonons originated.

Scripts for the visualization and evaluation of the artifacts and boundary heat
currents can be found within the directory visualization/ and scripts/.

Details on the commandline options can be found within src/main.cpp.

[1]: PÉRAUD, Jean-Philippe M.; LANDON, Colin D.; HADJICONSTANTINOU, Nicolas G. Monte Carlo methods for solving the Boltzmann transport equation. Annual Review of Heat Transfer, 2014, 17. Jg.

## Building

The compilation of boltzi requires external dependencies given in doc/libraries.md
that might require modification of the CMakeLists.txt
in the root directory to build successfully.
Most of those will be available through your package manager.
Compilation is only possible on GNU/Linux.
I recommend creating a directory called "build" in the root directory and then execute:

```bash
cd build
cmake ..
make
```

The final executable will then be within the build directory and called "boltzi".

Python3 scripts within scripts/ and visualization/ require the libraries
given in doc/libraries.md or requirements.txt,
which can be installed through pip with

```bash
pip3 install -r requirements.txt
```

## Compilation flags

There are several defines that change the behaviour of the simulation
and add tracing information.
A description of those defines can be found in doc/defines.md and
must be set through the appropriate entries in the CMakeLists.txt in the root directory.

## Tests

The majority of the tests in the tests/ directory
is outdated and must be modified to be used again.
