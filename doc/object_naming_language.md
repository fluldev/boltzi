# Object Naming Language

## Basics

Sample properties are defined in .obj (wavefront) files.
Those can be easily created with [blender](https://www.blender.org/) or
other 3D modelling software.
Surface properties and regions within the sample are set by the names
of the objects.
A name of an object is structured like a function call. It consists out
of the property name followed by parameters seperated by commata.
This specific string can be embedded within a different string which
identifies the object in the file.

## Types of Boundaries

- `isothermal_simple(temperature)`
- `vacuum_specular()`
- `rough_bk(std_deviation)`
- `rough_lambert()`
- `rough_even()`
- `measure()`

## Blender Export Parameters
In Blender: File->Export->Wavefront (.obj)
Include:
- [x] OBJ Objects

Transform:
- Scale = 1.0 (scaling should be done in the simulation configuration.)
- Forward = Y Forward
- UP = Z Up

Geometry:
- [x] Write Normals
- [x] Triangulate Faces
