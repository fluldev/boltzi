# Sample Invariants

1. Surface normals must point inside the sample.
2. Material faces must be somehow declared in the model and cause *tetgen* to
  not create tetrahedra which inhabit multiple materials. This is needed
  anyway since phonons get scattered when penetrating material boundaries.
3. Cell sizes must be neglectable in comparsion to the change of deviational temperature
   and resulting temperature gradients since temperatures are often just
   taken at the cell centers for performance.
4. Surface has no holes.
