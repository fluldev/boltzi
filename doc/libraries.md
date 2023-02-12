# Libraries

Libraries, required to build *boltzi*:

- [tetgen](https://github.com/libigl/tetgen)
- [Armadillo](http://arma.sourceforge.net/download.html)
- [OBJ-Loader](https://github.com/Bly7/OBJ-Loader)
- [yaml-cpp](https://github.com/jbeder/yaml-cpp)
- [exprtk](https://github.com/ArashPartow/exprtk)
Optional:
- [Tables-and-Graphs](https://github.com/tdulcet/Tables-and-Graphs) needed for `-DBOLTZI_DASHBOARD`
  compile as libtablengraph.a containing tables.o and graphs.o

Headerfiles and compiled static libraries should be put into `/usr/local/lib/`
and `/usr/local/include/` respectively in order for the compiler/linker to find them.
Some of them might be readily installable through the package manager of your distribution
or come with build instructions. Check how they are used in CMakeLists.txt if the library
does not contain building instructions.

Packages used in *python* scripts:

- [numpy](https://pypi.org/project/numpy/)
- [scipy](https://pypi.org/project/scipy/)
- [matplotlib](https://pypi.org/project/matplotlib/)
- [pyvista](https://pypi.org/project/pyvista/)
- [SciencePlots](https://pypi.org/project/SciencePlots/)

All of them can be installed with `pip` by typing:
`python3 -m pip install -r doc/requirements.txt`
from the git-repository root directory.
