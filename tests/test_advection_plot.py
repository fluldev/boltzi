import numpy as np
import pyvista as pv 

if __name__ == "__main__":
    data = np.genfromtxt("build/test_advection.out.csv", delimiter=",")
    plltr=pv.Plotter()
    plltr.add_mesh(pv.Cube((.5,.5,.5)))
    plltr.add_mesh(pv.lines_from_points(data[:, :3]))
    plltr.show()
