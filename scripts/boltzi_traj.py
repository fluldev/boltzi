#!/bin/env python3

import pyvista as pv
import numpy as np

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=str)
    parser.add_argument("traj", type=str)
    parser.add_argument("scale", type=float)
    args = parser.parse_args()

    reader = pv.get_reader(args.file)
    mesh = reader.read()

    traj = np.genfromtxt(args.traj, delimiter=",")

    plttr = pv.Plotter()
    plttr.add_mesh(pv.PolyData(traj[:1] / args.scale), render_points_as_spheres=True, color="red", point_size=10)
    plttr.add_mesh(pv.PolyData(traj[-1:] / args.scale), render_points_as_spheres=True, color="green", point_size=10)
    plttr.add_lines(traj / args.scale)
    plttr.add_mesh(mesh, culling="back")
    plttr.show()

    print(np.linalg.norm(traj[0] - traj[-1]))


