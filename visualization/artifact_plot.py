#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from packaging.version import parse 
import pyvista as pv 
from scipy.stats import linregress 
import argparse

if __name__ == "__main__":
    plt.style.use("science")
    parser = argparse.ArgumentParser()

    parser.add_argument("sample_model", type=str,
            help=".obj file of the sample geometry.")
    parser.add_argument("artifact_or_save", type=str,
            help="Artifact .boltzi_artifact or .boltzi_state simulation state file.")
    parser.add_argument("-n", "--isosurfaces", type=int, default=10,
            help="Amount of iso-surfaces of the temperature."),
    parser.add_argument("-l", "--line", type=float, nargs=7, help="Query Temperature along line from point a to b with resolution n.", action="append")
    parser.add_argument("-p", "--plane", type=float, nargs=10, help="Query Temperature on plane defined by surface center, normal, width, height and resolution nx,ny.")
    parser.add_argument("--solid", action="store_true", help="Show solid translucent block instead of isosurfaces.")
    parser.add_argument("--no_fit", action="store_true", help="Do not clip plotted volume.")
    parser.add_argument("--scatter", action="store_true")
    parser.add_argument("--linfit", action="store_true", help="Create linear fit for the line plot.")
    parser.add_argument("--export", type=str, default="")
    parser.add_argument("--scale", type=float, default=1)
    
    args = parser.parse_args()

    fname_obj = args.sample_model
    fname_artifact = args.artifact_or_save
    
    with open(fname_artifact, "r") as f:
        data_name = next(f).replace("\n","")
        if data_name != "simulation_state":
            data_time = float(next(f))
            n_cells = list(map(int,next(f).split(",")))
            if data_name in ["temperature","meanfreepath","meanfreetime","simulatedphonons"]:
                data = np.array([[
                    *list(map(float, k.removesuffix("\n").split(","))), 
                    *list(map(float, v.removesuffix("\n").split(",")))
                ] for k,v in zip(f,f)])
            elif data_name == "heatcurrent":
                data = []
                for k,v1,*v in zip(f,f,f,f):
                    data.append([
                        *list(map(float, k.removesuffix("\n").split(","))), 
                        *list(map(float, v1.removesuffix("\n").replace("\t","").split(","))),
                        *list(map(lambda x: float(x.removesuffix("\n").removeprefix("\t")), v))
                    ])
                    next(f)
                data = np.array(data)
            else:
                print(f"Unknown artifact \"{data_name}\"")
                exit(-1)
        # data format 0,1,2 = x,y,z data point and 3-n are the values at this data point
        else:
            data = np.array([[
                *list(map(float,l1.removesuffix("\n").split(","))), 
                *list(map(float,l2.removesuffix("\n").split(","))),
                *list(map(float,l2.removesuffix("\n").split(","))) 
                ]for l1, l2, l3 in zip(f,f,f)], dtype=float)

    # clear data of nan
    data = data[~np.any(np.isnan(data), axis=-1)]

    pvdata = pv.PolyData(data[:,:3])

    reader = pv.get_reader(fname_obj)
    mesh = reader.read()
    lower = mesh.bounds[::2]
    higher = mesh.bounds[1::2]
    
    plttr = pv.Plotter()
    plttr.add_mesh(mesh, culling="back")
    
    if data_name != "simulation_state":
        n_cells = np.array(n_cells)
        grid = pv.UniformGrid()
        grid.origin = lower
        grid.spacing = list((np.array(higher)-np.array(lower))/n_cells)
        grid.dimensions = list(n_cells+1)
    
    if data_name in ["temperature","meanfreetime","meanfreepath","simulatedphonons"]:
        pvdata[f"{data_name} t={data_time:.3e}"] = data[:,4]
        interp = grid.interpolate(pvdata, radius=1, sharpness=10, strategy='mask_points')
        # interp = interp.translate(list(np.array(grid.spacing)/2))
        if args.scatter:
            plttr.add_mesh(pvdata, render_points_as_spheres=True, point_size=10, cmap="hot", scalar_bar_args={"color" : "black"})
        else:
            if not args.no_fit:
                interp = interp.clip_surface(mesh, invert=False)
            if args.solid:
                plttr.add_mesh_clip_plane(interp, opacity=1, cmap="hot", scalar_bar_args={"color" : "black"})
            else:
                plttr.add_mesh(interp.contour(args.isosurfaces), opacity=.5, cmap="hot", scalar_bar_args={"color" : "black"})
        if args.line:
            for line in args.line:
                spls = interp.sample_over_line(line[:3], line[3:6], resolution=int(line[6]))
                fig, ax = plt.subplots(1,1)
                xs = np.linalg.norm(spls.points-np.array(line[:3]), axis=1)
                ys = np.array(spls[f"{data_name} t={data_time:.3e}"])
                mask = np.isnan(ys)
                xs = xs[~mask]
                ys = ys[~mask]
                ax.plot(xs, ys, c="black", label="$T$ simuliert")
                if args.linfit:
                    fit = linregress(xs, ys)
                    print(f"Slope: {fit.slope/args.scale}")
                    ax.plot(xs, fit.intercept + fit.slope*xs, c="red", ls="--", label=f"$ax+b, a={fit.slope/args.scale:.2e},~b={fit.intercept:.2e}$")
                    ax.legend()
                if args.export:
                    np.save(args.export, np.vstack([xs, ys]))
                plt.show()

                plttr.add_mesh(pv.Line(line[:3], line[3:6]), color="black")
        if args.plane is not None:
            center = np.array(args.plane[:3])
            normal = np.array(args.plane[3:6])
            plane = pv.Plane(center, normal, *args.plane[6:8], *list(map(int,args.plane[8:])))
            interp_s = plane.interpolate(pvdata, radius=1)
            if not args.no_fit:
                interp_s = interp_s.clip_surface(mesh, invert=False)
            plttr.add_mesh(interp_s, cmap="hot")
            plt.scatter(interp_s.points[:,0], interp_s.points[:,1], c=interp_s[f"{data_name} t={data_time:.3e}"], cmap="hot")
            plt.show()

    elif data_name == "heatcurrent":
        norms = np.linalg.norm(data[:,4:7],axis=1)
        print(f"Mean heatcurrent: {np.mean(data[:,4:7],axis=0)}")
        pvdata = pv.PolyData(data[:,:3])
        pvdata[f"vectors"] = data[:,4:7] / norms[:, np.newaxis]
        pvdata[f"{data_name} t={data_time:.3e}"] = norms
        plttr.add_mesh(
            pvdata.glyph(orient=f"vectors", scale=f"{data_name} t={data_time:.3e}", factor=min(grid.spacing)/np.max(norms), geom=pv.Arrow()),
            cmap="hot", scalar_bar_args={"color" : "black"}
        )
        #stream, src = mesh.streamlines(
        #    f"{data_name} t={data_time:.3e}", return_source=True, terminal_speed=0.0, n_points=200, source_radius=0.1
        #)
        #plttr.add_mesh(
        #    stream.tube(radius=0.0015),            
        #    #cmap="hot", scalar_bar_args={"color" : "black"}
        #)
    elif data_name == "simulation_state":
        print(f"n Phonons: {len(data)}")
        nppts = data[:,3:6]/ args.scale 
        mask = nppts[:,0] >= mesh.bounds[0] 
        mask &= nppts[:,0] <= mesh.bounds[1] 
        mask &= nppts[:,1] >= mesh.bounds[2] 
        mask &= nppts[:,1] <= mesh.bounds[3] 
        mask &= nppts[:,2] >= mesh.bounds[4] 
        mask &= nppts[:,2] <= mesh.bounds[5] 

        nppts = nppts[mask]
        pts = pv.PolyData(nppts)
        fig, axs = plt.subplots(3,1)
        axs[0].hist2d(nppts[:,0], nppts[:,1], bins=(10,10), density=True)
        axs[1].hist2d(nppts[:,1], nppts[:,2], bins=(10,10), density=True)
        axs[2].hist2d(nppts[:,2], nppts[:,0], bins=(10,10), density=True)
        plt.show()
        pts["wavevector magnitude * sign"] = np.linalg.norm(data[mask,:3],axis=-1) * data[mask,8]
        plttr.add_mesh(pts, cmap="hot", render_points_as_spheres=True, point_size=5, scalar_bar_args={"color" : "black"})
        plt.hist(np.linalg.norm(data[mask,:3], axis=1), bins=100, density=True)
        plt.show()

        dirs = data[mask,:3] / np.linalg.norm(data[mask,:3],axis=1)[:,np.newaxis]
        phis = np.arctan2(dirs[:,1], dirs[:,2])
        thetas = np.arccos(np.abs(dirs[:,0]))
        fig, (axphi,axtheta) = plt.subplots(1,2,subplot_kw={"projection":"polar"})
        axphi.hist(phis,bins=30, density=True)
        axphi.plot(np.linspace(0,2*np.pi,1000), np.ones(1000)/(2*np.pi), c="red")
        axtheta.hist(thetas,bins=30, density=True)
        pv.PolyData(dirs[::5]).plot()
        xs = np.linspace(0,np.pi/2,1000)
        axtheta.plot(xs, np.sin(xs)*np.cos(xs)*2, c="red")
        plt.show()

    plttr.set_background("white")
    plttr.show_grid(color="black")
    plttr.show()


