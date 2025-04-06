import sys
import os

import matplotlib

import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import pyvista as pv
import glob
matplotlib.use('Agg')

params = {'text.usetex' : False,
          'font.size'   : 12}
plt.rcParams.update(params)

from scipy.interpolate import griddata

assert len(sys.argv) == 3, "You need to enter the path to the solution!"
path = sys.argv[1]

example = sys.argv[2]

def exact_solution_u_ex1(t,x):
    return np.sin(t) * (1 - x) * x
def exact_solution_v_ex1(t,x):
    return np.cos(t) * (1 - x) * x

def exact_solution_u_ex2(t,x):
    return np.sin(np.pi * t) * np.sin(np.pi * x)
def exact_solution_v_ex2(t,x):
    return np.pi * np.cos(np.pi * t) * np.sin(np.pi * x)


def exact_solution_u_ex4(t,x):
    return np.sin(t) * np.sin(np.pi * x)
def exact_solution_v_ex4(t,x):
    return np.cos(t) * np.sin(np.pi * x)



# list of all VTK-Dateien sorted by time
vtk_files = sorted(glob.glob(path + "*.vtk"))
coordinates_t = []         # time

# Arrays for data
all_solutions_u = []  # u
all_solutions_v = []  # v
coordinates_x = None       # spatial coordinates

for file in vtk_files:
    # load die VTK-file
    data = pv.read(file)
    
    unique_points, indices = np.unique(data.points, axis=0, return_index=True)
    # extract spatial coordinate (1D)
    if coordinates_x is None:
        coordinates_x = data.points[:, 0][indices]
    
    coordinates_t.append(data["TIME"][0])
    

    # extract solution for u and v
    u = data['displacement'][indices] if 'displacement' in data.array_names else None
    v = data['velocity'][indices] if 'velocity' in data.array_names else None

    if u is not None:
        all_solutions_u.append(u[:, 0])
    if v is not None:
        all_solutions_v.append(v[:, 0])


solution_u = np.array(all_solutions_u)  # Shape: (time, space)
solution_v = np.array(all_solutions_v)  # Shape: (time, space)
coordinates_t = np.array(coordinates_t)

if coordinates_t.shape[0] > 1:

    # space-time solution vectors
    x_min, x_max = np.min(coordinates_x), np.max(coordinates_x)
    t_min, t_max = np.min(coordinates_t), np.max(coordinates_t)
    coordinates = np.vstack(
        (
            np.tensordot(coordinates_t, np.ones_like(coordinates_x), 0).flatten(),
            np.tensordot(np.ones_like(coordinates_t), coordinates_x, 0).flatten(),
        )
    ).T


    n_dofs = {"space": coordinates_x.shape[0], "time": coordinates_t.shape[0]}

    # grid for interpolation
    grid_t, grid_x = np.mgrid[t_min:t_max:200j, x_min:x_max:200j]

    # deformation
    values_u = griddata(coordinates, solution_u.flatten(), (grid_t, grid_x), method='linear')
    # velocity
    values_v = griddata(coordinates, solution_v.flatten(), (grid_t, grid_x), method='linear')


    # Check if interpolation produced valid results
    if values_u is None or values_v is None:
        print("Interpolation failed. Please check the input data.")
        sys.exit(1)

    # Fixing axis scaling issues (if necessary)
    values_u[np.isnan(values_u)] = 0  # Replace NaNs with 0
    values_v[np.isnan(values_v)] = 0
    if example == "1":
        # computing exact solution on grid
        exact_u = exact_solution_u_ex1(grid_t, grid_x)
        exact_v = exact_solution_v_ex1(grid_t, grid_x)
        # error
        error_u = np.abs(values_u - exact_u)
        error_v = np.abs(values_v - exact_v)

        titles = [r"$u_h$", r"$u$", r"$\Vert u_h - u\Vert $", r"$v_h$",  r"$v$", r"$\Vert v_h - v\Vert $"]
        data = [values_u, exact_u, error_u, values_v, exact_v, error_v]
        # plotting
        fig, axes = plt.subplots(2, 3, figsize=(10, 12))
    elif example == "2":
        # computing exact solution on grid
        exact_u = exact_solution_u_ex2(grid_t, grid_x)
        exact_v = exact_solution_v_ex2(grid_t, grid_x)
        # error
        error_u = np.abs(values_u - exact_u)
        error_v = np.abs(values_v - exact_v)

        titles = [r"$u_h$", r"$u$", r"$\Vert u_h - u\Vert $", r"$v_h$",  r"$v$", r"$\Vert v_h - v\Vert $"]
        data = [values_u, exact_u, error_u, values_v, exact_v, error_v]

        # plotting
        fig, axes = plt.subplots(2, 3, figsize=(10, 12))
    elif example == "4":
        # computing exact solution on grid
        exact_u = exact_solution_u_ex4(grid_t, grid_x)
        exact_v = exact_solution_v_ex4(grid_t, grid_x)
        # error
        error_u = np.abs(values_u - exact_u)
        error_v = np.abs(values_v - exact_v)

        titles = [r"$u_h$", r"$u$", r"$\Vert u_h - u\Vert $", r"$v_h$",  r"$v$", r"$\Vert v_h - v\Vert $"]
        data = [values_u, exact_u, error_u, values_v, exact_v, error_v]
        # plotting
        fig, axes = plt.subplots(2, 3, figsize=(10, 12))
    else:
        titles = [r"$u_h$", r"$v_h$"]
        data = [values_u, values_v]
        
        # plotting
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))



    for ax, title, data in zip(axes.flat, titles, data):
        im = ax.imshow(data, origin="lower", aspect='auto', extent=[x_min, x_max, t_min, t_max])
        ax.set_title(title)
        fig.colorbar(im, ax=ax)
        ax.set_xlabel('$x$')
        ax.set_ylabel('$t$')


    # Adjusting layout and saving plot
    plt.tight_layout()
    plt.savefig(path + 'solution_scaled.pdf', format='pdf', bbox_inches='tight')

    # 3D plotting
    fig = plt.figure(figsize=(14, 6))

    # 3D plot for displacements
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.plot_surface(grid_x, grid_t, values_u, cmap='viridis')
    ax1.set_title('Displacement')
    ax1.set_xlabel('$x$')
    ax1.set_ylabel('$t$')
    ax1.set_zlabel('$u$')

    # 3D plot for velocity
    ax2 = fig.add_subplot(122, projection='3d')
    ax2.plot_surface(grid_x, grid_t, values_v, cmap='viridis')
    ax2.set_title('Velocity')
    ax2.set_xlabel('$x$')
    ax2.set_ylabel('$t$')
    ax2.set_zlabel('$v$')

    # saving plot
    plt.savefig(path + 'solution_3d.svg', format='svg', bbox_inches='tight')

    # Create 3D plot for displacement (Deformation)
    fig_u = go.Figure(data=[go.Surface(z=values_u, x=grid_x, y=grid_t, colorscale='Viridis')])
    fig_u.update_layout(title='Displacement', scene=dict(
                        xaxis_title='x',
                        yaxis_title='t',
                        zaxis_title='u'),
                        autosize=False, width=800, height=800)

    # Save displacement plot as HTML
    fig_u.write_html(path + 'displacement_3d.html')

    # Create 3D plot for velocity
    fig_v = go.Figure(data=[go.Surface(z=values_v, x=grid_x, y=grid_t, colorscale='Viridis')])
    fig_v.update_layout(title='Velocity', scene=dict(
                        xaxis_title='x',
                        yaxis_title='t',
                        zaxis_title='v'),
                        autosize=False, width=800, height=800)

    # Save velocity plot as HTML
    fig_v.write_html(path + 'velocity_3d.html')




    # Create pyvista grid for VTK output
    grid = pv.StructuredGrid(grid_x, grid_t, np.zeros_like(grid_x))

    # Add displacement and velocity as point data
    grid.point_data["Displacement"] = values_u.flatten()
    grid.point_data["Velocity"] = values_v.flatten()

    # Save the grid to a VTK file
    vtk_file_path = path + "solution_data.vtk"
    grid.save(vtk_file_path)


