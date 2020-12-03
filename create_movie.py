import pyvista as pv
from glob import glob
from natsort import natsorted
from tqdm import tqdm

def create_movie(disloc_files, solutes_files, output_name):
    pv.set_plot_theme("document")
    plotter = pv.Plotter()
    plotter.open_movie(output_name)
    plotter.remove_scalar_bar()
    display = True
    plotter.view_isometric()
    plotter.view_xz()
    plotter.add_axes(interactive=True)

    for disloc_f, solutes_f in tqdm(zip(disloc_files, solutes_files)):
        # import pdb; pdb.set_trace()
        plotter.remove_actor(['disloc', 'solute'])

        dislocation = pv.read(disloc_f)
        solutes = pv.read(solutes_f)
        plotter.add_mesh(dislocation,
                         name = 'disloc',
                         show_scalar_bar=False,
                         render_lines_as_tubes=True,
                         line_width=10)
        # plotter.add_mesh(solutes, show_scalar_bar=False,name = 'solute',)
        plotter.add_mesh(solutes, 
                        # name = 'solute', 
                        show_scalar_bar=False, 
                        render_points_as_spheres=True,
                        point_size=20,)
        if display:
            print('Orient the view, then press "q" to close window and produce movie')
            plotter.show(auto_close=False)
            display = False
        plotter.write_frame()
    plotter.close()

def discover_files():
    dislocations = "dislocations_c.*.vtk"
    solutes = "pointdefects_c.*.vtk"
    return natsorted(glob(dislocations)), natsorted(glob(solutes))

files = discover_files()
create_movie(*files, 'movie.mp4')