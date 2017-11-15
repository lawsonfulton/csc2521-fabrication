import subprocess
import tempfile
import time

from solid import scad_render_to_file

open_scad_path = '/Applications/OpenSCAD.app/Contents/MacOS/OpenSCAD'

def save_to_scad_and_stl(geometry, basename, n_segments=48):
    start = time.time()

    scad_file = basename + '.scad'
    scad_render_to_file(geometry, scad_file, file_header='$fn = %s;' % n_segments)

    print('Converting to STL file:', basename)
    proc = subprocess.Popen(open_scad_path + ' -o "%s" "%s"' % (basename, scad_file), shell=True)
    while proc.poll() == None:
        time.sleep(0.1)

    print('Done. Took ', round(time.time() - start, 2), 's', sep='')
