import numpy as np


def load_ply(fname):
    """Loads ascii ply files. The standard probably allows/disallows
    things that are (not) accepted here, but it works well enough for
    my purposes, and it's rather fast."""
    # Accept file names and file handles
    if isinstance(fname, basestring):
        fh = open(fname)
    elif hasattr(fname, 'readline'):
        fh = fname
    else:
        raise ValueError('fname must be string or file handle')

    # This is a ply file
    assert fh.readline().strip() == "ply"
    assert fh.readline().strip() == "format ascii 1.0"
    # Read header
    while True:
        l = fh.readline().strip()
        if l.startswith('comment'):
            continue
        elif l.startswith('end_header'):
            break
        elif l.startswith('element vertex'):
            n_vertex = int(l.split()[-1])
        elif l.startswith('element face'):
            n_faces = int(l.split()[-1])
    # Load vertex
    vertex = np.fromfile(fh, dtype=np.float,
                         count=3*n_vertex, sep=' ').reshape((n_vertex, 3))
    # Load faces, the first column should contain nothing but 3's, but
    # I don't check it
    faces = np.fromfile(fh, dtype=np.int,
                        count=4*n_faces, sep=' ').reshape((n_faces, 4))[:, 1:]
    fh.close()
    return vertex, faces

if __name__ == '__main__':
    v, f = load_ply('Armadillo_ascii_converted_using_meshlab.ply')
