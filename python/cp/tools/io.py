import numpy as np


def load_ply(fname):
    """Loads ascii and binary ply files. The standard probably
    allows/disallows things that are (not) accepted here, but it works
    well enough for my purposes, and it's rather fast (even more
    reading binary)."""
    # Accept file names and file handles
    if isinstance(fname, basestring):
        fh = open(fname)
    elif hasattr(fname, 'readline'):
        fh = fname
    else:
        raise ValueError('fname must be string or file handle')

    # This is a ply file
    try:
        assert fh.readline().strip() == "ply"
    except AssertionError:
        raise ValueError("Not a ply file")
    # Figure out the format
    file_format = fh.readline().strip().split()
    assert file_format[0] == "format" and file_format[2] == "1.0"
    if file_format[1] == 'ascii':
        binary = False
        sep = ' '
        faces_multiplier = 4
        vertex_dtype = np.float32
        faces_dtype = np.int32
    else:
        binary = True
        sep = ''
        faces_multiplier = 1
        if file_format[1] == 'binary_little_endian':
            # This is not tested (haven't found a file with this
            # format), but it should work anyway
            vertex_dtype = np.dtype('<f4')
            faces_dtype = np.dtype('<u2, 3<i4')
        elif file_format[1] == 'binary_big_endian':
            vertex_dtype = np.dtype('>f4')
            faces_dtype = np.dtype('>u2, 3>i4')
        else:
            raise ValueError("Format not understood")
    
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
    # np.fromfile is very low level, and doesn't accept gzipped
    # files. If needed, np.fromstring(fh.read(... could be used
    vertex = np.fromfile(fh, dtype=vertex_dtype,
                         count=3*n_vertex, sep=sep).reshape((n_vertex, 3))
    # Load faces, the first column should contain nothing but 3's, but
    # I don't check it
    faces = np.fromfile(fh, dtype=faces_dtype, count=faces_multiplier * n_faces,
                        sep=sep)
    fh.close()
    if binary:
        faces = faces['f1']  # Second field
    else:
        faces = faces.reshape((n_faces, 4))[:, 1:]
    return vertex, faces

if __name__ == '__main__':
    v, f = load_ply('Armadillo_ascii_converted_using_meshlab.ply')
