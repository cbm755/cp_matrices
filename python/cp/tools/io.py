import numpy as np


def load_ply(fname):
    """Loads ascii and binary ply files.

    Input
    -----
    fname : file name or file handle

    Output
    ------
    vertices : array of shape (npoints, 3)
               The coordinates of each vertex. Has dtype float64
    faces : array of shape (nfaces, 3)
            Contains the vertices' indices of each face. That is,
            vertices[faces[i]] gives you the coordinates of the
            vertices in the i-th face. Has dtype int32

    Notes
    -----
    The standard (dis)allows things that are (not) accepted here, but
    it works well enough for my purposes, and it's rather fast (even
    more reading binary).

    If your file is not correctly loaded, open it in MeshLab, and save
    it as a .ply file (binary is faster, more compact, ascii is easier
    to read, seems to lose some precision), and check *None* of the
    additional parameters.

    More information about the format:
    http://paulbourke.net/dataformats/ply/"""
    
    # Accept file names and file handles
    if isinstance(fname, basestring):
        fh = open(fname)
    elif hasattr(fname, 'readline'):
        fh = fname
    else:
        raise ValueError('fname must be string or file handle')

    def readline():
        l = fh.readline().strip()
        while l.startswith('comment'):
            l = fh.readline().strip()
        return l
        
    # This is a ply file
    try:
        assert readline() == "ply"
    except AssertionError:
        raise ValueError("Not a ply file")
    # Figure out the format
    file_format = readline().split()
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
            #endiannes = '<'
            vertex_dtype = np.dtype('<f4')
            faces_dtype = np.dtype('<u1, 3<i4')
        elif file_format[1] == 'binary_big_endian':
            #endiannes = '>'
            vertex_dtype = np.dtype('>f4')
            faces_dtype = np.dtype('>u1, 3>i4')
        else:
            raise ValueError("Format not understood")
        
        #name2dtype = {'char': '{0}i1'.format(endiannes),
        #              'uchar': '{0}u1'.format(endiannes),
        #              'short': '{0}i2'.format(endiannes),
        #              'ushort': '{0}u2'.format(endiannes),
        #              'int': '{0}i4'.format(endiannes),
        #              'uint': '{0}u4'.format(endiannes),
        #              'float': '{0}f4'.format(endiannes),
        #              'double': '{0}f8'.format(endiannes),
        #              }

    # Read header
    while True:
        l = readline()
        if l.startswith('comment'):
            continue
        elif l.startswith('end_header'):
            break
        elif l.startswith('element vertex'):
            n_vertex = int(l.split()[-1])
            for i in xrange(3):
                assert readline().startswith('property float')
        elif l.startswith('element face'):
            n_faces = int(l.split()[-1])
            assert readline().startswith('property list uchar int')

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
    if not vertex.dtype.isnative:
        vertex = vertex.byteswap().newbyteorder()
    if not faces.dtype.isnative:
        faces = faces.byteswap().newbyteorder()

    return vertex.astype(np.float64), faces.astype(np.int32)

if __name__ == '__main__':
    v, f = load_ply('Armadillo_ascii_converted_using_meshlab.ply')
