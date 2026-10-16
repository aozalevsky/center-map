"""center_map: recenter crystallographic maps on a selection.

    pip install center-map
    echo "import center_map" >> ~/.pymolrc
"""
import os
import tempfile

import numpy as np
from pymol import cmd, CmdException


def center_map(map='', selection='all', margin=5.0, quiet=1, _self=cmd):
    '''
DESCRIPTION

    Recenter crystallographic map(s) on a selection, like the old EDS/RCSB
    DSN6 maps did. Maps from "fetch code, type=2fofc" cover one unit cell
    starting at the origin, so the model often sits partly outside them.

    The map must cover a full unit cell of the selection's crystal, or
    already contain the new box (e.g. centering again on fewer atoms).
    Values are copied from the periodic lattice, not interpolated, and keep
    the contour scale of the original map.

USAGE

    center_map [ map [, selection [, margin ]]]

ARGUMENTS

    map = str: map object {default: all maps}

    selection = str: atoms the map should cover {default: all}

    margin = float: padding around the atoms, in Angstrom {default: 5}

EXAMPLE

    fetch 3tj3
    fetch 3tj3, type=2fofc
    center_map
    isomesh mesh, 3tj3_2fofc, 1.0
    '''
    margin, quiet = float(margin), int(quiet)
    xyz = _self.get_coords(selection)
    sym = _self.get_symmetry(selection)
    if xyz is None or not sym:
        raise CmdException('"%s" has no atoms or no unit cell' % selection)
    for name in [map] if map else _self.get_names_of_type('object:map'):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'map.ccp4')
            _self.save(path, name)
            with open(path, 'rb') as f:
                raw = f.read()
        h = np.frombuffer(raw, np.int32, 256).copy()  # CCP4 header words
        cell = h[10:16].view(np.float32)
        if not np.allclose(cell, sym[:6], rtol=1e-3):
            raise CmdException('%s: unit cell differs from "%s"' % (name, selection))
        n, start, axes = h[0:3], h[4:7], h[16:19] - 1  # column, row, section
        N = h[7:10][axes]  # cell sampling along the file axes
        data = np.frombuffer(raw, np.float32, int(n.prod()), 1024 + h[23]).reshape(n[::-1])

        # box in grid units around the atoms; margin along each axis is
        # margin * |reciprocal axis|, i.e. measured perpendicular to the cell faces
        a, b, c = cell[:3]
        ca, cb, cg = np.cos(np.radians(cell[3:]))
        sg = np.sqrt(1 - cg * cg)
        v = np.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg)
        frac = np.linalg.inv([[a, b * cg, c * cb],
                              [0, b * sg, c * (ca - cb * cg) / sg],
                              [0, 0, c * v / sg]])[axes]
        f = xyz @ frac.T
        pad = margin * np.linalg.norm(frac, axis=1)
        lo = np.floor((f.min(0) - pad) * N).astype(np.int32)
        hi = np.ceil((f.max(0) + pad) * N).astype(np.int32) + 1

        # crystal density is periodic: grid point i is grid point i mod N; that
        # needs a full cell along an axis, unless the box lies inside the map
        if h[3] != 2 or not ((n >= N) | ((lo >= start) & (hi <= start + n))).all():
            raise CmdException('%s does not cover a full unit cell, nor the box '
                               'around "%s"' % (name, selection))
        idx = [(np.arange(l, u) - s) % m for l, u, s, m in zip(lo, hi, start, N)]
        box = data[np.ix_(*idx[::-1])]

        h[0:3], h[4:7], h[23] = hi - lo, lo, 0
        tmpname = _self.get_unused_name('_center_map')
        normalize = _self.get('normalize_ccp4_maps')
        _self.set('normalize_ccp4_maps', 0)  # sigma of the box != sigma of the cell
        try:
            _self.load_raw(h.tobytes() + box.tobytes(), 'ccp4', tmpname, quiet=1)
        finally:
            _self.set('normalize_ccp4_maps', normalize)
        _self.delete(name)
        _self.set_name(tmpname, name)
        if not quiet:
            print(' center_map: %s covers "%s" + %g A' % (name, selection, margin))


def __init_plugin__(app=None):
    pass  # lets the Plugin Manager install this file; the command is registered on import


cmd.extend('center_map', center_map)
cmd.auto_arg[0]['center_map'] = cmd.auto_arg[0]['map_double']
cmd.auto_arg[1]['center_map'] = cmd.auto_arg[1]['map_trim']
