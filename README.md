# center-map

`center_map` recenters crystallographic density maps on the model, the way
the old EDS/RCSB DSN6 maps were. It is a command for PyMOL.

Maps from `fetch 3tj3, type=2fofc` (PDBe) cover one unit cell starting at
the origin, so much of the model sits outside them. `center_map` cuts a new
box around the selection out of the periodic crystal.

## Install

Into the Python environment the viewer runs in:

    pip install center-map

then add one line to `~/.pymolrc`:

    import center_map

Without pip (e.g. bundled installers): *Plugin > Plugin Manager > Install
New Plugin* and choose `center_map.py`. The only dependency is numpy,
which is already bundled.

## Use

    fetch 3tj3
    fetch 3tj3, type=2fofc
    center_map                          # all maps, all atoms, 5 A margin
    isomesh mesh, 3tj3_2fofc, 1.0

    center_map 3tj3_2fofc, chain A      # again, on one chain of one map

`center_map` replaces the map object in place; make meshes after it. A
centered map can be centered again on atoms inside its box; for a bigger
box, fetch the map again. `help center_map` lists the arguments.

## Why the density is exact

* Crystal density is periodic, so a map covering one full unit cell already
  holds every value: grid point *i* is grid point *i* mod *N*. The box is
  filled by copying grid values; nothing is interpolated or resampled, also
  for non-orthogonal cells.
* The original contour scale is kept. Maps are normally rescaled on load to
  the sigma of the loaded grid, which for a box around the protein differs
  from the sigma of the cell; that is switched off while the box is loaded.
* It refuses maps it cannot handle exactly: maps that do not cover a full
  unit cell (e.g. asymmetric-unit maps from CCP4 `fft`), and maps whose cell
  differs from the selection's (cryo-EM maps, maps made with `map_new`,
  other crystal forms).

## Limitations

One map state; object matrices (moved objects) are ignored.

## Tests

    pip install pytest gemmi
    python -m pytest tests      # downloads a few PDBe maps on the first run

gemmi is used only by the tests, as an independent reference: the boxes
must match gemmi's periodic lookup bit for bit (P 21 21 21, P 1, P 21 and
I 2 2 2 entries), land on the original map when the viewer places them,
and cover every atom plus the margin.

## License

MIT. PyMOL is a trademark of Schrödinger, LLC. This project is not
affiliated with or endorsed by Schrödinger.
