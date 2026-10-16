"""Tests for center_map. Need network on first run (maps are cached in .pytest_cache).

gemmi is used only here, as an independent reference implementation.
"""
import numpy as np
import pytest
from pymol import cmd, CmdException

from center_map import center_map  # importing registers the PyMOL command

gemmi = pytest.importorskip('gemmi')

SITE = '3tj3 and chain A and resi 100-110'

ENTRIES = {
    '3tj3': 'P 21 21 21, model sticks out of the cell',
    '2vb1': 'P 1, triclinic',
    '1ejg': 'P 21, monoclinic',
    '1tqn': 'I 2 2 2, centred',
}


@pytest.fixture(scope='session')
def fetch_dir(pytestconfig):
    return str(pytestconfig.cache.mkdir('pdbe'))


def fetch(fetch_dir, code, *types):
    cmd.reinitialize()
    cmd.set('fetch_path', fetch_dir)
    for t in ('cif',) + types:
        if cmd.fetch(code, type=t, async_=0) == -1:
            pytest.skip('cannot fetch %s %s' % (code, t))


def as_gemmi(name, tmp_path):
    path = str(tmp_path / (name + '.ccp4'))
    cmd.save(path, name)
    return gemmi.read_ccp4_map(path)


def box_xyz(m):
    """Cell sampling, box start and box data of an exported map, in x,y,z order."""
    N = np.array([m.header_i32(8 + i) for i in range(3)])
    start = np.rint(np.array(m.get_extent().minimum.tolist()) * N).astype(int)
    m.setup(float('nan'), gemmi.MapSetup.ReorderOnly)
    return N, start, np.array(m.grid, copy=True)


def far_corner(ref, code):
    """Pseudoatom at lattice point (2, -3, 1): a box around it wraps on every axis."""
    corner = ref.grid.unit_cell.orthogonalize(gemmi.Fractional(2, -3, 1))
    cmd.pseudoatom('spot', pos=corner.tolist())
    cmd.set_symmetry('spot', *cmd.get_symmetry(code))
    return 'spot'


@pytest.mark.parametrize('code', ENTRIES)
@pytest.mark.parametrize('where', ['model', 'far_corner'])
def test_values_identical_to_gemmi_periodic_lookup(code, where, fetch_dir, tmp_path):
    fetch(fetch_dir, code, '2fofc')
    ref = as_gemmi(code + '_2fofc', tmp_path)
    ref.setup(float('nan'))
    assert not np.isnan(np.array(ref.grid)).any()
    selection = far_corner(ref, code) if where == 'far_corner' else code

    center_map(code + '_2fofc', selection)

    N, start, ours = box_xyz(as_gemmi(code + '_2fofc', tmp_path))
    assert ours.shape != tuple(N) or start.any(), 'map was not moved'
    expected = ref.grid.get_subarray(start.tolist(), list(ours.shape))
    assert np.array_equal(ours, expected)  # bit-identical, so no renormalization either
    assert cmd.get_setting_boolean('normalize_ccp4_maps')


@pytest.mark.parametrize('code', ENTRIES)
def test_recentering_an_already_centered_map(code, fetch_dir, tmp_path):
    fetch(fetch_dir, code, '2fofc')
    ref = as_gemmi(code + '_2fofc', tmp_path)
    ref.setup(float('nan'))
    # a box wider than the cell on every axis, no longer starting at the origin
    center_map(code + '_2fofc', 'first (%s and q > 0)' % code, 0.55 * max(cmd.get_symmetry(code)[:3]))
    N, start, data = box_xyz(as_gemmi(code + '_2fofc', tmp_path))
    assert start.any() and (np.array(data.shape) >= N).all()

    center_map(code + '_2fofc', far_corner(ref, code))

    _, start, ours = box_xyz(as_gemmi(code + '_2fofc', tmp_path))
    assert np.array_equal(ours, ref.grid.get_subarray(start.tolist(), list(ours.shape)))


def test_centering_again_inside_a_centered_map(fetch_dir, tmp_path):
    fetch(fetch_dir, '3tj3', '2fofc')
    ref = as_gemmi('3tj3_2fofc', tmp_path)
    ref.setup(float('nan'))
    center_map('3tj3_2fofc', '3tj3', 5)  # box narrower than the cell along x and z

    center_map('3tj3_2fofc', SITE, 3)

    _, start, ours = box_xyz(as_gemmi('3tj3_2fofc', tmp_path))
    assert np.array_equal(ours, ref.grid.get_subarray(start.tolist(), list(ours.shape)))


@pytest.mark.parametrize('code', ENTRIES)
def test_box_covers_every_atom_plus_margin(code, fetch_dir, tmp_path):
    fetch(fetch_dir, code, '2fofc')
    margin = 3.0
    center_map(code + '_2fofc', code, margin)

    m = as_gemmi(code + '_2fofc', tmp_path)
    cell = m.grid.unit_cell
    N, start, data = box_xyz(m)
    frac = np.array([cell.fractionalize(gemmi.Position(*p)).tolist()
                     for p in cmd.get_coords(code)])
    rec = cell.reciprocal()
    pad = margin * np.array([rec.a, rec.b, rec.c])
    assert ((frac - pad) * N >= start).all()
    assert ((frac + pad) * N <= start + np.array(data.shape) - 1).all()


@pytest.mark.parametrize('code', ENTRIES)
def test_pymol_places_centered_map_on_original(code, fetch_dir, tmp_path):
    fetch(fetch_dir, code, '2fofc')
    cmd.copy('orig', code + '_2fofc')
    before = cmd.get_volume_field('orig')

    center_map(code + '_2fofc', code)
    _, start, data = box_xyz(as_gemmi(code + '_2fofc', tmp_path))
    assert data.shape != before.shape or start.any(), 'map was not moved'

    # PyMOL interpolates the centered map at the original (in-cell) grid nodes;
    # exact up to float32 rounding, while a 1-voxel shift would differ by ~1 sigma
    cmd.map_set('orig', 'copy', code + '_2fofc')
    after = cmd.get_volume_field('orig')
    inner = [(np.arange(n) > s) & (np.arange(n) < s + d - 1)
             for n, s, d in zip(before.shape, start, data.shape)]
    mask = np.ix_(*inner)
    assert before[mask].size > 1000
    np.testing.assert_allclose(after[mask], before[mask], rtol=0, atol=2e-3)


def test_default_centers_all_maps_on_all_atoms(fetch_dir):
    fetch(fetch_dir, '3tj3', '2fofc', 'fofc')
    cmd.do('center_map')
    lo, hi = np.array(cmd.get_extent('3tj3'))
    for name in ('3tj3_2fofc', '3tj3_fofc'):
        mlo, mhi = np.array(cmd.get_extent(name))
        assert (mlo <= lo - 5).all() and (mhi >= hi + 5).all(), name


def test_refuses_map_not_covering_full_cell(fetch_dir):
    fetch(fetch_dir, '3tj3', '2fofc')
    inside = '3tj3 and x > 20 and y > 10 and z > 10'  # map_trim ignores atoms off the map
    assert cmd.count_atoms(inside)
    cmd.map_trim('3tj3_2fofc', inside, 2.0)
    assert cmd.get_volume_field('3tj3_2fofc').size < 64 * 80 * 108
    with pytest.raises(CmdException, match='unit cell'):
        center_map('3tj3_2fofc', '3tj3')


def test_refuses_map_from_another_cell(fetch_dir):
    fetch(fetch_dir, '3tj3')
    if cmd.fetch('2vb1', type='2fofc', async_=0) == -1:
        pytest.skip('cannot fetch 2vb1 2fofc')
    with pytest.raises(CmdException, match='cell'):
        center_map('2vb1_2fofc', '3tj3')


def test_refuses_map_computed_in_pymol(fetch_dir):
    fetch(fetch_dir, '3tj3')
    cmd.map_new('gmap', 'gaussian', 1.0, SITE, 2)
    assert 'gmap' in cmd.get_names()
    with pytest.raises(CmdException, match='cell'):
        center_map('gmap', '3tj3')
