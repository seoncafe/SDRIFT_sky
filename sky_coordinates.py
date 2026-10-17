"""Fixed equatorial-axis / Galactic rotations shared with the browser.

Matrix generated from Astropy's ICRS -> Galactic Cartesian basis transform.
The stored geocentric apparent directions are rotated in the same fixed axes;
this does not introduce another time-dependent aberration correction.
"""
import numpy as np

COORD_SYSTEMS = ('equatorial', 'galactic')
EQ_TO_GAL = np.array([
    [-0.054875657712591654, -0.8734370519556157, -0.48383507361671546],
    [0.49410943719272676, -0.4448297212232951, 0.7469821839866675],
    [-0.8676661375596576, -0.19807633727300059, 0.4559838136873016],
])


def check_coords(coords):
    if coords not in COORD_SYSTEMS:
        raise ValueError('Coordinates must be equatorial or galactic')
    return coords


def from_equatorial(vectors, coords='equatorial'):
    check_coords(coords)
    v = np.asarray(vectors, dtype=float)
    return v @ EQ_TO_GAL.T if coords == 'galactic' else v


def to_equatorial(vectors, coords='equatorial'):
    check_coords(coords)
    v = np.asarray(vectors, dtype=float)
    return v @ EQ_TO_GAL if coords == 'galactic' else v
