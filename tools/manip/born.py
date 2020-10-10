"""Calculate solvation energies and forces using Born model."""
import argparse
import logging
import numpy as np
from pdb2pqr.io import read_pqr


_LOGGER = logging.getLogger()


def build_parser():
    """Build an argument parser.

    :returns:  argument parser
    :rtype:  argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        description=(
            "Calculate solvation energies and forces using the Born model"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--verbose", action="store_true", help="provide per-atom information"
    )
    parser.add_argument(
        "--forces", action="store_true",
        help="provide forces as well as energies"
    )
    parser.add_argument(
        "--solvent-dielectric", help="solvent dielectric constant",
        default=78.54, type=float
    )
    parser.add_argument("pqr_path", help="path to PQR-format molecule file")
    return parser


def distance2_matrix(coords):
    """Generate a matrix of distances squared.

    :param np.ndarray coords:  coords
    :returns:  matrix of distances squared
    :rtype:  np.ndarray
    """
    matrix = []
    natom = coords.shape[0]
    for iatom in range(natom):
        coord1 = coords[iatom, :]
        disp2 = np.square(coord1 - coords)
        dist2 = np.sum(disp2, axis=1)
        matrix.append(dist2)
    matrix = np.array(matrix)
    return matrix


def scale_distance2(distances2, radii):
    """Generate a matrix of distances squared scaled by radii.

    :param np.ndarray distances2:  matrix of distances squared between atoms
        (e.g., as returned by :func:`distance2_matrix`)
    :param np.ndarray radii:  vector of radii for atoms
    :returns:  matrix of distances squared scaled by radii
    :rtype:  np.ndarray
    """
    return distances2 / radii[:, None] / radii[None, :]


def main():
    """Main driver function."""
    parser = build_parser()
    args = parser.parse_args()
    _LOGGER.info("Verbose:  %s", args.verbose)
    _LOGGER.info("Calculate forces:  %s", args.forces)
    _LOGGER.info("Solvent dielectric:  %s", args.solvent_dielectric)
    _LOGGER.info("PQR path:  %s", args.pqr_path)
    with open(args.pqr_path, "r") as pqr_file:
        atom_list = read_pqr(pqr_file)
    coords = []
    charges = []
    radii = []
    natom = len(atom_list)
    for atom in atom_list:
        coords.append([atom.x, atom.y, atom.z])
        charges.append(atom.charge)
        radii.append(atom.radius)
    coords = np.array(coords)
    charges = np.array(charges)
    radii = np.array(radii)
    print(radii)
    distances2 = distance2_matrix(coords)
    print(distances2)
    scaled2 = scale_distance2(distances2, radii)
    print(scaled2)


if __name__ == "__main__":
    main()
