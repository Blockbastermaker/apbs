"""Calculate solvation energies and forces using Born model."""
import argparse
import logging
import numpy as np
from pdb2pqr.io import read_pqr


_LOGGER = logging.getLogger()
ELECTRON_CHARGE = 1.6021773e-19
AVOGADRO_NUMBER = 6.0221367e+23
VACUUM_PERMIT = 8.8541878e-12
SCALING = (
    (1e-3)*(1e10)*ELECTRON_CHARGE*ELECTRON_CHARGE*AVOGADRO_NUMBER/(
        4*np.pi*VACUUM_PERMIT
    )
)


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


def born_distances(coords, radii):
    """Generate a matrix of effective Born distances.

    :param np.ndarray coords:  array of atomic coordinates
    :param np.ndarray radii:  array of atomic radii
    :returns:  matrix of effective Born distances.
    :rtype:  np.ndarray
    """
    distance2 = distance2_matrix(coords)
    scaled2 = scale_distance2(distance2, radii)
    radii_mat = np.outer(radii, radii)
    distance2 = distance2 + np.multiply(radii_mat, np.exp(-0.25*scaled2))
    return np.sqrt(distance2)


def born_distance(pos1, pos2, radius1, radius2):
    """Generate an effective Born distance.

    :param np.ndarray pos1:  cooordinates of first atom
    :param np.ndarray pos2:  coordinates of second atom
    :param float radius1:  radius of first atom
    :param float radius2:  radius of second atom
    :returns:  effective Born distance
    :rtype:  float
    """
    disp = pos1 - pos2
    dist2 = np.dot(disp, disp)
    dist2 = dist2 + radius1*radius2*np.exp(-0.25*dist2/radius1/radius2)
    return np.sqrt(dist2)


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
    # distances = born_distances(coords, radii)
    # charges = np.outer(charges, charges)
    # energies = np.divide(charges, distances)
    # np.fill_diagonal(energies, 0.0)
    # energies = 0.5*SCALING*energies
    # print(energies)
    # print(np.sum(energies, axis=0))
    for iatom in range(natom):
        charge1 = charges[iatom]
        radius1 = radii[iatom]
        pos1 = coords[iatom, :]
        energy = 0
        for jatom in range(iatom+1, natom):
            charge2 = charges[jatom]
            radius2 = radii[jatom]
            pos2 = coords[jatom, :]
            print(pos1, pos2)
            dist = born_distance(pos1, pos2, radius1, radius2)
            energy += SCALING*charge1*charge2/dist
            print(energy)


if __name__ == "__main__":
    main()
