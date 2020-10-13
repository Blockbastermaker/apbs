"""Test Born model functionality."""
import argparse
from collections import deque
from subprocess import run
from pathlib import Path
import numpy as np
import pandas as pd


BORN_BINARY = "C:\\APBS-3.0.0\\share\\apbs\\tools\\bin\\Release\\born.exe"
DIELECTRIC = 78.54
DATA_PATH = Path("data/born-data_apbs-c.csv")


def icosahedron():
    """Generate points on a tetrahedron.

    Tetrahedron has edge length 2 and circumradius of :func:`np.sqrt`(phi + 2),
    where `phi` is the golden ratio.

    :returns:  vertices of tetrahedron
    :rtype:  np.ndarray
    """
    vertices = []
    phi = 0.5*(1.0 + np.sqrt(5.0))
    x = 0
    for y in (-1, 1):
        for z in (-phi, phi):
            point = deque([x, y, z])
            for i in range(len(point)):
                vertices.append(list(point))
                point.rotate()
    return np.array(vertices)


def dump_pqr(pqr_file, coords, charges, radii):
    """Dump information to PQR-like file.

    :param pqr_file:  file object for PQR data (ready for writing)
    :type pqr_file:  file
    :param coords:  DataFrame of [x, y, z] coordinates
    :type coords:  pd.DataFrame
    :param charges:  DataFrame of charges
    :type charges:  pd.DataFrame
    :param radii:  DataFrame of radii
    :type radii:  pd.DataFrame
    """
    for iatom, coord in coords.iterrows():
        pqr_file.write(
            "ATOM {num} DUM DUM {num} "
            "{coords[0]:4.3f} {coords[1]:4.3f} {coords[2]:4.3f} "
            "{charge} {radius}\n".format(
                num=iatom+1, coords=coord, charge=charges[iatom],
                radius=radii[iatom]
            )
        )


def parse_output(out_file):
    """Parse output from old APBS `born` program.

    :param out_file:  file object ready for reading
    :type out_file:  file
    :returns:  list of results
    :rtype: list(dict)
    """
    results = {}
    for line in out_file:
        words = line.strip().split()
        if words:
            if words[0] == "Atom":
                key = "atom %d" % int(words[1][:-1])
                what = words[2].lower()
                try:
                    value = float(words[4])
                except ValueError:
                    value = None
                if key in results:
                    results[key][what] = value
                else:
                    results[key] = {what: value}
            elif words[0] == "GB":
                key = "total GB"
                what = words[2].lower()
                try:
                    value = float(words[4])
                except ValueError:
                    value = None
                if key in results:
                    results[key][what] = value
                else:
                    results[key] = {what: value}
    results_ = []
    for key, val in results.items():
        val["what"] = key
        results_.append(val)
    return results_


def generate_molecules():
    """Generate molecule data.

    :return:  DataFrame with molecule data
    :rtype:  pd.DataFrame
    """
    vertices = icosahedron()
    radii = deque([0.5, 1.0, 2.0, 3.0])
    charges = deque([-2.0, -1.0, 1.0, 2.0])
    point0 = np.array([0, 0, 0])
    atoms = []
    imol = 0
    for ipoint, point1 in enumerate(vertices):
        for point2 in vertices[ipoint+1:]:
            positions = np.array([point0, point1, point2])
            for iatom in range(3):
                atom = {"molecule": imol, "atom": iatom+1}
                atom["radius"] = list(radii)[iatom]
                atom["charge"] = list(charges)[iatom]
                atom["x"] = positions[iatom][0]
                atom["y"] = positions[iatom][1]
                atom["z"] = positions[iatom][2]
                atoms.append(atom)
            radii.rotate()
            charges.rotate()
            imol += 1
    results = pd.DataFrame(atoms)
    print(results)
    return results


def generate_old_results(dielectric, born_binary, csv_path):
    """Generate results using old APBS C program.

    :param dielectric:  dielectric value
    :type dielectric:  float
    :param born_binary:  path to binary for Born program
    :type born_binary:  str
    :param csv_path:  path to CSV file for writing results
    :type csv_path:  str
    """
    molecules = generate_molecules()
    results = []
    for mol_id, molecule in molecules.groupby("molecule"):
        pqr_path = "ion{:0>2}.pqr".format(mol_id)
        print("Generating PQR file:  %s" % pqr_path)
        with open(pqr_path, "wt") as pqr_file:
            dump_pqr(
                pqr_file, molecule[["x", "y", "z"]], molecule["charge"],
                molecule["radius"]
            )
        out_path = "ion{:0>2}.out".format(mol_id)
        print("Generating output:  %s" % out_path)
        with open(out_path, "wt") as out_file:
            born_cmd = args.apbs_c_binary
            born_args = "-v -f %4.2f %s" % (DIELECTRIC, pqr_path)
            run([born_cmd, born_args], stdout=out_file)
        print("Parsing output:  %s" % out_path)
        with open(out_path, "rt") as out_file:
            for result in parse_output(out_file):
                result["pqr"] = pqr_path
                results.append(result)
    results = pd.DataFrame(results)
    results = results[
        [
            "pqr", "what", "energy", "x-force", "y-force", "z-force",
            "ri-force", "rj-force"
        ]
    ]
    print("Writing results to %s." % DATA_PATH)
    results.to_csv(DATA_PATH, index=False)
    print(results.to_string())


def build_parser():
    """Build an argument parser.

    :return:  argument parser
    :rtype:  argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        description="Generate data files for Born model",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--what", help="which data file to generate", choices=["apbs-c"],
        default="apbs-c"
    )
    parser.add_argument(
        "--apbs-c-binary", help="path to APBS C binary", default=BORN_BINARY
    )
    parser.add_argument(
        "--solvent-dielectric", help="solvent dielectric", default=DIELECTRIC
    )
    parser.add_argument(
        "--csv-path", help="path for CSV output", default=DATA_PATH
    )
    return parser


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    print(args)
    if args.what == "apbs-c":
        generate_old_results(
            dielectric=args.solvent_dielectric, born_binary=args.apbs_c_binary,
            csv_path=args.csv_path
        )
    else:
        raise ValueError("Unknown --what option: %s" % args.what)
