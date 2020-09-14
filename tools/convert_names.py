"""Convert atom and residue naming schemes between programs.

.. todo:: This could be a much stronger tool if it converted structure files
   to and from well-defined formats to resolve residue and atom names.
"""
import logging
from pdb2pqr.pdb import read_pdb, ATOM, HETATM
raise NotImplementedError(
    "This is on hold until PDB2PQR has write_pdb function"
)

_LOGGER = logging.getLogger()


class NameMapper:
    """A class that maps between residue/atom naming schemes."""
    def __init__(self, scheme1_name, scheme2_name):
        """Initialize with scheme names.

        :param scheme1_name:  first naming scheme
        :type scheme1_name:  str
        :param scheme2_name:  second naming scheme
        :type scheme2_name:  str
        """
        self.scheme1_name = scheme1_name
        self.scheme2_name = scheme2_name
        self.mappings = {}

    def add_mapping(self, res1, atom1, res2, atom2):
        """Add a mapping between the naming schemes.

        :param res1: residue name from scheme 1 or None if applies to all
            residues
        :type res1:  str
        :param atom1:  atom name from scheme 1 or None if applies to all atoms
            within the specified residue
        :type atom1:  str
        :param res2: residue name from scheme 2 or None if applies to all
            residues
        :type res2:  str
        :param atom2:  atom name from scheme 2 or None if applies to all atoms
            within the specified residue
        :type atom2:  str
        """
        if res1 in self.mappings:
            self.mappings[res1][atom1] = (res2, atom2)
        else:
            self.mappings[res1] = {atom1: (res2, atom2)}

    def add_mappings(self, mappings):
        """Add multiple mappings between naming schemes.

        :param mappings:  a list of 4-tuples of the form
            ``(res1, atom2, res2, atom2)`` where

            ``res1``
                residue name from scheme 1 or None if applies to all residues
            ``atom1``
                atom name from scheme 1 or None if applies to all atoms within
                the specified residue
            ``res2``
                residue name from scheme 1 or None if applies to all residues
            ``atom2``
                atom name from scheme 2 or None if applies to all atoms within
                the specified residue
        :type mappings:  [(str, str, str, str)]
        """
        for res1, atom1, res2, atom2 in mappings:
            self.add_mapping(res1, atom1, res2, atom2)

    def get_mapping(self, res=None, atom=None):
        """Map from scheme 1 to scheme 2 for the specified residue and atom.

        :param res:  residue name or None if applies to all residues
        :type res:  str
        :param atom:  atom name or None if applies to all atoms within the 
            specified residue
        :returns:  (residue, atom) tuple for mapping
        :rtype:  (str, str)
        :raises KeyError:  if mapping isn't found
        """
        res_dict = self.mappings[res]
        return res_dict[atom]

    def convert_pdb(self, lines):
        """


AMBER_TO_CHARMM = dict(
    [
        ("HIE", "HSE"), 
        ("HID", "HSD"), 
        ("2HB ", "HB1 "), 
        ("3HB ", "HB2 "), 
        ("4HB ", "HB3 "), 
        ("2HG ", "HG1 "), 
        ("3HG ", "HG2 "), 
        ("2HA ", "HA1 "), 
        ("3HA ", "HA2 "), 
        ("2HE ", "HE1 "), 
        ("3HE ", "HE2 "), 
        ("4HE ", "HE3 "), 
        ("2HZ ", "HZ1 "), 
        ("3HZ ", "HZ2 "), 
        ("4HZ ", "HZ3 "), 
        ("2HD ", "HD1 "), 
        ("3HD ", "HD2 "), 
        ("4HD ", "HD3 "), 
        ("2HE2", "HE21"), 
        ("3HE2", "HE22"), 
        ("2HH1", "HH11"), 
        ("3HH1", "HH12"), 
        ("2HH2", "HH21"), 
        ("3HH2", "HH22"), 
        ("2HD2", "HD21"), 
        ("3HD2", "HD22"), 
        ("4HD2", "HD23"), 
        ("2HD1", "HD11"), 
        ("3HD1", "HD12"), 
        ("4HD1", "HD13"), 
        ("2HG2", "HG21"), 
        ("3HG2", "HG22"), 
        ("4HG2", "HG23"), 
        ("2HG1", "HG11"), 
        ("3HG1", "HG12"), 
        ("4HG1", "HG13"), 
        ("HG  SER", "HG1 SER"), 
        (" H ", " HN"), 
    ]
)
CHARMM_TO_AMBER = dict(
    [(c, a) for (a, c) in AMBER_TO_CHARMM.items()]
)
WHATIF_TO_AMBER = []
# Convert all bases to ribo form (i.e., not deoxy)
WHATIF_TO_AMBER += [
    (" U ", " RU "),
    (" G ", " RG "),
    (" C ", " RC "),
    (" A ", " RA ")
]
# Convert U 'H5' to H5''
WHATIF_TO_AMBER += [
    (" H5 RU ", " H5' RU "),
    (" HO2 RU ", " H2' RU ")
]
# phosphate index G
WHATIF_TO_AMBER += [
    (" OP1 RG ", " O1P RG "),
    (" OP2 RG ", " O2P RG "),
]
# (52) Convert G 'H5' to H5'' 
WHATIF_TO_AMBER += [
    (" H5 RG ", " H5' RG "),
    (" HO2 RG ", " H2' RG "),
]
# phosphate index A
WHATIF_TO_AMBER += [
    (" OP1 RA ", " O1P RA "),
    (" OP2 RA ", " O2P RA "),
]
# (119) Convert A 'H5' to H5'' 
WHATIF_TO_AMBER += [
    (" H5 RA ", " H5' RA "),
    (" HO2 RA ", " H2' RA "), 
]
# (233) phosphate index U
WHATIF_TO_AMBER += [
    (" OP1 RU ", " O1P RU "),
    (" OP2 RU ", " O2P RU "),
]
# (420) phosphate index C
WHATIF_TO_AMBER += [
    (" OP1 RC ", " O1P RC "),
    (" OP2 RC ", " O2P RC "),
]
# Convert C 'H5' to H5''
WHATIF_TO_AMBER += [
    (" 'H5' RC ", " H5'' RC "),
    (" 'HO2 RC ", " H2'' RC "),
]
# (48740) 3' end terminal H
WHATIF_TO_AMBER += [
    (" 'HO3 RA ", " H3T TE3 "),
    (" 'HO3 RC ", " H3T TE3 "),
    (" 'HO3 RG ", " H3T TE3 "),
    (" 'HO3 RU ", " H3T TE3 "),
]
# (not sure how to deal with terminal O3')
# (need to find all terminal 'HO3 and modify all matching O3' by hand)
# (48743) protein label changes VAL
WHATIF_TO_AMBER += [
    (" 1HG1 VAL ", " HG11 VAL "),
    (" 2HG1 VAL ", " HG12 VAL "),
    (" 3HG1 VAL ", " HG13 VAL "),
    (" 1HG2 VAL ", " HG21 VAL "),
    (" 2HG2 VAL ", " HG22 VAL "),
    (" 3HG2 VAL ", " HG23 VAL "),
    (" RC VAL ", " C VAL "),
    (" H1 VAL ", " HT1 NTE "),
    (" H2 VAL ", " HT2 NTE "),
    (" H3 VAL ", " HT3 NTE "),
]
# (need to change N VAL to N NTE by hand for all terminal VAL's)
# nonterminal VAL H
WHATIF_TO_AMBER += [
    (" H VAL ", " HN VAL "),
]
# (48761) protein label changes LYS
WHATIF_TO_AMBER += [
    (" RC LYS ", " C LYS "),
    (" H LYS ", " HN LYS "),
    (" 2HB LYS ", " HB1 LYS "),
    (" 3HB LYS ", " HB2 LYS "),
    (" 2HG LYS ", " HG1 LYS "),
    (" 3HG LYS ", " HG2 LYS "),
    (" 2HD LYS ", " HD1 LYS "),
    (" 3HD LYS ", " HD2 LYS "),
    (" 2HE LYS ", " HE1 LYS "),
    (" 3HE LYS ", " HE2 LYS "),
]
# (48783) protein label changes GLU
WHATIF_TO_AMBER += [
    (" RC GLU ", " C GLU "),
    (" H GLU ", " HN GLU "),
    (" 2HB GLU ", " HB1 GLU "),
    (" 3HB GLU ", " HB2 GLU "),
    (" 2HG GLU ", " HG1 GLU "),
    (" 3HG GLU ", " HG2 GLU "),
]
# (48798) protein label changes LEU
WHATIF_TO_AMBER += [
    (" RC LEU ", " C LEU "),
    (" H LEU ", " HN LEU "),
    (" 2HB LEU ", " HB1 LEU "),
    (" 3HB LEU ", " HB2 LEU "),
    (" 1HD1 LEU ", " HD11 LEU "),
    (" 2HD1 LEU ", " HD12 LEU "),
    (" 3HD1 LEU ", " HD13 LEU "),
    (" 1HD2 LEU ", " HD21 LEU "),
    (" 2HD2 LEU ", " HD22 LEU "),
    (" 3HD2 LEU ", " HD23 LEU "),
]
# (48851) protein label changes ALA
WHATIF_TO_AMBER += [
    (" RC ALA ", " C ALA "),
    (" H ALA ", " HN ALA "),
    (" 1HB ALA ", " HB1 ALA "),
    (" 2HB ALA ", " HB2 ALA "),
    (" 3HB ALA ", " HB3 ALA "),
]
# (48861) protein label changes GLY
WHATIF_TO_AMBER += [
    (" RC GLY ", " C GLY "),
    (" H GLY ", " HN GLY "),
    (" 2HA GLY ", " HA1 GLY "),
    (" 3HA GLY ", " HA2 GLY "),
]
# (48882) protein label changes HIS
# tricky because there are 3 protonation states of histidine; amber
# dinstiguishes between them and WHATIF doesn't.  probably have to 
# figure out by hand what's what, and change HIS to HSD/HSE/HSP.
# shortcut (hopefully): change all HIS to HSE, then look for "HD1 HSE"
# which will be a part of any HIS for which this will be incorrect, then
# change only the incorrect ones.
# (already used) (" HIS ", " HSE "),
WHATIF_TO_AMBER += [
    (" RC HSD ", " C HSD "),
    (" H HSD ", " HN HSD "),
    (" 2HB HSD ", " HB1 HSD "),
    (" 3HB HSD ", " HB2 HSD "),
    (" RC HSE ", " C HSE "),
    (" H HSE ", " HN HSE "),
    (" 2HB HSE ", " HB1 HSE "),
    (" 3HB HSE ", " HB2 HSE "),
    (" RC HSP ", " C HSP "),
    (" H HSP ", " HN HSP "),
    (" 2HB HSP ", " HB1 HSP "),
    (" 3HB HSP ", " HB2 HSP "),
]
# (48901) protein label changes PHE
WHATIF_TO_AMBER += [
    (" RC PHE ", " C PHE "),
    (" H PHE ", " HN PHE "),
    (" 2HB PHE ", " HB1 PHE "),
    (" 3HB PHE ", " HB2 PHE "),
]
# (48960) protein label changes ARG
WHATIF_TO_AMBER += [
    (" RC ARG ", " C ARG "),
    (" H ARG ", " HN ARG "),
    (" 2HB ARG ", " HB1 ARG "),
    (" 3HB ARG ", " HB2 ARG "),
    (" 2HG ARG ", " HG1 ARG "),
    (" 3HG ARG ", " HG2 ARG "),
    (" 2HD ARG ", " HD1 ARG "),
    (" 3HD ARG ", " HD2 ARG "),
    (" 1HH1 ARG ", " HH11 ARG "),
    (" 2HH1 ARG ", " HH12 ARG "),
    (" 1HH2 ARG ", " HH21 ARG "),
    (" 2HH2 ARG ", " HH22 ARG "),
]
# (49028) protein label changes TRP
WHATIF_TO_AMBER += [
    (" RC TRP ", " C TRP "),
    (" H TRP ", " HN TRP "),
    (" HB3 TRP ", " HB1 TRP "),
]
# (49052) protein label changes ASN
WHATIF_TO_AMBER += [
    (" RC ASN ", " C ASN "),
    (" H ASN ", " HN ASN "),
    (" 2HB ASN ", " HB1 ASN "),
    (" 3HB ASN ", " HB2 ASN "),
    (" 1HD2 ASN ", " HD21 ASN "),
    (" 2HD2 ASN ", " HD22 ASN "),
]
# (49068) protein label changes PRO
WHATIF_TO_AMBER += [
    (" RC PRO ", " C PRO "),
    (" 2HB PRO ", " HB1 PRO "),
    (" 3HB PRO ", " HB2 PRO "),
    (" 2HG PRO ", " HG1 PRO "),
    (" 3HG PRO ", " HG2 PRO "),
    (" 2HD PRO ", " HD1 PRO "),
    (" 3HD PRO ", " HD2 PRO "),
]
# (49156) protein label changes TYR
WHATIF_TO_AMBER += [
    (" RC TYR ", " C TYR "),
    (" H TYR ", " HN TYR "),
    (" HB3 TYR ", " HB1 TYR "),
]
# (49156) protein label changes ILE
WHATIF_TO_AMBER += [
    (" RC ILE ", " C ILE "),
    (" H ILE ", " HN ILE "),
    (" CD1 ILE ", " CD ILE "),
    (" 2HG1 ILE ", " HG11 ILE "),
    (" 3HG1 ILE ", " HG12 ILE "),
    (" 1HG2 ILE ", " HG21 ILE "),
    (" 2HG2 ILE ", " HG22 ILE "),
    (" 3HG2 ILE ", " HG23 ILE "),
    (" 1HD1 ILE ", " HD1 ILE "),
    (" 2HD1 ILE ", " HD2 ILE "),
    (" 3HD1 ILE ", " HD3 ILE "),
]
# (49361) protein label changes ASP
WHATIF_TO_AMBER += [
    (" RC ASP ", " C ASP "),
    (" H ASP ", " HN ASP "),
    (" 2HB ASP ", " HB1 ASP "),
    (" 3HB ASP ", " HB2 ASP "),
]
# (49394) protein label changes GLN
WHATIF_TO_AMBER += [
    (" RC GLN ", " C GLN "),
    (" H GLN ", " HN GLN "),
    (" 2HB GLN ", " HB1 GLN "),
    (" 3HB GLN ", " HB2 GLN "),
    (" 2HG GLN ", " HG1 GLN "),
    (" 3HG GLN ", " HG2 GLN "),
    (" 1HE2 GLN ", " HE21 GLN "),
    (" 2HE2 GLN ", " HE22 GLN "),
]
# (49431) protein label changes THR
WHATIF_TO_AMBER += [
    (" RC THR ", " C THR "),
    (" H THR ", " HN THR "),
    (" 1HG2 THR ", " HG21 THR "),
    (" 2HG2 THR ", " HG22 THR "),
    (" 3HG2 THR ", " HG23 THR "),
]
# (49445) protein label changes MET
WHATIF_TO_AMBER += [
    (" RC MET ", " C MET "),
    (" H MET ", " HN MET "),
    (" 2HB MET ", " HB1 MET "),
    (" 3HB MET ", " HB2 MET "),
    (" 2HG MET ", " HG1 MET "),
    (" 3HG MET ", " HG2 MET "),
    (" 1HE MET ", " HE1 MET "),
    (" 2HE MET ", " HE2 MET "),
    (" 3HE MET ", " HE3 MET "),
]
# (50445) protein label changes SER
WHATIF_TO_AMBER += [
    (" RC SER ", " C SER "),
    (" H SER ", " HN SER "),
    (" 2HB SER ", " HB1 SER "),
    (" 3HB SER ", " HB2 SER "),
    (" HG SER ", " HG1 SER "),
]
# (56041) protein label changes CYS
WHATIF_TO_AMBER += [
    (" RC CYS ", " C CYS "),
    (" H CYS ", " HN CYS "),
    (" 2HB CYS ", " HB1 CYS "),
    (" 3HB CYS ", " HB2 CYS "),
    (" HG CYS ", " HG1 CYS "),
]
# (52623) errors on O' and O''      
# these seen to occur at the ends of amino acid chains (?)
# so they may be the unlinked O's of the peptide bond.
# can change O' to OT2 etc. but C atoms may need to be changed by hand.
# terminal COO group properties are independent of amino acid, with the
# exception of MET.  but there aren't any C-terminal MET's, so i'm leaving
# that out for now.
WHATIF_TO_AMBER += [
    (" RC CTE ", " C CTE "),
    (" O' ARG ", " OT1 CTE "),
    (" O' GLU ", " OT1 CTE "),
    (" O' VAL ", " OT1 CTE "),
    (" O' GLY ", " OT1 CTE "),
    (" O' ALA ", " OT1 CTE "),
    (" O' TRP ", " OT1 CTE "),
    (" O' THR ", " OT1 CTE "),
    (" O' SER ", " OT1 CTE "),
    (" O' LYS ", " OT1 CTE "),
    (" O'' O ", " OT2 CTE "),
    (" O'' O2 ", " OT2 CTE "),
    (" O'' ARG ", " OT2 CTE "),
    (" O'' ALA ", " OT2 CTE "),
    (" O'' TRP ", " OT2 CTE "),
    (" O'' SER ", " OT2 CTE "),
    (" O'' LYS ", " OT2 CTE "),
    (" O'' GLY ", " OT2 CTE "),
]
# N-terminal end: replace N with N NTE (prob. have to do by hand)
# replace H1 with HT1 NTE with exceptions PRO, SER.
WHATIF_TO_AMBER += [
    (" H1 GLY ", " HT1 NTE "),
    (" H2 GLY ", " HT2 NTE "),
    (" H3 GLY ", " HT3 NTE "),
    (" H1 ASP ", " HT1 NTE "),
    (" H2 ASP ", " HT2 NTE "),
    (" H3 ASP ", " HT3 NTE "),
    (" H1 MET ", " HT1 NTE "),
    (" H2 MET ", " HT2 NTE "),
    (" H3 MET ", " HT3 NTE "),
    (" H1 ALA ", " HT1 NTE "),
    (" H2 ALA ", " HT2 NTE "),
    (" H3 ALA ", " HT3 NTE "),
    (" H1 GLU ", " HT1 NTE "),
    (" H2 GLU ", " HT2 NTE "),
    (" H3 GLU ", " HT3 NTE "),
    (" H1 LYS ", " HT1 NTE "),
    (" H2 LYS ", " HT2 NTE "),
    (" H3 LYS ", " HT3 NTE "),
    (" H1 ARG ", " HT1 NTE "),
    (" H2 ARG ", " HT2 NTE "),
    (" H3 ARG ", " HT3 NTE "),
]
# weird glutamate protonation 2HE (55204)
# just leave it for now and add fake param to amber.
# N-terminal PRO (73791)
# reconfig to PRON (by hand)
WHATIF_TO_AMBER += [
    (" H2 PRO ", " HN1 PRON "),
    (" H3 PRO ", " HN2 PRON "),
]
WHATIF_TO_AMBER = dict(WHATIF_TO_AMBER)
AMBER_TO_WHATIF = dict(
    (a, w) for w, a in WHATIF_TO_AMBER.items()
)
WHATIF_TO_CHARMM = {}
for whatif, amber in WHATIF_TO_AMBER.items():
    try:
        charmm = AMBER_TO_CHARMM[amber]
        WHATIF_TO_CHARMM[whatif] = charmm
        print(whatif, charmm)
    except KeyError:
        print("Skipping %s" % whatif)