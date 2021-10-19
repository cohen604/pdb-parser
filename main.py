from typing import Optional, Union, List
import argparse

THREE_SYMBOLS_TO_ONE = {
    "Ala": "A",
    "Asx": "B",
    "Cys": "C",
    "Asp": "D",
    "Glu": "E",
    "Phe": "F",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Lys": "K",
    "Leu": "L",
    "Met": "M",
    "Asn": "N",
    "Pro": "P",
    "Gln": "Q",
    "Arg": "R",
    "Ser": "S",
    "Thr": "T",
    "Val": "V",
    "Trp": "W",
    "Tyr": "Y"
}


class AtomRecord:
    # Took it from https://www.wwpdb.org/documentation/file-format-content/format23/sect9.html#ATOM
    ATOM_RECORD_FIELDS = {
        "Record": slice(0, 6),
        "Serial": slice(6, 11),
        "Atom Name": slice(12, 16),
        "Residue Name": slice(17, 20),
        "ChainID": slice(21, 22),
        "Residue Sequence": slice(22, 26),
        "x": slice(30, 38),
        "y": slice(38, 46),
        "z": slice(46, 54),
        "Occupancy": slice(54, 60),
        "Temperature Factor": slice(60, 66)
    }

    record: str
    serial: int
    atom_name: str
    residue_name: str
    chain_ID: chr
    residue_sequence: int
    x: float
    y: float
    z: float
    occupancy: float
    temperature_factor: float

    def __init__(self, data: str):
        self._parse_data(data.strip())
        self.residue_name_one_letter = THREE_SYMBOLS_TO_ONE[self.residue_name.capitalize()]

    def _parse_data(self, data: str):
        self._parse_record(data[self.ATOM_RECORD_FIELDS["Record"]])
        self._parse_serial(data[self.ATOM_RECORD_FIELDS["Serial"]])
        self._parse_atom_name(data[self.ATOM_RECORD_FIELDS["Atom Name"]])
        self._parse_residue_name(data[self.ATOM_RECORD_FIELDS["Residue Name"]])
        self._parse_chain_id(data[self.ATOM_RECORD_FIELDS["ChainID"]])
        self._parse_residue_sequence(data[self.ATOM_RECORD_FIELDS["Residue Sequence"]])
        self._parse_x(data[self.ATOM_RECORD_FIELDS["x"]])
        self._parse_y(data[self.ATOM_RECORD_FIELDS["y"]])
        self._parse_z(data[self.ATOM_RECORD_FIELDS["z"]])
        self._parse_occupancy(data[self.ATOM_RECORD_FIELDS["Occupancy"]])
        self._parse_temperature_factor(data[self.ATOM_RECORD_FIELDS["Temperature Factor"]])

    def _parse_record(self, record: str):
        self.record = record.strip()

    def _parse_serial(self, serial: str):
        self.serial = int(serial.strip())

    def _parse_atom_name(self, atom_name: str):
        self.atom_name = atom_name.strip()

    def _parse_residue_name(self, residue_name: str):
        self.residue_name = residue_name.strip()

    def _parse_chain_id(self, chain_id: str):
        self.chain_ID = chain_id.strip()

    def _parse_residue_sequence(self, residue_sequence: str):
        self.residue_sequence = int(residue_sequence.strip())

    def _parse_x(self, x: str):
        self.x = float(x.strip())

    def _parse_y(self, y: str):
        self.y = float(y.strip())

    def _parse_z(self, z: str):
        self.z = float(z.strip())

    def _parse_occupancy(self, occupancy: str):
        self.occupancy = float(occupancy.strip())

    def _parse_temperature_factor(self, temperature_factor: str):
        self.temperature_factor = float(temperature_factor.strip())


class Chain:
    atoms: List[AtomRecord]
    chain_id: Optional[str]

    def __init__(self, atoms: List[AtomRecord] = []):
        self.atoms = atoms
        self.chain_id = None
        self._set_chain_id()

    def change_occupancy(self, occupancy: float):
        for atom in self.atoms:
            atom.occupancy = occupancy

    def append_atom(self, atom: AtomRecord):
        self.atoms.append(atom)
        if self.chain_id is None:
            self._set_chain_id()

    def _set_chain_id(self):
        if len(self.atoms) > 0 and self.chain_id is None:
            self.chain_id = self.atoms[0].chain_ID


class PDB:
    chains: List[Chain]

    def __init__(self, chains: List[Chain] = []):
        self.chains = chains

    def append_chain(self, chain: Chain):
        self.chains.append(chain)

    def change_occupancy(self, occupancy: float, chain_id: str = "ALL"):
        for chain in self.chains:
            if chain_id == chain.chain_id or chain_id == "ALL":
                chain.change_occupancy(occupancy)

    def write_to_pdb_file(self):
        pass


def parse_pdb_from_file(pdb_path: str) -> PDB:
    pdb = PDB()
    with open(pdb_path, 'r', encoding="UTF-8") as fhand:
        chain = Chain([].copy())
        for line in fhand:
            if not line.startswith("ATOM"):
                if line.startswith("TER"):
                    pdb.append_chain(chain)
                    chain = Chain([].copy())
                    continue
                continue
            atom = AtomRecord(line)
            chain.append_atom(atom)

    return pdb


def extract_amino_chain(chains: List[List[AtomRecord]]) -> str:
    unique_amino_acids_chains = []
    for chain in chains:
        unique_amino_atom = set()
        for atom in chain:
            unique_amino_atom.add((atom.residue_sequence, atom.residue_name_one_letter))
        sorted_unique_amino_atom = sorted(unique_amino_atom, key=lambda x: x[0])
        unique_amino_acids_chains.append(sorted_unique_amino_atom)

    amino_chain_string_list = set()
    for amino_chain in unique_amino_acids_chains:
        amino_chain_string = ''.join([value for key, value in amino_chain])
        amino_chain_string_list.add(amino_chain_string)

    return str(amino_chain_string_list)


def change_occupancy(pdb_file, output_file: str, occupancy: List[str], chain_id: Union[List[str], str, None] = "ALL"):
    ATOM_RECORD_FIELDS = {
        "ChainID": slice(21, 22),
        "Occupancy": slice(54, 60),
    }
    changes = {}
    if chain_id is None:
        chain_id = "ALL"
        changes[chain_id] = occupancy[0]
    else:
        for key, value in zip(chain_id, occupancy):
            changes[key] = value

    with open(output_file, "w") as results:
        for line in pdb_file:
            if line.startswith("ATOM"):
                atom_chain_id = line[ATOM_RECORD_FIELDS["ChainID"]].strip()
                num_white_spaces = (ATOM_RECORD_FIELDS["Occupancy"].stop - ATOM_RECORD_FIELDS["Occupancy"].start)
                if atom_chain_id in changes.keys():
                    occ = f'{changes[atom_chain_id]: >{num_white_spaces}}'
                elif chain_id == "ALL":
                    occ = f'{changes[chain_id]: >{num_white_spaces}}'
                line = line[:ATOM_RECORD_FIELDS["Occupancy"].start] + occ + line[ATOM_RECORD_FIELDS["Occupancy"].stop:]
            results.write(line)


def main():
    parser = argparse.ArgumentParser(description='PDB Tools')
    sub_parsers = parser.add_subparsers(help="help", dest='command')

    change_occupancy_parser = sub_parsers.add_parser('change-occupancy', help="help")
    change_occupancy_parser.add_argument('-i', '--input', dest='input_file_name', help="File to parse",required=True)
    change_occupancy_parser.add_argument('-o', '--output', dest='output_file_name', help="File output", required=True)
    change_occupancy_parser.add_argument('--occupancy', dest='occupancy', type=str, nargs='*', help="Occupancy number", required=True)
    change_occupancy_parser.add_argument('--chain-id', dest='chain_id', type=str, nargs='*', help="Chain ID to change occupancy to")

    args = parser.parse_args()

    if args.command == "change-occupancy":
        with open(args.input_file_name, "r") as fh:
            change_occupancy(pdb_file=fh, output_file=args.output_file_name, occupancy=args.occupancy, chain_id=args.chain_id)


if __name__ == '__main__':
    main()
