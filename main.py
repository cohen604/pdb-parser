from typing import Optional

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


def main():
    file_path = "Bfr2_align122resample.pdb"

    chains = []
    with open(file_path, 'r', encoding="UTF-8") as fhand:
        chain = []
        for line in fhand:
            if not line.startswith("ATOM"):
                if line.startswith("TER"):
                    chains.append(chain)
                    chain = []
                    continue
                continue
            atom = AtomRecord(line)
            chain.append(atom)

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

    print(amino_chain_string_list)


if __name__ == '__main__':
    main()
