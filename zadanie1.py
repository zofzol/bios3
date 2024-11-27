import os
import math
from Bio import PDB
import matplotlib.pyplot as plt
import numpy as np

# Tworzenie obiektu do pobierania struktur
pdb_list = PDB.PDBList()
pdb_id = '4YWO'

# Pobieranie pliku PDB (jeśli jeszcze nie został pobrany)
if not os.path.exists(f"{pdb_id}.cif"):
    pdb_list.retrieve_pdb_file(pdb_id, file_format='mmCif', pdir='.')

# Wczytywanie struktury za pomocą MMCIFParser
mmcif_parser = PDB.MMCIFParser(QUIET=True)
structure = mmcif_parser.get_structure(pdb_id, f"{pdb_id}.cif")

# Pobieranie współrzędnych atomów węgla alfa
coordinates = []
for model in structure:
    for chain in model:
        for residue in chain:
            # Upewnienie się, że reszta zawiera atomy
            if 'CA' in residue:
                ca_atom = residue['CA']
                coordinates.append(list(ca_atom.get_coord()))

# Obliczanie macierzy kontaktów
threshold = 8.0  # Próg odległości (w Angstromach)
matrix = []
for c in coordinates:
    tmp = []
    for d in coordinates:
        euk_dis = math.sqrt((c[0] - d[0])**2 + (c[1] - d[1])**2 + (c[2] - d[2])**2)
        tmp.append(1 if euk_dis < threshold else 0)
    matrix.append(tmp)

# Konwersja macierzy na numpy array dla łatwej manipulacji
contact_map = np.array(matrix)

# Wizualizacja macierzy kontaktów
plt.figure(figsize=(10, 8))
plt.title(f"Contact Map for {pdb_id}", fontsize=16)
plt.imshow(contact_map, cmap='Greys', interpolation='none')
plt.colorbar(label='Contact (1=Yes, 0=No)')
plt.xlabel('Residue Index')
plt.ylabel('Residue Index')
plt.savefig("contact_map_for_4YWO.png", dpi=300, bbox_inches='tight')  # dpi zwiększa jakość, bbox_inches zmniejsza marginesy
plt.show()

