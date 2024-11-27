import argparse
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import Select
import os

# Ładowanie struktury PDB za pomocą PDBList
pdb_list = PDB.PDBList()

# Ustawienie parsera, który umożliwia czytanie plików PDB
parser = PDB.PDBParser(QUIET=True)

# Wczytanie struktury z lokalnego pliku PDB
# Używamy ścieżki bezpośrednio do pliku na komputerze, możesz to zmienić na odpowiednią lokalizację.
structure = parser.get_structure('430d', r'C:\Users\zosia\OneDrive\Pulpit\sztudia\laby\kody\bs3\430d.pdb')

# Klasa do selekcji atomów, które mają być zapisane w reprezentacji gruboziarnistej
class Gzik(Select):
    def __init__(self):
        self.first_c4_removed = False  # Flaga, która śledzi, czy pierwszy atom C4' został już usunięty
        self.first_n9_removed = False  # Flaga, która śledzi, czy pierwszy atom N9 został już usunięty
        self.first_c6_removed = False  # Flaga, która śledzi, czy pierwszy atom C6 został już usunięty
        self.first_c2_removed = False  # Flaga, która śledzi, czy pierwszy atom C2 został już usunięty

    def accept_atom(self, atom):
        """
        Funkcja sprawdzająca, które atomy mają być uwzględnione w reprezentacji gruboziarnistej.
        Jeśli napotkamy pierwszy atom C4', N9, C6 lub C2, to je pomijamy.
        """
        # Definicja atomów, które będą zachowane
        puryny = {"N9", "C2", "C6"}  # Atomy purynowe
        pirymidyny = {"N1", "C2", "C4"}  # Atomy pirymidynowe
        szkielet = {"P", "C4'"}  # Atomy szkieletowe

        # Pobieramy nazwę atomu i nazwę reszty (np. A, G, C, U)
        nazwa_atomu = atom.get_name()
        reszta = atom.get_parent().get_resname()

        # Sprawdzamy, czy napotkaliśmy pierwszy raz atom C4', N9, C6 lub C2
        if nazwa_atomu == "C4'" and not self.first_c4_removed:
            self.first_c4_removed = True  # Zmieniamy flagę, aby kolejny atom C4' był zachowany
            return False  # Usuwamy ten pierwszy atom C4'
        elif nazwa_atomu == "N9" and not self.first_n9_removed:
            self.first_n9_removed = True  # Zmieniamy flagę, aby kolejny atom N9 był zachowany
            return False  # Usuwamy ten pierwszy atom N9
        elif nazwa_atomu == "C6" and not self.first_c6_removed:
            self.first_c6_removed = True  # Zmieniamy flagę, aby kolejny atom C6 był zachowany
            return False  # Usuwamy ten pierwszy atom C6
        elif nazwa_atomu == "C2" and not self.first_c2_removed:
            self.first_c2_removed = True  # Zmieniamy flagę, aby kolejny atom C2 był zachowany
            return False  # Usuwamy ten pierwszy atom C2

        # Sprawdzamy, czy atom należy do puryn lub pirymidyn, czy jest częścią szkieletu
        if reszta == "A" or reszta == "G":  # Puryny (Adenina, Guanina)
            return nazwa_atomu in puryny or nazwa_atomu in szkielet
        elif reszta == "C" or reszta == "U":  # Pirymidyny (Cytozyna, Uracyl)
            return nazwa_atomu in pirymidyny or nazwa_atomu in szkielet
        return False  # Jeśli atom nie należy do żadnej z grup, nie jest uwzględniany

# Nazwa pliku wyjściowego, do którego zapisujemy gruboziarnistą reprezentację
output_file = "430d_cg.pdb"

# Obiekt do zapisywania struktury
io = PDBIO()

# Ustawiamy strukturę, którą chcemy zapisać
io.set_structure(structure)

# Zapisujemy strukturę do pliku, filtrując atomy za pomocą klasy Gzik
io.save(output_file, select=Gzik())

print(f"Reprezentacja gruboziarnista zapisana do pliku: {output_file}")
