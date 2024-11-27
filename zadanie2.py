from Bio import PDB
import numpy as np
import matplotlib.pyplot as plt

# Funkcja obliczająca kąty phi i psi dla danego obiektu struktury białka
def calculate_phi_psi(structure):
    phi_psi_angles = []  # Lista do przechowywania kątów phi i psi
    for model in structure:  # Iteracja po modelach w strukturze (dla plików z wieloma modelami)
        for chain in model:  # Iteracja po łańcuchach w modelu
            residues = list(chain)  # Pobranie wszystkich reszt w łańcuchu
            for i in range(1, len(residues) - 1):  # Iteracja po resztach (pomijając pierwszą i ostatnią)
                
                # Obliczanie kąta phi: sprawdzamy, czy wymagane atomy są obecne
                if residues[i - 1].has_id("C") and residues[i].has_id("N") and residues[i].has_id("CA") and residues[i].has_id("C"):
                    try:
                        # Obliczenie kąta torsyjnego phi przy użyciu czterech atomów
                        phi = PDB.calc_dihedral(
                            residues[i - 1]["C"].get_vector(),  # C poprzedniego resztka
                            residues[i]["N"].get_vector(),      # N bieżącego resztka
                            residues[i]["CA"].get_vector(),     # Cα bieżącego resztka
                            residues[i]["C"].get_vector(),      # C bieżącego resztka
                        )
                    except Exception:
                        phi = None  # Jeśli obliczenie się nie powiedzie, ustawiamy phi na None
                else:
                    phi = None

                # Obliczanie kąta psi: sprawdzamy, czy wymagane atomy są obecne
                if residues[i].has_id("N") and residues[i].has_id("CA") and residues[i].has_id("C") and residues[i + 1].has_id("N"):
                    try:
                        # Obliczenie kąta torsyjnego psi przy użyciu czterech atomów
                        psi = PDB.calc_dihedral(
                            residues[i]["N"].get_vector(),      # N bieżącego resztka
                            residues[i]["CA"].get_vector(),     # Cα bieżącego resztka
                            residues[i]["C"].get_vector(),      # C bieżącego resztka
                            residues[i + 1]["N"].get_vector(),  # N następnego resztka
                        )
                    except Exception:
                        psi = None  # Jeśli obliczenie się nie powiedzie, ustawiamy psi na None
                else:
                    psi = None

                # Dodanie pary (phi, psi) do listy, jeśli obydwa kąty są poprawnie obliczone
                if phi is not None and psi is not None:
                    phi_psi_angles.append((np.degrees(phi), np.degrees(psi)))  # Konwersja z radianów na stopnie

    return phi_psi_angles  # Zwracamy listę par (phi, psi)

# Główna część programu
if __name__ == "__main__":
    pdb_file = "4YWO.pdb"  # Nazwa pliku PDB, który zawiera strukturę białka
    parser = PDB.PDBParser(QUIET=True)  # Ustawiamy QUIET=True, aby uniknąć komunikatów ostrzegawczych
    structure = parser.get_structure("4YWO", pdb_file)  # Parsowanie pliku PDB i załadowanie struktury

    # Obliczanie kątów phi i psi dla struktury
    phi_psi_angles = calculate_phi_psi(structure)

    # Tworzenie wykresu Ramachandrana
    if phi_psi_angles:  # Jeśli mamy jakieś wartości phi i psi
        phi, psi = zip(*phi_psi_angles)  # Rozdzielenie listy par (phi, psi) na dwie listy: phi i psi
        plt.figure(figsize=(8, 8))  # Ustawienie rozmiaru wykresu
        plt.scatter(phi, psi, alpha=0.7, s=10, c="blue")  # Rysowanie punktów na wykresie
        plt.xlim(-180, 180)  # Ustawienie zakresu osi X (phi)
        plt.ylim(-180, 180)  # Ustawienie zakresu osi Y (psi)
        plt.xlabel("Phi (°)")  # Etykieta osi X
        plt.ylabel("Psi (°)")  # Etykieta osi Y
        plt.title("Wykres Ramachandrana dla 4YWO")  # Tytuł wykresu
        plt.grid(True)  # Włączenie siatki na wykresie
        plt.savefig("ramachandran_4YWO.png")  # Zapisanie wykresu do pliku PNG
        plt.show()  # Wyświetlenie wykresu
    else:
        print("Nie udało się obliczyć kątów phi i psi.")
