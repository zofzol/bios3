from Bio.PDB.PDBList import PDBList
from Bio import PDB
import Bio.PDB.Residue
from Bio.PDB import PDBIO, Atom, Residue, Chain, Model, Structure
import Bio.PDB.PDBIO
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

pdb_parser = PDB.PDBParser(QUIET=True)
cg_structure = pdb_parser.get_structure('430D_cg', r'C:\Users\zosia\OneDrive\Pulpit\sztudia\laby\kody\bs3\430d_cg.pdb')

adenina = pdb_parser.get_structure('ade', r'C:\Users\zosia\OneDrive\Pulpit\sztudia\laby\kody\bs3\templates\a.pdb')
guanina = pdb_parser.get_structure('gua', r'C:\Users\zosia\OneDrive\Pulpit\sztudia\laby\kody\bs3\templates\g.pdb')
uracyl = pdb_parser.get_structure('ura', r'C:\Users\zosia\OneDrive\Pulpit\sztudia\laby\kody\bs3\templates\u.pdb')
cytozyna = pdb_parser.get_structure('cyt', r'C:\Users\zosia\OneDrive\Pulpit\sztudia\laby\kody\bs3\templates\c.pdb')
backbone1 = pdb_parser.get_structure('rib', r'C:\Users\zosia\OneDrive\Pulpit\sztudia\laby\kody\bs3\templates\r1.pdb')
backbone9 = pdb_parser.get_structure('rib', r'C:\Users\zosia\OneDrive\Pulpit\sztudia\laby\kody\bs3\templates\r9.pdb')

templates = {
    "A": adenina,
    "G": guanina,
    "U": uracyl,
    "C": cytozyna,
    "B1": backbone1,
    "B9": backbone9
}

restored_structure = Bio.PDB.Structure.Structure('nowa')
model = Bio.PDB.Model.Model(0)
restored_structure.add(model)
chain = Bio.PDB.Chain.Chain("A")
model.add(chain)

lista_reszt = [adenina, guanina, uracyl, cytozyna]

#iteracja po atomach wczytanej struktury
def save_atoms_to_list(cg_structure):
    list_coords = []
    side_chain_atoms = []

    for model in cg_structure:
        for chain in model:
            for res in chain:
                list_atoms = []
                side_chain = []
                for atom in res.get_atoms():
                    if atom.get_name() in ['P', 'C4\'']: #tutaj mam atomy szkieletu (backbone) dla danej reszty
                        list_atoms.append(atom)
                        print("Znaleziono rybke w ZOO")
                    else:
                        side_chain.append(atom)
                        if res.resname in ['A', 'G']:
                            if atom.get_name()=="N9":
                                list_atoms.append(atom)
                        if res.resname in ['C', 'U']:
                            if atom.get_name()=="N1":
                                list_atoms.append(atom)     
                        print("nie ma rybki :c")
                list_coords.append(list_atoms) # Listę atomów zapisuję do kolejnej listy. Dzięki temu łatwo mogę wybrać atomy z konkretnej reszty, która mnie interesuje.
                side_chain_atoms.append(side_chain)
    # print(list_coords)
    # print(side_chain_atoms)
    return list_coords, side_chain_atoms


def save_residue_atoms_to_list(templates): 
    list_coords = []

    for structure in templates:
        for model in structure:
            for chain in model:
                for res in chain:
                    list_atoms = []
                    for atom in res.get_atoms():
                        list_atoms.append(atom)
                    list_coords.append(list_atoms)
    print(list_coords)
    print(f"Number of atoms added: {len(list_coords)}")

    return list_coords

#ekstrakcja atomów do macierzy superimmpotencji
def save_cg_residue_atoms_to_list(templates):
    super = []
    for structure in templates:
        for model in structure:
            for chain in model:
                for resi in chain:
                    if resi.resname == "A" or resi.resname == "G":  # Purine (Adenine and Guanine)
                        atoms_to_S = []
                        for atom in resi.get_atoms():
                            if atom.get_name() in ["N9", "C2", "C6"]:  # Specific atoms for Purines
                                atoms_to_S.append(atom)
                        super.append(atoms_to_S)
                    elif resi.resname == "C" or resi.resname == "U":  # Pyrimidine (Cytosine and Uracil)
                        atoms_to_S = []
                        for atom in resi.get_atoms():
                            if atom.get_name() in ["N1", "C2", "C4"]:  # Specific atoms for Pyrimidines
                                atoms_to_S.append(atom)
                        super.append(atoms_to_S)
                    
    print(f"Number of atoms added: {len(super)}")
    return super



St_430d_backbone_part, St_430d_residue_part = save_atoms_to_list(cg_structure) # podział cd_structure na backbone i nie backbone

cg_backbone_1_list, smieciowe = save_atoms_to_list(backbone1) #wyodrębnienie części cg z backbone +N1 
cg_backbone_9_list, smieciowe = save_atoms_to_list(backbone9)#wyodrębnienie części cg z backbone +N9
list_of_atoms_coords_template = save_residue_atoms_to_list(templates.values())
cg_res_template_list = save_cg_residue_atoms_to_list(lista_reszt) #wyodrębnienie części cg z reszt 


#supepepermint mimoser maxxing
superimposer = PDB.Superimposer()


nr_reszty = 0
nr_atomu = 0
# Iteracja po fragmentach szkieletu
for i in range(len(St_430d_backbone_part)):
   
    # Tworzenie nowej reszty
    reszta = St_430d_backbone_part[i][0].get_parent()  # Pobranie reszty
    nr_reszty += 1
    residue = Bio.PDB.Residue.Residue((" ", nr_reszty, " "), reszta.resname, " ")
    chain.add(residue)


    # Wybór template dla różnych reszt
    match reszta.resname:
        case "A":
            pobrany_template = list_of_atoms_coords_template[0]
            do_przesuniecia = cg_res_template_list[0]
        case "G":
            pobrany_template = list_of_atoms_coords_template[1]
            do_przesuniecia = cg_res_template_list[1]
        case "U":
            pobrany_template = list_of_atoms_coords_template[2]
            do_przesuniecia = cg_res_template_list[2]
        case "C":
            pobrany_template = list_of_atoms_coords_template[3]
            do_przesuniecia = cg_res_template_list[3]

    #szkielet = list_of_atoms_coords_template[4]  # Wybór szkieletu 4 --> backbone 1 ; 5 --> bacckbone 9
    jaka_to_reszta = St_430d_backbone_part[i][0].get_parent()

    if jaka_to_reszta.resname == "U" or jaka_to_reszta.resname == "C":
        szkielet = list_of_atoms_coords_template[4]
        superimposer.set_atoms(St_430d_backbone_part[i], cg_backbone_1_list[0])  # Superimpozycja

    elif jaka_to_reszta.resname == "A" or jaka_to_reszta.resname=="G":
        szkielet = list_of_atoms_coords_template[5]
        superimposer.set_atoms(St_430d_backbone_part[i], cg_backbone_9_list[0])  # Superimpozycja
    
    superimposer.apply(szkielet)

    superimposer.set_atoms(fixed=St_430d_residue_part[i], moving=do_przesuniecia)
    superimposer.apply(pobrany_template)

     # Dodawanie atomów do reszty
    for skeleton_atom in szkielet:
        if skeleton_atom.get_name() == "N1" or skeleton_atom.get_name() == "N9" :
            continue #skip
        nr_atomu += 1
        atom_name = skeleton_atom.get_name()
        atom_coordinates = skeleton_atom.get_coord()
        element = atom_name[0]
        
        # Tworzenie nowego atomu na podstawie danych z szkieletu
        new_atom = Bio.PDB.Atom.Atom(
            name=atom_name,                    # Nazwa atomu, np. "P", "C4'"
            coord=atom_coordinates,            # Współrzędne atomu (lista [x, y, z])
            bfactor=0.0,                       # Współczynnik temperaturowy (domyślny)
            occupancy=1.0,                     # Zajętość atomu (domyślnie 1.0)
            altloc=" ",                        # Alternatywne lokacje (domyślnie brak)
            fullname=" " + atom_name,          # Pełna nazwa atomu z wyrównaniem
            serial_number=nr_atomu,            # Unikalny numer atomu
            element=element                    # Symbol pierwiastka
        )
        residue.add(new_atom)


    # Dodanie atomów z szablonu (template) do reszty
    for template_atom in pobrany_template:
        nr_atomu += 1
        atom_name = template_atom.get_name()
        atom_coordinates = template_atom.get_coord()  
        element = atom_name[0]  

        # Tworzenie nowego atomu na podstawie danych z szablonu
        new_atom = Bio.PDB.Atom.Atom(
            name=atom_name,                    # Nazwa atomu
            coord=atom_coordinates,            # Współrzędne atomu (lista [x, y, z])
            bfactor=0.0,                       # Współczynnik temperaturowy (domyślny)
            occupancy=1.0,                     # Zajętość atomu (domyślnie 1.0)
            altloc=" ",                        # Alternatywne lokacje (domyślnie brak)
            fullname=" " + atom_name,          # Pełna nazwa atomu z wyrównaniem
            serial_number=nr_atomu,            # Unikalny numer atomu
            element=element                    # Symbol pierwiastka
        )
        residue.add(new_atom)         # Dodanie nowego atomu do aktualnej reszty


# Tworzenie obiektu PDBIO
io = Bio.PDB.PDBIO()
io.set_structure(restored_structure)  # Ustawienie struktury do zapisania

# Zapisanie do pliku w formacie PDB
with open("430d_res.pdb", "w", encoding="utf-8") as f:
    io.save(f)
