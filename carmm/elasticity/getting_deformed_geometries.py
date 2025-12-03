def generate_deformed_strutures(atoms_object, norm_strains = [0.01, 0.01], shear_strains = [0.01, 0.01], write_strain=True):

    import pickle
    import numpy as np
    from pymatgen.analysis.elasticity import Strain, Deformation
    from pymatgen.io.ase import AseAtomsAdaptor

    eq_structure = atoms_object
    structure = AseAtomsAdaptor.get_structure(eq_structure)

    deformations: list[Deformation] = []

    strain_list = []
    for ind in [(0, 0), (1, 1), (2, 2)]:
        for amount in norm_strains:
            strain = Strain.from_index_amount(ind, amount)
            strain_list.append(strain)
            deformations.append(strain.get_deformation_matrix())
    for ind in [(0, 1), (0, 2), (1, 2)]:
        for amount in shear_strains:
            strain = Strain.from_index_amount(ind, amount)
            strain_list.append(strain)
            deformations.append(strain.get_deformation_matrix())

    if write_strain:
        strain_tensor = np.array(strain_list)
        with open('strain_tensor.pkl', 'wb') as fp:
            pickle.dump(strain_tensor, fp)

    return eq_structure, structure, deformations

def create_files_and_directories(eq_structure, structure, deformations, copy_input_and_submission=False):
    import os
    import shutil
    from pymatgen.io.ase import AseAtomsAdaptor
    deformed_struct = [defo.apply_to_structure(structure) for defo in deformations]

    for i, def_struc in enumerate(deformed_struct):
        dir_no = i + 1
        eq_structure_copy = eq_structure.copy()
        directory = f'defor_{dir_no}'
        parent_dir = os.getcwd()
        path_final = os.path.join(parent_dir, directory)
        if os.path.exists(path_final):
            shutil.rmtree(path_final)
        os.mkdir(path_final)
        atoms = AseAtomsAdaptor.get_atoms(def_struc)
        eq_structure_copy.set_cell(atoms.get_cell(), scale_atoms=True)
        eq_structure_copy.write(path_final + '/geometry.in')
        if copy_input_and_submission:
            shutil.copy(parent_dir + '/input.py', path_final + '/input.py')
            shutil.copy(parent_dir + '/submission.script', path_final + '/submission.script')





