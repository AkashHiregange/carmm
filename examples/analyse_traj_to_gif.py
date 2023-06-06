'''
Example script showing how to convert an example trajectory file (.traj)
into a gif visualised in povray.
'''

def test_traj_to_gif():
    from carmm.analyse.traj_to_gif import traj_to_gif, gifmaker, atom_subs
    from ase.io import read

    file = 'data/NH3-H3O_traj/nh3-h3o.traj'

    # Call the overall function like this, the rest of the script just tests the functions individually
    # traj_to_gif(file, frames_per_second=10, pause_time=1, atom_subs=[['N', 'C']], keep_temp_files=True)

    atoms = read(file)
    frame_atoms = atom_subs(atoms, [['N', 'C']], 41, 'nh3-h3o', 'traj')
    assert len(frame_atoms) == 41
    assert frame_atoms.symbols[4] == 'C'

    traj_atoms = read('nh3-h3o-povray.traj')
    assert len(traj_atoms) == 41
    assert traj_atoms.symbols[4] == 'C'

    filenames, delay = gifmaker(41, 'nh3-h3o', 'traj', frames_per_second=10, pause_time=1, keep_png_files=False)
    assert delay == '1x10'
    assert len(filenames) == 59
