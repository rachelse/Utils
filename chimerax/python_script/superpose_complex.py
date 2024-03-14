'''
File: superpose_complex.py
Project: chimerax
File Created: 7th Mar 2024
Author: Rachel Seongeun Kim (seamustard52@gmail.com)
-----
Copyright: Rachel Seongeun Kim
'''

"""
How to use
----------
open superpose_complex.py
# open query
open #1
#open target
open #2
superpose_complex #2 0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012
"""

def rotate(session, atoms, ut):
    if len(ut) != 12:
        raise ValueError("ut must be a 12-element list")
    xyz_rotated = atoms.coords
    xyz_rotated[:,0] = (ut[0]*atoms.coords[:,0] + ut[1]*atoms.coords[:,1] + ut[2]*atoms.coords[:,2]) + ut[9]
    xyz_rotated[:,1] = (ut[3]*atoms.coords[:,0] + ut[4]*atoms.coords[:,1] + ut[5]*atoms.coords[:,2]) + ut[10]
    xyz_rotated[:,2] = (ut[6]*atoms.coords[:,0] + ut[7]*atoms.coords[:,1] + ut[8]*atoms.coords[:,2]) + ut[11]
    atoms.coords = xyz_rotated

def register_command(logger):
    from chimerax.core.commands import CmdDesc, register, FloatsArg
    from chimerax.atomic import AtomsArg
    desc = CmdDesc(
        required = [('atoms', AtomsArg),
                    ('ut', FloatsArg)
                    # ('scorecomplex_report', OpenFileNameArg),
                    ],
        synopsis = 'Rotate the structure with the given unitary matrix and translation vector.'
    )
    register('superpose_complex', desc, rotate, logger=logger)

register_command(session.logger)