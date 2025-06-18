'''
File: visualize_foldseek_cluster.py
Project: python_script
File Created: 13th Oct 2024
Author: Rachel Seongeun Kim (seamustard52@gmail.com)
-----
Copyright: Rachel Seongeun Kim
'''

"""
How to use
----------
open visualize_foldseek_cluster.py

# 1. show the cluster
show_cluster </path/to/cluster_result_file> </path/to/cluster_report_file> pdb_path </path/to/pdb>

# 2. apply rotate and translate with superposition matrix
rotate_translate <model identifier> <u,t matrix>
e.g. rotate_translate #2 -0.542,0.828,-0.138,-0.822,-0.557,-0.115,-0.172,0.050,0.983,128.937,343.236,-8.978

----------------------
Additional Information
----------------------
# 1. use the cluster as a unit e.g. show specific cluster (clu_7soy_1)
hide cartoon
show clu_7soy_1 cartoon
"""

import sys
def _pdb_handler(query_target):
    if len(query_target) ==1:
        if ":" in query_target[0]:
            query, target = query_target[0].split(":")
        else:
            print("Please provide the query and target")
    elif len(query_target) == 2:
        query = query_target[0]
        target = query_target[1]
    else:
        print("Please provide the query and target")
        return
    if ".pdb" not in query:
        query += ".pdb"
    if ".pdb" not in target:
        target += ".pdb"

    return query, target

def rotate_translate(session, atoms, ut):
    if len(ut) != 12:
        raise ValueError("ut must be a 12-element list")
    xyz_rotated = atoms.coords
    xyz_rotated[:,0] = (ut[0]*atoms.coords[:,0] + ut[1]*atoms.coords[:,1] + ut[2]*atoms.coords[:,2]) + ut[9]
    xyz_rotated[:,1] = (ut[3]*atoms.coords[:,0] + ut[4]*atoms.coords[:,1] + ut[5]*atoms.coords[:,2]) + ut[10]
    xyz_rotated[:,2] = (ut[6]*atoms.coords[:,0] + ut[7]*atoms.coords[:,1] + ut[8]*atoms.coords[:,2]) + ut[11]
    atoms.coords = xyz_rotated

def parse_report(cluster_info, report): # filtermultimer_report : query, target, qCov, tCov, qTm, tTM, ilddt, u, t
    cluster = {}
    ut = {}
    entry_id = {}
    last_id = 1
    REP_IDX, MEM_IDX=1,2
    # QC_IDX, TC_IDX=3,4 # matched REP chains (query), matched MEM chains (target)
    # QCOV_IDX, TCOV_IDX=5,6
    # QTM_IDX, TTM_IDX=7,8
    # QCTM_IDX, TCTM_IDX=9,10
    # ILDDT_IDX=11
    U_IDX, T_IDX= 12, 13

    with open(report, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line_splitted = line.split("\t")
            rep = line_splitted[REP_IDX]
            mem = line_splitted[MEM_IDX]
            if rep not in ut:
                ut[rep] = {}
            ut[rep][mem] = ",".join([line_splitted[U_IDX], line_splitted[T_IDX]])
            # else :
            #     if line_splitted[1] not in ut[line_splitted[0]]:
            #         ut[line_splitted[0]][line_splitted[1]] = ",".join([line_splitted[7], line_splitted[8]])

    with open(cluster_info, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line_splitted = line.strip().split("\t")
            rep, mem = line_splitted[0], line_splitted[1]
            if rep not in cluster:
                cluster[rep] = {}
            print(rep, mem)#, ut[rep])
            if mem in ut[rep].keys():
                cluster[rep][mem] = ut[rep][mem]
            elif rep in ut[mem].keys():
                # TODO: how will we handle this?
                continue
                # cluster[line_splitted[0]][line_splitted[1]] = ut[line_splitted[1]][line_splitted[0]]
            # cluster[rep][mem] = ut[rep][mem]

            if mem not in entry_id:
                entry_id[mem] = last_id
                last_id +=1

    return cluster, entry_id

def open_structure(session, entry_id, pdb_path):
    from chimerax.core.fetch import fetch_file
    from chimerax.core.commands import run
    if pdb_path != "":
        run(session, "cd " + pdb_path)
    for k, v in entry_id.items():
        run(session, "open " + k + ".pdb id " + str(v))
    run(session, "hide atoms")
    run(session, "show cartoon")

def move_member(session, cluster, entry_id):
    from chimerax.core.commands import run
    for rep, mem in cluster.items():
        for m in mem.keys():
            rep_id = entry_id[rep]
            mem_id = entry_id[m]
            if rep_id != mem_id:
                print(mem_id, cluster[rep][m])
                run(session, "rotate_translate #" + str(mem_id) + " " + cluster[rep][m])

def rename_cluster(session, cluster, entry_id):
    from chimerax.core.commands import run
    for rep, mem in cluster.items():
        rep_id = entry_id[rep]
        rep_entry = rep.split(".pdb")[0]
        for m in mem.keys():
            mem_id = entry_id[m]
            mem_entry = m.split(".pdb")[0]
            if rep_id != mem_id:
                run(session, "rename #"+str(mem_id)+" "+rep_entry+":"+mem_entry)
            else:
                run(session, "rename #"+str(mem_id)+" "+rep_entry)

def color_cluster(session, cluster, entry_id):
    from chimerax.core.commands import run
    for rep, mem in cluster.items():
        rep_id = entry_id[rep]
        for m in mem.keys():
            mem_id = entry_id[m]
            if rep_id != mem_id:
                run(session, "color #"+str(mem_id)+" white")
                run(session, "transparency #"+str(mem_id)+" 50 cartoon")
            else:
                run(session, "color #"+str(mem_id)+" bychain")

def assign_cluster(session, cluster, entry_id):
    from chimerax.core.commands import run
    for rep, mem in cluster.items():
        rep_entry = rep.split(".pdb")[0]
        if len(mem) == 1:
            run(session, "name clu_"+rep_entry+" #"+str(entry_id[list(mem.keys())[0]]))
        mem_ids = ",".join([str(entry_id[str(m)]) for m in list(mem.keys())])
        run(session, "name clu_"+rep_entry+" #"+mem_ids)

def align_cluster(session, cluster_info, complexcluster_report, pdb_path=""):
    cluster, entry_id = parse_report(cluster_info, complexcluster_report)
    # cluster, entry_id = parse_cluster(complexcluster_result)
    open_structure(session, entry_id, pdb_path)
    move_member(session,  cluster, entry_id)
    rename_cluster(session, cluster, entry_id)
    color_cluster(session, cluster, entry_id)
    assign_cluster(session, cluster, entry_id)
    
def register_command(logger):
    from chimerax.core.commands import CmdDesc, register, OpenFileNameArg, StringArg, ListOf, FloatsArg, IntArg
    from chimerax.atomic import AtomsArg
    rt = CmdDesc(
        required = [('atoms', AtomsArg),
                    ('ut', FloatsArg)
                    ],
        synopsis = 'Rotate the structure with the given unitary matrix and translation vector.'
    )
    register('rotate_translate', rt, rotate_translate, logger=logger)

    show_clu = CmdDesc(
        required = [ ('cluster_info', OpenFileNameArg),
                    ('complexcluster_report', OpenFileNameArg),
                    ],
        optional=[ ("pdb_path", StringArg),
                # ("usecols", IntArg),
                ],
        synopsis = 'parse complexcluster result and show the cluster in ChimeraX',
    )

    register('show_cluster', show_clu, align_cluster, logger=logger)

print(sys.argv)
if len(sys.argv) <2 :
    register_command(session.logger)
else: 
    if sys.argv[1] == "show_cluster" and len(sys.argv) == 6:
        align_cluster(session, sys.argv[2], sys.argv[3], sys.argv[5])
