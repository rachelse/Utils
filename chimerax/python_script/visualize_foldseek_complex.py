'''
File: visualize_foldseek_complex.py
Project: chimerax
File Created: 13th Mar 2024
Author: Rachel Seongeun Kim (seamustard52@gmail.com)
-----
Description: This script is to parse the scorecomplex report and extract the tmscore, unitary matrix and translation vector.
Copyright: Rachel Seongeun Kim
'''

"""
How to use
----------
open visualize_foldseek_complex.py
# 1. rotate and translate the complex by providing the ut matrix
# open query
open #1
#open target
open #2
rotate_translate #2 0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012

# 2. print u,t matrix and superpose the complex
superpose_complex /path/to/scorecomplex_report query,target pdb_path /path/to/pdb qt_id 1,2

# 3. show the cluster
open show_cluster.py
show_cluster /path/to/scorecomplex_report /path/to/complexcluster_result pdb_path /path/to/pdb
print_tmscore /path/to/scorecomplex_report query,target
print_tmscore /path/to/scorecomplex_report query:target

----------
Additional Information
----------------------
# 1. use the cluster as a unit e.g. show specific cluster (clu_7soy_1)
hide cartoon
show clu_7soy_1 cartoon
"""
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

def save_info(info, dict):
    dict[info[0]][info[1]] = {
        "qTM-score": float(info[4]),
        "tTM-score": float(info[5]),
        "ut": ",".join([info[6], info[7]])
    }

def parse_report(scorecomplex_report):
    sc_result = {}
    with open(scorecomplex_report, 'r') as f:
        lines = f.readlines()
        for line in lines:
            # #1: query #2: target #5: qTM-score #6: tTM-score #7: u (rotation) #8: t (translation)
            line_splitted = line.split("\t")

            if line_splitted[0] not in sc_result:
                sc_result[line_splitted[0]] = {}
                save_info(line_splitted, sc_result)
            else :
                if line_splitted[1] not in sc_result[line_splitted[0]]:
                    save_info(line_splitted, sc_result)
                elif sc_result[line_splitted[0]][line_splitted[1]]["tTM-score"] < float(line_splitted[5]):
                    save_info(line_splitted, sc_result)
    return sc_result

def parse_cluster(cluster_result):
    cluster = {}
    entry_id = {}
    last_id = 1
    with open(cluster_result, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line_splitted = line.strip().split("\t")
            if line_splitted[0] not in cluster:
                cluster[line_splitted[0]] = [line_splitted[1]]
            else:
                cluster[line_splitted[0]].append(line_splitted[1])
            if line_splitted[1] not in entry_id:
                entry_id[line_splitted[1]] = last_id
                last_id += 1
    return cluster, entry_id

def open_structure(session, entry_id, pdb_path):
    from chimerax.core.fetch import fetch_file
    from chimerax.core.commands import run
    if pdb_path != "":
        run(session, "cd " + pdb_path)
    for k, v in entry_id.items():
        run(session, "open " + k + " id " + str(v))
    run(session, "hide atoms")
    run(session, "show cartoon")

def move_member(session, report, cluster, entry_id):
    from chimerax.core.commands import run
    for rep, mem in cluster.items():
        if rep not in report:
            continue

        for m in mem:
            rep_id = entry_id[rep]
            mem_id = entry_id[m]
            if rep_id != mem_id:
                if rep not in report or m not in report[rep]:
                    continue
                print(mem_id, report[rep][m]["ut"])
                run(session, "rotate_translate #" + str(mem_id) + " " + report[rep][m]["ut"])

def rename_cluster(session, cluster, entry_id):
    from chimerax.core.commands import run
    for rep, mem in cluster.items():
        rep_id = entry_id[rep]
        rep_entry = rep.split(".pdb")[0]
        for m in mem:
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
        for m in mem:
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
            run(session, "name clu_"+rep_entry+" #"+str(entry_id[mem[0]]))
        mem_ids = ",".join([str(entry_id[m]) for m in mem])
        run(session, "name clu_"+rep_entry+" #"+mem_ids)

def align_cluster(session, scorecomplex_report, complexcluster_result, pdb_path=""):
    report = parse_report(scorecomplex_report)
    cluster, entry_id = parse_cluster(complexcluster_result)
    open_structure(session, entry_id, pdb_path)
    move_member(session, report, cluster, entry_id)
    rename_cluster(session, cluster, entry_id)
    color_cluster(session, cluster, entry_id)
    assign_cluster(session, cluster, entry_id)

def print_tmscore(session, scorecomplex_report, query_target):
    query, target = _pdb_handler(query_target)
    report = parse_report(scorecomplex_report)
    print(f"qTM-score: {report[query][target]['qTM-score']}, tTM-score: {report[query][target]['tTM-score']}")

def print_ut(session, scorecomplex_report, query_target):
    query, target = _pdb_handler(query_target)
    report = parse_report(scorecomplex_report)
    print(f"ut: {report[query][target]['ut']}")

def print_report(session, scorecomplex_report, query_target):
    query, target = _pdb_handler(query_target)
    with open(scorecomplex_report, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line_splitted = line.split("\t")
            if line_splitted[0] == query and line_splitted[1] == target:
                print(f"query: {query}, target: {target}, qTM-score: {line_splitted[4]}, tTM-score: {line_splitted[5]}, ut: {line_splitted[6]},{line_splitted[7]}")

def superpose_complex(session, scorecomplex_report, query_target, pdb_path="", qt_id=[1,1]):
    from chimerax.core.commands import run

    report = parse_report(scorecomplex_report)
    query, target = _pdb_handler(query_target)
    entry_id = {}
    entry_id[query] = qt_id[0]
    entry_id[target] = qt_id[1]
    open_structure(session, entry_id, pdb_path)
    run(session, "rotate_translate #"+str(entry_id[target])+" "+report[query][target]["ut"])
    
def register_command(logger):
    from chimerax.core.commands import CmdDesc, register, OpenFileNameArg, StringArg, ListOf, FloatsArg, IntArg
    from chimerax.atomic import AtomsArg

    rt = CmdDesc(
        required = [('atoms', AtomsArg),
                    ('ut', FloatsArg)
                    # ('scorecomplex_report', OpenFileNameArg),
                    ],
        synopsis = 'Rotate the structure with the given unitary matrix and translation vector.'
    )
    register('rotate_translate', rt, rotate_translate, logger=logger)

    show_clu = CmdDesc(
        required = [ ('scorecomplex_report', OpenFileNameArg),
                    ('complexcluster_result', OpenFileNameArg),
                    ],
        optional=[ ("pdb_path", StringArg),
                # ("usecols", IntArg),
                ],
        synopsis = 'parse complexcluster result and show the cluster in ChimeraX',
    )

    register('show_cluster', show_clu, align_cluster, logger=logger)

    tm = CmdDesc(
        required = [ ('scorecomplex_report', OpenFileNameArg),
                    ('query_target', ListOf(StringArg))
                    ],
        synopsis = 'parse complexcluster result and print the result of the query and target',
    )
    ut = CmdDesc(
        required = [ ('scorecomplex_report', OpenFileNameArg),
                    ('query_target', ListOf(StringArg))
                    ],
        synopsis = 'parse complexcluster result and print the result of the query and target',
    )
    result = CmdDesc(
        required = [ ('scorecomplex_report', OpenFileNameArg),
                    ('query_target', ListOf(StringArg))
                    ],
        synopsis = 'parse complexcluster result and print the result of the query and target',
    )
    register('print_tmscore', tm, print_tmscore, logger=logger)
    register('print_ut', ut, print_ut, logger=logger)
    register('print_report', result, print_report, logger=logger)

    superpose = CmdDesc(
        required = [ ('scorecomplex_report', OpenFileNameArg),
                    ('query_target', ListOf(StringArg))
                    ],
        optional=[ ("pdb_path", StringArg),
                ("qt_id", ListOf(IntArg)),
                # ("usecols", IntArg),
                ],
        synopsis = 'superpose the complex with the given ut matrix',
    )
    register('superpose_complex', superpose, superpose_complex, logger=logger)


register_command(session.logger)