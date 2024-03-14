# Utility
- Author: Rachel Seongeun Kim
- Table of Contents
  1. [Tutorial](#tutorials)
  2. [ChimeraX Functions](#ChimeraX)
     1. [Foldseek Complex](#functions-specialized-for-foldseek-complex)

--- 

## Tutorials
- [ChimeraX tutorial(2023-12-01)](https://github.com/rachelse/lab_231201)
---

## ChimeraX
### Functions specialized for Foldseek Complex
 > Helper functions to visualize the result of complexsearch and complexcluster 
 [source code](./chimerax/visualize_foldseek_complex.py)   
 [example data](./chimerax/data/)
1. Functions
   - rotate_translate : rotate and translate the complex by providing the ut matrix
   - superpose_complex : open query and target, and superpose the complex by providing scorecomplex report file
   - show_cluster : open all pdb files in the cluster, and superpose members on the representative by providing scorecomplex report file and complexcluster result file
   - print_tmscore : print the TM-score of the complex by providing scorecomplex report file, query and target
   - print_ut : print the UT matrix of the complex by providing scorecomplex report file, query and target
   - print_report : print the report of the complex by providing scorecomplex report file, query and target

2. Features
   - cluster alias : use the cluster as a unit e.g. show specific cluster (clu_repname)

3. How to use
   ```sh
   open visualize_foldseek_complex.py

   # 1. rotate_translate
   open query.pdb
   open target.pdb
   rotate_translate #2 0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,100,200,300

   # 2. superpose_complex
   superpose_complex /path/to/scorecomplex_report query,target pdb_path /path/to/pdb_dir qt_id 1,2

   # 3. show_cluster
   show_cluster /path/to/scorecomplex_report /path/to/complexcluster_result pdb_path /path/to/pdb

   # 4. print functions
   print_tmscore /path/to/scorecomplex_report query,target
   print_ut /path/to/scorecomplex_report query,target
   print_report /path/to/scorecomplex_report query,target

   ----------
   ## Additional information
   # 1. cluster alias e.g. clu_7soy_1
   hide cartoon
   show clu_7soy_1 cartoon
   ```
4. Visualized examples from ChimeraX<br>
   show_cluster<br>
   <img src = ./image/showcluster.png width=50% height=50%>   
   show clu_7soy_1 cartoon<br>
   <img src = ./image/showcluster2.png width=50% height=50%>   
   superpose_complex<br>
   <img src = ./image/superpose_complex.png width=50% height=50%>   