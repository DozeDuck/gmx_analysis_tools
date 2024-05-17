#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:15:48 2023

@author: dozeduck
"""
import getopt
import sys
import re
import pandas as pd
# import plotly
import plotly.graph_objs as go
import plotly.io as pio
# for PCA
import numpy as np
# from rpy2.robjects import r
# for histogram dist plot
import plotly.figure_factory as ff
import argparse
# for free energy
import io
# for bool values
import ast
# metal restraints adding
import math
# for contact map
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import csv
import matplotlib.pyplot as plt
          
parser = argparse.ArgumentParser(description='Version: 1.6  \n'
                                 'Usage: python gmx_tool_box.py can generate figures based on *.xvg files; adding restraints for metals; merge lig_GMX.gro to rec.gro and modify topol.top automaticlly \n' \
                                  './gmx_tool_box -m <file1> <file2> <file3>  -o <output_name> -r <true/fault> -a <2/3/true> -x <xaxis_name> -y <yaxis_name> -c <rdf_cutoff> -p <plot_title> \n ' \
                                  './gmx_tool_box -m rmsd1.xvg rmsd2.xvg rmsd3.xvg -o rmsd123.png -a true \n -vi true'\
                                  './gmx_tool_box -m rmsf1.xvg rmsf2.xvg rmsf3.xvg -rn true -o rmsd123.png -a true \n '\
                                  './gmx_tool_box -m sasa1.xvg sasa2.xvg sasa3.xvg -o sasa123.png -rn true -xn "Time (ns)" -yn "SASA (nm<sup>2</sup>)" \n ' \
                                  './gmx_tool_box -m rdf1.xvg rdf2.xvg rdf3.xvg -rc 4.7 -o rdf123.png \n ' \
                                  './gmx_tool_box -m sasa1.xvg sasa2.xvg sasa3.xvg sasa4.xvg sasa5.xvg sasa6.xvg sasa7.xvg sasa8.xvg sasa9.xvg -o multy_sasa_9.png -a true -ma 10 -ml true \n '\
                                  './gmx_tool_box -m resarea1.xvg resarea2.xvg resarea3.xvg -o test_resare123.png -a true -r true \n ' \
                                  './gmx_tool_box -m 2dproj1.xvg -pn 1000 -pz 500 -o pca_myname.png \n ' \
                                  './gmx_tool_box -m dihedral-ave.xvg -o dihedral-ave.png -pn 1 -hg true -t "Dihedral F38@C-T39@N-T39@CA-T39@C" \n ' \
                                  './gmx_tool_box -m rmsd1.xvg rmsd2.xvg rmsd3.xvg -o error_bar.png -eb "error band" -rp 3 -tr 0.2 \n ' \
                                  '##################### Adding restraints to metal atoms ##################### \n' \
                                  './gmx_tool_box -mtgro rec.gro -mtml MG MN ZN CA -mtrl HIS GLU ASP ASN CYS LYS TYR -mtal ND1 OE1 OE2 OD1 OD2 ND2 SG NZ OH -mtd 0.4 -mtn 3 -mtbr 200000 -mtar 10000 \n ' \
                                  '##################### Merge receptor.gro and ligand_GMX.gro and edit topol.top ##################### \n' \
                                  './gmx_tool_box -gmgro rec.gro -gmlgro lig_GMX.gro -gmlitp lig_GMX.itp -gmtop topol.top \n'\
                                  '##################### Contact map detection ####################### \n' \
                                  './gmx_tool_box -cmp md_noPBC.pdb -cmf md_noPBC.pdb -cml LIG -cmo contact.png -cmd 3.5 \n' \
                                  '##################### Peptide format to ligand format #################\n' \
                                  './gmx_tool_box -plp pep.pdb -pln PEP\n' \
                                  '##################### gmx dssp ploting ###############\n' \
                                  './gmx_tool_box -dsf dssp.dat -dst md_noPBC_dt1000.pdb -dso dssp.png -dsx false -dsc true\n' \
                                  '##################### renumber extreme1.pdb MODEL count ###############\n' \
                                  './gmx_tool_box -rnf extreme1.pdb\n' \
                                  '##################### add new residue to gromacs forcefield ###############\n' \
                                  'gmx_tool_box -arlt lig_GMX.itp -arat N CT C O CT S CT CT N3 C O CT CT NX NX O O O O H1 H1 H1 H1 H1 H1 H1 H1 H1 H -arfb $HOME/workspace/gromacs2023.2-plumed/share/gromacs/top/amber99sb-ildn.ff/ffbonded.itp  -arfa amber99sb-ildn-cr1-pga-tpo-sep-lyseba-atp-adp-crystalWater.ff/atomtypes.atp  -arfnb amber99sb-ildn-cr1-pga-tpo-sep-lyseba-atp-adp-crystalWater.ff/ffnonbonded.itp',
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-m', '--multi_files', nargs='+', required=False, default=[0], help='input gro file')
parser.add_argument('-o', '--output_name', default='lack-of-name.png', help='output_name.png, suitable file format for output: png, jpg, jpeg, webp, svg, pdf, eps, json')
parser.add_argument('-rn', '--renumber', default='false', help='true: default is fault, renumber the residues, mainly for the duplicated residue numbers while work on rmsf and sasa per residue')
parser.add_argument('-rc', '--rdf_cutoff', default=0, help='number: rdf distance cutoff value')
parser.add_argument('-a', '--average',  default='false', help='true: default is fault, output the average score for each replica')
parser.add_argument('-t', '--plot_name', default=0, help='plot_name: the title showed in plot')
parser.add_argument('-pn', '--nbin', default=1, help='represent "nbin", mainly used for pca ploting, default value=1, the larger the smoother, can set to 1000, however, if there is white line in your PCA_Density plot, please change the value(-n 500) to get a cleaner plot')
parser.add_argument('-pz', '--size', nargs='+',  default=500, help='represent "size", mainly used for pca ploting, default value=500, the larger the higher resolution, if there is white line in your PCA Density plot, please change the value to get a cleaner plot')
parser.add_argument('-ma', '--move_average', default=0, help='the window size of moving average analysis for SASA, can be set to 2,3,4,5,6,...')
parser.add_argument('-ml', '--mean_value', default='false', help='whether user want to draw a straight line to represent the mean value for each trace. default value = false; user can change it to true')
parser.add_argument('-hg', '--histogram', default='false', help='whether user want to generate histogram plot, default value = false')
parser.add_argument('-xn', '--xaxis_name', default=0, help='User indicate the name of X axis, e.g: Time (ns)')
parser.add_argument('-yn', '--yaxis_name', default=0, help='User indicate the name of Y axis, e.g: SASA (nm<sup>2</sup>)')
parser.add_argument('-xs', '--xaxis_size', default=800, help='User indicate the size of width, default = 800')
parser.add_argument('-ys', '--yaxis_size', default=600, help='User indicate the size of Y height, default = 600')
parser.add_argument('-xyf', '--xy_font', default=40, help='User indicate the size of x and y font, default = 40')

parser.add_argument('-xl', '--xaxis_low', default=0, help='User indicate the range of xaxis_low, default = 0')
parser.add_argument('-xh', '--xaxis_high', default=0, help='User indicate the range of xaxis_high, default = 0')
parser.add_argument('-yl', '--yaxis_low', default=0, help='User indicate the range of yaxis_low, default = 0')
parser.add_argument('-yh', '--yaxis_high', default=0, help='User indicate the range of yaxis_high, default = 0')

parser.add_argument('-eb', '--error_bar', default='false', help='whether use "error bar" or "error band", default value = false')
parser.add_argument('-rp', '--replica_number', default=3, help='Number of replicas, default = 3')
parser.add_argument('-tr', '--transparency', default=0.2, help='Number of transparent (for error band), default = 0.2')

parser.add_argument('-tf', '--title_font', default=24, help='User indicate the size of plot title font, default = 24')
parser.add_argument('-ls', '--legend_show', default='True', help='Whether user want to show the legend, default = True')
parser.add_argument('-lf', '--legend_font', default=30, help='Define the legend font size, defaul = 30')
parser.add_argument('-fm', '--font_family', default='Arial', help='Font family default = Arial')
parser.add_argument('-fc', '--font_color', default='black', help='Font color default = black')
parser.add_argument('-gs', '--grid_show', default='True', help='Whether show the grid line, defaul = True')
parser.add_argument('-vi', '--violin', default='False',required=False, help='Whether show the data in violin style, defaul = False')
# define the parameters for metal restraints adding
parser.add_argument('-mtgro', '--mtgros', required=False, default=0, help='input gro file')
parser.add_argument('-mtml', '--mtmetallist', nargs='+', default=["MG","MN","ZN","CA"], help='the list of Metals,default value is MG,MN,ZN,CA')
parser.add_argument('-mtrl', '--mtreslist', nargs='+', default=["HIS","GLU","ASP","ASN", "CYS", "LYS", "TYR"], help='the list of residues,default value is "HIS","GLU","ASP","ASN", "CYS", "LYS", "TYR"')
parser.add_argument('-mtal', '--mtatomlist', nargs='+', default=["ND1","OE1","OE2", "OD1","OD2","ND2","SG", "NZ", "OH"], help='the list of atom name,default value is "ND1","OE1","OE2", "OD1","OD2","ND2","SG", "NZ", "OH"')
parser.add_argument('-mtd', '--mtdistance', default=0.4, help='the distance default = 0.4 nm')
parser.add_argument('-mtn', '--mtneighbours', default=3, help='the number of neighbours default value is 3 neighbours')
parser.add_argument('-mtbr', '--mtbond', default=200000, help='give the path for your em.mdp file, or this script will use the default_em.mdp')
parser.add_argument('-mtar', '--mtangle', default=10000, help='give the path for your nvt.mdp file, or this script will use the default_nvt.mdp')
# define the parameters for gromerger
parser.add_argument('-gmgro', '--gmreceptor_gro', required=False, default=0, help='The receptor gro file generated from gmx pdb2gmx')
parser.add_argument('-gmlgro', '--gmligand_gro', required=False, default='lig_GMX.gro', help='The ligand gro file generated from acpype')
parser.add_argument('-gmlitp', '--gmligand_itp', required=False, default='lig_GMX.itp', help='The ligand itp file generated from acpype')
parser.add_argument('-gmtop', '--gmreceptor_top', required=False, default='topol.top', help='The protein topol file, default = topol.top')
# define the parameters for ligand contact map detect
parser.add_argument('-cmp', '--cm_topol', required=False, default=0, help='The topology file, can be in md_noPBC.pdb, topol.top, md.tpr file')
parser.add_argument('-cmf', '--cm_traj', required=False, default='md_noPBC.pdb', help='The trajectory file, canbe md_noPBC.pdb, md.xtc, md.trr')
parser.add_argument('-cml', '--cm_lig', required=False, default='LIG', help='The residue name of ligand')
parser.add_argument('-cmo', '--cm_output', required=False, default='ligand_contact.png', help='The output name for the contact map')
parser.add_argument('-cmd', '--cm_distance', required=False, default=3.5, help='The distance threshold for contact map, in angstrom, default value is 3.5 A')
# define the parameters for petide to ligand
parser.add_argument('-plp', '--pl_pdb', required=False, default=0, help='The peptide alone PDB file')
parser.add_argument('-pln', '--pl_name', required=False, default='PEP', help='The resname for peptide after converting to ligand')
# define the parameter for gmx dssp anaysis, usage: python gmx_tool_box.py -dsf dssp.dat -dst md_noPBC_dt1000.pdb -dso dssp.png -dsx false -dsc true
parser.add_argument('-dsf', '--ds_dat', required=False, default=0, help='input dssp.dat file')
parser.add_argument('-dst', '--ds_traj', required=False, default='md_noPBC_dt1000.pdb', help='the traj  file')
parser.add_argument('-dso', '--ds_output', required=False, default='dssp.png', help='output name')
parser.add_argument('-dsx', '--ds_original', required=False, default='false', help='whether use original structure types')
parser.add_argument('-dsc', '--ds_color', required=False, default='true', help='whether generate a unique color bar')
# define the parameer for renumber the MODEL count in extrem1.pdb 
parser.add_argument('-rnf',  '--rn_file', required=False, default=0, help='input extreme1.pdb file')
# define the parameter for adding new residues to forcefield
parser.add_argument('-arlt', '--arlig_itp', default=0, required=False, help='the file generated by acpype -i lig.mol2 -c user; lig_GMX.itp')
parser.add_argument('-arlg', '--arlig_gro', default='lig_GMX.gro', required=False, help='the file generated by acpype -i lig.mol2 -c user; lig_GMX.gro')
parser.add_argument('-arfr', '--arff_rtp', default='aminoacids.rtp', required=False, help='the file from $GMX/share/gromacs/top/amber99sb-ildn/aminoacids.rtp')
parser.add_argument('-arfh', '--arff_hdb', default='aminoacids.hdb', required=False, help='the file from $GMX/share/gromacs/top/amber99sb-ildn/aminoacids.hdb')
parser.add_argument('-arfb', '--arff_bonded',  default='ffbonded.itp', required=False, help='the file from $GMX/share/gromacs/top/amber99sb-ildn/ffbonded.itp')
parser.add_argument('-arfnb', '--arff_nonbonded',  default='ffnonbonded.itp', required=False, help='the file from $GMX/share/gromacs/top/amber99sb-ildn/ffnonbonded.itp')
parser.add_argument('-arfa', '--arff_atomtypes',  default='atomtypes.atp', required=False, help='the file from $GMX/share/gromacs/top/amber99sb-ildn/atomtypes.atp')
parser.add_argument('-arfrt', '--arff_restype',  default='residuetypes.dat', required=False, help='the file from $GMX/share/gromacs/top/residuetypes.dat')
parser.add_argument('-aria', '--arisamino_acid', default='true', required=False, help='whether treat it as an aminoacid residue, default value = true')
parser.add_argument('-arat', '--aratom_types', nargs='+', required=False, default=['a'], help='mapped atomtypes')
parser.add_argument('-aro', '--aroutput_name', default='output', help='the output file name')


args = parser.parse_args()

# define the parameters
multi_files     = [x for x in args.multi_files]
output_name     = str(args.output_name)
renumber        = str(args.renumber)
rdf_cutoff      = int(args.rdf_cutoff)
average         = str(args.average)
plot_name       = str(args.plot_name)
nbin            = int(args.nbin)
size            = int(args.size)
move_average    = int(args.move_average)
mean_value      = str(args.mean_value)
histogram       = str(args.histogram)
xaxis_name      = str(args.xaxis_name)
yaxis_name      = str(args.yaxis_name)

x_low           = float(args.xaxis_low)
x_high          = float(args.xaxis_high)
y_low           = float(args.yaxis_low)
y_high          = float(args.yaxis_high)

error_bar       = str(args.error_bar)
replica_number  = int(args.replica_number)
transparency    = float(args.transparency)

xaxis_size      = int(args.xaxis_size)
yaxis_size      = int(args.yaxis_size)
xy_font         = int(args.xy_font)
title_font      = int(args.title_font)
legend_show     = ast.literal_eval(args.legend_show)
legend_font     = int(args.legend_font)
font_family     = str(args.font_family)
font_color      = str(args.font_color)
grid_show       = ast.literal_eval(args.grid_show)
violin          = str(args.violin)
# define the parameters for metal restraints adding
mr_file             = args.mtgros
mr_metal_list       = [str(x) for x in args.mtmetallist]
mr_residue_list     = [str(x) for x in args.mtreslist]
mr_atom_list        = [str(x) for x in args.mtatomlist]
mr_distance_value   = args.mtdistance
mr_num_neighbours   = int(args.mtneighbours)
mr_bond_strength    = args.mtbond
mr_angle_strength   = args.mtangle
# define the parameters for gromerger
gm_receptor_gro         = args.gmreceptor_gro
gm_ligand_gro           = str(args.gmligand_gro)
gm_ligand_itp           = str(args.gmligand_itp)
gm_receptor_top         = str(args.gmreceptor_top)
# define the parameters for ligand contact map
cm_topol            = args.cm_topol
cm_traj             = str(args.cm_traj)
cm_lig              = str(args.cm_lig)
cm_output           = str(args.cm_output)
cm_distance         = float(args.cm_distance)
# define the parameters for petide to ligand
pl_pdb              = args.pl_pdb
pl_name             = str(args.pl_name)
# define the parameter for gmx dssp figure
ds_dssp = args.ds_dat
ds_traj = args.ds_traj
ds_output = args.ds_output
ds_original = args.ds_original
ds_color = args.ds_color
# define the parameter for renumber the MODEL count in extreme1.pdb
rn_file = args.rn_file
# define the parameter for adding new residues to forcefield
ar_lig_itp         = str(args.arlig_itp)
ar_lig_gro         = str(args.arlig_gro)
ar_ff_rtp          = str(args.arff_rtp)
ar_ff_hdb          = str(args.arff_hdb)
ar_ff_bonded       = str(args.arff_bonded)
ar_ff_nonbonded    = str(args.arff_nonbonded)
ar_ff_aomtypes     = str(args.arff_atomtypes)
ar_ff_restype      = str(args.arff_restype)
ar_isamino_acid      = str(args.arisamino_acid)
ar_atom_types      = [x for x in args.aratom_types]
ar_output_name     = str(args.aroutput_name)



class plotly_go():
    flag = ''
    sasa_flag = ''
    pca_flag = ''
    time1 = []
    values1 = []
    sd1 = []
    time2 = []
    values2 = []
    sd2 = []
    time3 = []
    values3 = []
    sd3 = []
    max_value = []
    average_value = []
    multi_flag = ''

    def __init__(self, multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, violin, x_low, x_high, y_low, y_high, error_bar, replica_number, transparency):

        if len(multi_files) >=1:
            # print(multi_files)
            file1 = multi_files[0]
            self.flag_recognizer(file1)
            if self.pca_flag != 1 and self.flag != 'pca' and self.flag != 'free energy':
                self.plotly_multy(multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, self.flag, violin, x_low, x_high, y_low, y_high, error_bar, replica_number, transparency)
            elif self.pca_flag == 1:
                self.plotly_pca(multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, self.flag, x_low, x_high, y_low, y_high, error_bar, replica_number, transparency)
            elif self.flag == 'pca':
                self.plotly_pca(multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, self.flag, x_low, x_high, y_low, y_high, error_bar, replica_number, transparency)
            elif self.flag == 'free energy':
                self.plotly_free_energy(multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, self.flag, violin, x_low, x_high, y_low, y_high, error_bar, replica_number, transparency)





    def flag_recognizer(self,file1):                                                   # first method to be called in __main__, used for creating object and charactors.
        flags_map = {
            'rms,': 'rmsd',
            'rmsf,': 'rmsf',
            'sasa,': 'sasa',
            'gyrate,': 'gyrate',
            'dipoles,': 'dipoles',
            'distance,': 'distance',
            'rdf,': 'rdf',
            'convergence': 'convergence',
            'anaeig,': 'pca',
            'angle,': 'angle',
            'free': 'free energy'
        }
        if file1.endswith(".xvg"):
            with open(file1, 'r') as f:
                lines = f.readlines()
                if len(lines) >= 3:
                    try:
                        flag = lines[2].split()[5]
                        self.flag = flags_map.get(flag, flag)
                    except IndexError:
                        pass
                if len(lines) >= 9 and '-or' in lines[8]:
                    self.sasa_flag = '-or'
    
                if 'pca' in str(file1).lower() or '2dproj' in str(file1):
                    self.pca_flag = 1
                print("I know you are plotting " + self.flag + " figures!")
        elif file1.endswith(".csv"):
            found = False
            for key in flags_map:
                if key.strip(',') in plot_name.lower():  
                    found = True
                    self.flag = flags_map[key]  
                    break  
        
        elif file1.endswith(".dat"):
            with open(file1, 'r') as f:
                lines = f.readlines()            
                first_line = lines[0]
                flag = 'free' if 'free' in first_line else None   
                self.flag = flags_map.get(flag, flag)   
            print("I know you are plotting " + self.flag + " figures!")


    def consist(self,x_values):
        seq1 = [x_values[0]]
        seq2 = []
        # seq3 = []
        for i in range(1, len(x_values)):
            # find the break point
            if x_values[i] <= x_values[i-1]+2 and seq2 == []:
                seq1.append(x_values[i])
            else:
                seq2.append(x_values[i])
        return seq1, seq2


    def read_data(self, file, x_name, renumber):
        # 从文件中读取数据
        x_data, y_data, sd_data = [], [], []
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("#") or line.startswith("@"):
                    continue
                else:
                    # 解析数据行
                    split_line = line.split()
                    x_value = float(split_line[0])
                    y_value = float(split_line[1])

                    if x_name == 'Time (ps)':  # 将时间从ps转换为ns
                        x_value /= 1000

                    if x_name == 'Residue' and renumber == 'true':
                        x_value = len(x_data) + 1

                    x_data.append(x_value)
                    y_data.append(y_value)

                    # 读取标准差（如果存在）
                    try:
                        sd_data.append(float(split_line[2]))
                    except IndexError:
                        pass
        return x_data, y_data, sd_data
    
    def read_data_dat(self, file_name):
        with open(file_name, 'r') as file:
            lines = file.readlines()
        # Identifying the line with column names
        columns_line = [line for line in lines if line.startswith('#! FIELDS')][0]
        # Extracting column names
        column_names = columns_line.strip().split()[2:]  # Skip '#! FIELDS'
        # 判断free_energy data 在第几列
        index_of_free_energy = [index for index, name in enumerate(column_names) if 'free' in name][0]
        # Extracting data lines (those not starting with '#')
        data_lines = [line for line in lines if not line.startswith('#')]
        # Converting data lines to a pandas DataFrame
        df = pd.read_csv(io.StringIO('\n'.join(data_lines)), delim_whitespace=True, names=column_names)
        
        if index_of_free_energy == 2:
            x_data = df[df.columns[0]].tolist()  # 第一列作为X轴数据
            y_data = df[df.columns[1]].tolist()  # 第2列作为y轴数据
            z_data = df[df.columns[2]].tolist()  # 第3列作为z轴数据
        elif index_of_free_energy == 1:
            x_data = df[df.columns[0]].tolist()  # 第一列作为X轴数据
            y_data = df[df.columns[1]].tolist()  # 第2列作为y轴数据
            z_data = []

        # CSV文件不包含标准差数据，因此sd_data保持为空 
                   
        return x_data, y_data, z_data, df, index_of_free_energy, column_names
    
    def extract_plot_details(self, multi_files, plot_name, xaxis_name, yaxis_name, flag, histogram):
        regex = r"\[|\]|'"
        
        # 提取或设置图表标题
        if plot_name == '0':
            with open(multi_files[0], "r") as f:
                plot_title = re.sub(regex, "", str(re.findall('"([^"]*)"', f.readlines()[13])))
        else:
            plot_title = str(plot_name)

        # 提取或设置X轴名称
        if xaxis_name == '0':
            with open(multi_files[0], "r") as f:
                x_name = re.sub(regex, "", str(re.findall('"([^"]*)"', f.readlines()[14])))
        else:
            x_name = xaxis_name

        # 提取或设置Y轴名称
        if yaxis_name == '0':
            with open(multi_files[0], "r") as f:
                y_name = re.sub(regex, "", str(re.findall('"([^"]*)"', f.readlines()[15])))
            if plot_title in ['Solvent Accessible Surface', 'Area per residue over the trajectory']:
                y_name = 'Area (nm<sup>2</sup>)'
            elif flag == 'dihedral_distribution' and histogram == 'true':
                y_name = 'Probability'
        else:
            y_name = yaxis_name

        return plot_title, x_name, y_name

    def define_trace(self, x_data, y_data, file_name, colour, violine='False', flag=0, labels=0):
        # 创建并返回迹线
        if flag == 'pca':
            trace = go.Scatter(
                x=x_data,
                y=y_data,
                mode='markers',
                marker=dict(
                    color=labels,  # 设置颜色为标签的数值
                    colorscale=colour,  # 颜色映射，你可以根据需要选择不同的颜色映射
                    colorbar=dict(title='Label Range'),  # 添加颜色条
                ),
            )
        elif violine != 'False':
            trace = go.Violin(x0=str(file_name).split('.')[0], y=y_data, line=dict(color='black'), fillcolor=colour, name=str(file_name).split('.')[0], box_visible=True, meanline_visible=True, opacity=0.6)
        else:
            trace = go.Scatter(x=x_data, y=y_data, line=dict(color=colour), name=str(file_name).split('.')[0])
        return trace

    def calculate_for_error_bar_or_band(self, multi_files, x_name, replica_number, uploaded_filenames):
        df_data = pd.DataFrame()
        df_average = pd.DataFrame()
        df_sd   = pd.DataFrame()
        count = 1
        if multi_files[0].endswith(".xvg"):
            for i, file in enumerate(multi_files):
                x_data, y_data, _ = self.read_data(file, x_name, renumber)
                if i == 0:
                    x_datas = x_data
                df_data[f'y_data_{i+1}'] = y_data
                if i == (count * replica_number) - 1:  # 检查是否达到组内文件数量
                    # 计算当前df_data的所有列的平均值和标准差
                    mean_vals = df_data.mean(axis=1)
                    std_vals = df_data.std(axis=1)
                    # 将计算得到的平均值和标准差添加到相应的DataFrame中
                    df_average[uploaded_filenames[(count-1)*3]] = mean_vals
                    df_sd[uploaded_filenames[(count-1)*3]] = std_vals
                    # 重置df_data以便下一组的使用，并更新计数器
                    df_data = pd.DataFrame()
                    count += 1
        elif multi_files[0].endswith(".csv"):
            for i, file in enumerate(multi_files):
                x_data, y_data, _ = self.read_data_csv(file, x_name, renumber) 
                if i == 0:
                    x_datas = x_data
                df_data[f'y_data_{i+1}'] = y_data
                if i == (count * replica_number) - 1:  # 检查是否达到组内文件数量
                    # 计算当前df_data的所有列的平均值和标准差
                    mean_vals = df_data.mean(axis=1)
                    std_vals = df_data.std(axis=1)
                    # 将计算得到的平均值和标准差添加到相应的DataFrame中
                    df_average[uploaded_filenames[(count-1)*3]] = mean_vals
                    df_sd[uploaded_filenames[(count-1)*3]] = std_vals
                    # 重置df_data以便下一组的使用，并更新计数器
                    df_data = pd.DataFrame()
                    count += 1
        elif multi_files[0].endswith(".dat"):
            for i, file in enumerate(multi_files):
                x_data, y_data, z_data, df, index_of_free_energy, column_names = self.read_data_dat(file, x_name, renumber) 
                if i == 0:
                    x_datas = x_data
                df_data[f'y_data_{i+1}'] = y_data
                if i == (count * replica_number) - 1:  # 检查是否达到组内文件数量
                    # 计算当前df_data的所有列的平均值和标准差
                    mean_vals = df_data.mean(axis=1)
                    std_vals = df_data.std(axis=1)
                    # 将计算得到的平均值和标准差添加到相应的DataFrame中
                    df_average[uploaded_filenames[(count-1)*3]] = mean_vals
                    df_sd[uploaded_filenames[(count-1)*3]] = std_vals
                    # 重置df_data以便下一组的使用，并更新计数器
                    df_data = pd.DataFrame()
                    count += 1
        return df_average, df_sd, x_datas

    def hex_to_rgba(self, hex_color, alpha=0.2):
        # 移除可能的 "#" 符号
        hex_color = hex_color.lstrip('#')
        # 通过列表推导从十六进制字符串中提取并转换为RGB整数值
        rgb = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
        # 将RGB整数值和透明度alpha组合成RGBA字符串
        return f"rgba({rgb[0]}, {rgb[1]}, {rgb[2]}, {alpha})"

    def define_trace_for_error_bands(self, error_bar, df_average, df_sd, x_data, transparency):
        Plotly = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
        colors = ['rgb(0,100,80)', 'rgb(255,0,0)']  # 不同组使用不同颜色
        traces = []
        if error_bar =='error band':
            for idx, (col_name_avg, col_name_sd) in enumerate(zip(df_average.columns, df_sd.columns)):
                y = df_average[col_name_avg]
                y_std = df_sd[col_name_sd]
                y_upper = y + y_std
                y_lower = y - y_std
                fill_color = self.hex_to_rgba(Plotly[idx], alpha=transparency)
                
                traces.append(go.Scatter(
                    x=x_data,
                    y=y,
                    line=dict(color=Plotly[idx]),
                    mode='lines',
                    name=col_name_avg  # 使用列名作为轨迹名称
                ))
                traces.append(go.Scatter(
                    x=list(x_data) + list(x_data)[::-1],
                    y=list(y_upper) + list(y_lower)[::-1],
                    fill='toself',
                    fillcolor=fill_color,
                    line=dict(color='rgba(255,255,255,0)'),
                    hoverinfo="skip",
                    showlegend=False
                ))
        elif error_bar == 'error bar':
            for idx, (col_name_avg, col_name_sd) in enumerate(zip(df_average.columns, df_sd.columns)):
                y = df_average[col_name_avg]
                y_std = df_sd[col_name_sd]
                traces.append(go.Scatter(
                    x=x_data,
                    y=df_average[col_name_avg],
                    error_y=dict(type='data', array=df_sd[col_name_sd], visible=True)
                ))
        return traces


    def setup_layout(self, plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show, x_low, x_high, y_low, y_high, violine='False', flag=0):
        # 设置布局
        if flag == 'pca':
            legend_show = False
        if violine != 'False':
            x_name = ''
        if y_low == y_high == 0:
            y_autorange = True
        else:
            y_autorange = False
        if x_low == x_high == 0:
            x_autorange = True
        else:
            x_autorange = False
            
        layout = go.Layout(
            title=plot_title, title_x=0.5, title_y=1, font=dict(size=title_font, color=font_color),
            xaxis=dict(title=x_name, titlefont=dict(size=xy_font, color=font_color, family=font_family), zeroline=False, autorange=x_autorange, range=[x_low,x_high],
                       showgrid=grid_show, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
            yaxis=dict(title=y_name, titlefont=dict(size=xy_font, color=font_color, family=font_family), zeroline=False, autorange=y_autorange, range=[y_low,y_high],
                       showgrid=grid_show, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
            legend=dict(x=1, y=1, orientation='v', font=dict(size=legend_font, color=font_color)), showlegend=legend_show,
            plot_bgcolor='rgba(255, 255, 255, 0.1)',
            paper_bgcolor='rgba(255, 255, 255, 0.2)',
            width=xaxis_size, height=yaxis_size
        )
        return layout
   
    def plot_graph(self, data, layout, output_file_name):
        # 使用数据和布局绘制图形
        fig = go.Figure(data=data, layout=layout)
        pio.write_image(fig, output_file_name)

    def plot_histogram(self, histogram_data, group_labels, plot_title, output_file_name, colors, nbin):
        # 处理直方图
        fig_hist = ff.create_distplot(histogram_data, group_labels, colors=colors, bin_size=nbin, curve_type='normal')
        fig_hist.update_layout(title_text=plot_title)
        pio.write_image(fig_hist, "histogram_" + output_file_name)

    def calculate_average(self, multi_files, xaxis_name, renumber):
        # 计算平均值
        sum_data = None
        for file in multi_files:
            x_data, y_data, _ = self.read_data(file, xaxis_name, renumber)
            if sum_data is None:
                sum_data = np.array(y_data)
            else:
                sum_data += np.array(y_data)
        return sum_data / len(multi_files)
    
    def output_average_file(self, output_file_name, average_value, multi_files, xaxis_name, renumber, x_data):
        with open(output_file_name[:-4]+"_average.xvg", 'w') as f:
            with open(multi_files[0], "r") as a:
                lines = a.readlines()
                for num in range(len(lines)):
                    if lines[num].startswith("#") or lines[num].startswith("@"):
                        f.write(lines[num])
                    else:
                        pass
            for num in range(len(average_value)):
                average_line = "{}     {}\n".format(x_data[num], average_value[num])
                f.write(average_line)

    def moving_average(self, y_data, window_size):
        # 计算移动平均
        return np.convolve(y_data, np.ones(window_size) / window_size, mode='valid')


    def plotly_multy(self, multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, flag, violin, x_low, x_high, y_low, y_high, error_bar, replica_number, transparency):
        Plotly = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
        data, histogram_data, group_labels = [], [], []

        # 读取plot_title, x_name, y_name
        plot_title, x_name, y_name = self.extract_plot_details(multi_files, plot_name, xaxis_name, yaxis_name, flag, histogram)

        # 读取数据并创建迹线
        for i, file in enumerate(multi_files):
            x_data, y_data, _ = self.read_data(file, xaxis_name, renumber)
            trace = self.define_trace(x_data, y_data, file, Plotly[i % len(Plotly)], violine=violin)
            data.append(trace)

            # 添加直方图数据
            if histogram == 'true':
                histogram_data.append(y_data)
                group_labels.append(str(file).split('.')[0])

        # 设置布局
        layout = self.setup_layout(plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show, x_low, x_high, y_low, y_high, violine=violin)

        # 绘制图形
        self.plot_graph(data, layout, output_name)

        # 处理直方图
        if histogram == 'true':
            self.plot_histogram(histogram_data, group_labels, plot_title, output_name, Plotly, nbin)

        # 处理平均值
        if average == 'true':
            average_data = self.calculate_average(multi_files, xaxis_name, renumber)
            average_trace = self.define_trace(x_data, average_data, "Average", 'black')
            # data.append(average_trace)
            data = average_trace
            self.plot_graph(data, layout, "Average_" + output_name)
            self.output_average_file(output_name, average_data, multi_files, xaxis_name, renumber, x_data)

        # 处理移动平均
        if move_average != 0:
            ma_data = []
            for file in multi_files:
                _, y_data, _ = self.read_data(file, xaxis_name, renumber)
                ma_y_data = self.moving_average(y_data, move_average)
                ma_trace = self.define_trace(x_data[move_average - 1:], ma_y_data, str(file).split('.')[0], Plotly[0])
                ma_data.append(ma_trace)
            self.plot_graph(ma_data, layout, "MovingAverage_" + output_name)
        
        # 处理error band or error bar
        if error_bar != 'false':
            df_average, df_sd, x_data = self.calculate_for_error_bar_or_band(multi_files, x_name, replica_number, multi_files)
            error_data = self.define_trace_for_error_bands(error_bar, df_average, df_sd, x_data, transparency)
            error_layout = self.setup_layout(plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show, x_low, x_high, y_low, y_high, violine=violin)
            self.plot_graph(error_data, error_layout, "error_bar_" + output_name)

        # # 调用上述方法
        # for i, file in enumerate(multi_files):
        #     x_data, y_data, sd_data = self.read_data(file, xaxis_name, renumber)
        #     trace = self.define_trace(x_data, y_data, file, Plotly[i % len(Plotly)])
        #     data.append(trace)

        # layout = self.setup_layout(plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show)
        # self.plot_graph(data, layout, output_name)


         
            
    def plotly_pca(self, multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, flag, x_low, x_high, y_low, y_high):
        data = []
        color = ['rainbow']
        # labels = []
        # 使用 extract_plot_details 方法获取图表标题和轴标签
        plot_title, x_name, y_name = self.extract_plot_details(multi_files, plot_name, xaxis_name, yaxis_name, flag, histogram)

        # 处理 PCA 数据
        for i, file in enumerate(multi_files):          
            x_data, y_data, _ = self.read_data(file, "PC1", renumber)  # 假设 "PC1" 和 "PC2" 是合适的轴名称
            labels = [x for x in range(len(y_data))]
            
            # 使用 define_trace 创建迹线
            trace = self.define_trace(x_data, y_data, file, 'rainbow', flag, labels=labels)  # 假设使用 'rainbow' 作为颜色
            data.append(trace)

        # 使用 setup_layout 设置布局
        layout = self.setup_layout(plot_title, title_font, 'PC1 (nm)', 'PC2 (nm)', xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show,  x_low, x_high, y_low, y_high, flag=flag)

        # 使用 plot_graph 绘制图形
        self.plot_graph(data, layout, "Scatter_" + output_name)

    def plotly_free_energy(self, multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, flag, violin, x_low, x_high, y_low, y_high, error_bar, replica_number, transparency):
        Plotly = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
        data, histogram_data, group_labels = [], [], []
        plot_title = 'Free Energy Surface'
        for i, file in enumerate(multi_files):
            x_data, y_data, z_data, df, index_of_free_energy, column_names = self.read_data_dat(file)
            # 如果有3列，则为phi psi 自由能
            if index_of_free_energy == 2:
                if 'phi' in column_names or 'psi' in column_names:
                    x_values = np.degrees(np.unique(x_data))
                    y_values = np.degrees(np.unique(y_data))
                else:
                    x_values = np.unique(x_data)
                    y_values = np.unique(y_data)                   
                x_grid, y_grid = np.meshgrid(x_values, y_values)
                z_data_array = np.array(z_data)
                free_energy_grid = z_data_array.reshape(len(y_values), len(x_values))
                # Plot the FES
                plt.contourf(x_grid, y_grid, free_energy_grid, levels=100)
                if xaxis_name == '0' and yaxis_name == '0':
                    xaxis_name = column_names[0]
                    yaxis_name = column_names[1]
                    plt.xlabel(xaxis_name)
                    plt.ylabel(yaxis_name)
                else:
                    plt.xlabel(xaxis_name)
                    plt.ylabel(yaxis_name)                
                plt.colorbar(label='Free energy / kJ mol^-1')
                plt.title('Free Energy Surface')
                plt.savefig(str(i) + "_" + output_name)
                # self.streamlit_download_file_plotly(str(i) + "_" + output_name, "/tmp/" + str(i) + "_" + output_name)
                # 3D plot
                fig=plt.figure(figsize=(24,14), dpi=600)
                ax = plt.axes(projection='3d')
                surf = ax.plot_surface(x_grid, y_grid, free_energy_grid, cmap = 'jet', rstride=1, cstride=1, alpha=None, linewidth=0, antialiased=True)
                # Set axes label
                ax.set_title('Free Energy Surface')
                ax.set_xlabel(xaxis_name, labelpad=5)
                ax.set_ylabel(yaxis_name, labelpad=5)
                ax.set_zlabel('Free energy / kJ mol^-1', labelpad=5)
                fig.colorbar(surf, shrink=0.7, aspect=15)
                plt.savefig(str(i) + "_3D_" + output_name)
                # self.streamlit_download_file_plotly(str(i) + "_3D_" + output_name, "/tmp/" + str(i) + "_3D_" + output_name)


            # 如果有2列，则为distance 自由能
            elif index_of_free_energy == 1:
                if xaxis_name == 'auto detect':
                    x_name = column_names[0]
                else:
                    x_name = xaxis_name
                if yaxis_name == 'auto detect':
                    y_name = 'Free energy / kJ mol<sup>-1</sup>'
                else:
                    y_name = yaxis_name
                # 使用 define_trace 创建迹线
                trace = self.define_trace(x_data, y_data, multi_files[i], Plotly[i % len(Plotly)], violine=violin)
                data.append(trace)
        if data != []:
           # 设置布局
            layout = self.setup_layout(plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show, x_low, x_high, y_low, y_high, violine=violin)
            # 绘制图形
            self.plot_graph(data, layout, output_name)




##########################################################################################################################################################################################
class mr():
    head = ''
    total_atom = 0
    resid = []
    resname= []
    atomname = []
    index = []
    x = []
    y = []
    z = []
    xyz = []
    last = ''
    metals = []
    coordinators = []
    metal1 = 0
    metal2 = 0
    metal3 = 0
    metal4 = 0
    metal5 = 0
    metal6 = 0
    metal7 = 0
    atom1 = []
    atom2 = []
    atom3 = []
    atom4 = []
    atom5 = []
    atom6 = []
    atom7 = []    
    
    def __init__(self, gro,num_neighbours, distance_value, atom_list, metal_list, residue_list, bond_strength, angle_strength):
        self.GROreader(gro)
        self.MetalMiner(metal_list)
        self.coordinator(num_neighbours, distance_value, atom_list, metal_list, residue_list)
        # self.bond_cal(atom6,bond_strength)
        # self.pair_cal()
        # self.angle_cal(angle_strength)
        self.bond_cal( self.metal1, self.atom1, self.metal2, self.atom2, self.metal3, self.atom3, self.metal4, self.atom4, self.metal5, self.atom5, self.metal6, self.atom6, bond_strength)
        self.pair_cal( self.metal1, self.atom1, self.metal2, self.atom2, self.metal3, self.atom3, self.metal4, self.atom4, self.metal5, self.atom5, self.metal6, self.atom6)
        self.angle_cal( self.metal1, self.atom1, self.metal2, self.atom2, self.metal3, self.atom3, self.metal4, self.atom4, self.metal5, self.atom5, self.metal6, self.atom6, angle_strength)

    def GROreader(self,gro): 
        with open(gro, 'r') as file:
            lines = file.readlines()
            
        # extra lines
        self.head = lines[0].strip()
        self.total_atom = int(lines[1])
        self.last = lines[-1]
        
        # 忽略前两行和最后一行
        lines = lines[2:-1]
        
        # 逐行解析内容
        for line in lines:
            line = line.strip()  # 去除首尾空格和换行符
            match = re.match(r'(\d+)([A-Za-z]{2,})', line)
            if match:
                self.resid.append(int(match.group(1)))
                self.resname.append(str(match.group(2)))
            self.atomname.append(str(line.split()[1]))                        # The 3rd column is the atom name C CA CD1 CD2 and so on
            self.index.append(int(line.split()[2]))                   # Column 4 is the residue name TYR ALA etc.
            self.x.append(float(line.split()[3]))                         # The 5th column is the name of the chain it is on
            self.y.append(float(line.split()[4]))               # The sixth column is the residue number
            self.z.append(float(line.split()[5]))                   # Column 7 is the x-coordinate of the atom
            self.xyz.append([float(line.split()[3]),float(line.split()[4]),float(line.split()[5])])
    
    def MetalMiner(self, metal_list):
        print(metal_list)
        for i in range(len(self.atomname)):
            if self.atomname[i] in metal_list and self.resname[i] in metal_list:
                self.metals.append(self.index[i])
        # the index in list should -1   
        for i in range(len(self.metals)):
            try:
                setattr(self, f'metal{i+1}', self.metals[i] - 1)
            except IndexError:
                pass

        metals_name = [self.atomname[i-1] for i in self.metals]
        sentence = "The metals atom index are: {}".format(list(zip(metals_name, self.metals)))
        print(sentence)
        # print(self.metal1,self.metal2,self.metal3)
            
        
        
        
    def coordinator(self, num_neighbours, distance_value, atom_list, metal_list, residue_list):
        # find the atom index
        if self.metal1 != 0:
            for i in range(len(self.index)):
                if self.resname[i] in residue_list and self.atomname[i] in atom_list:
                    if self.distance(self.metal1, i) <= distance_value:
                        self.atom1.append(i)
            # find the neighbour atoms which are from the same residue.
            index_to_remove = []
            if len(self.atom1) > num_neighbours:
                for i in range(len(self.atom1)-1):
                    if abs(self.atom1[i] - self.atom1[i+1]) == 1:
                        if self.distance(self.metal1, self.atom1[i]) < self.distance(self.metal1, self.atom1[i+1]):
                            index_to_remove.append(i+1)
                        else:
                            index_to_remove.append(i)
                try:
                    sorted_indexes_to_remove = sorted(index_to_remove, reverse=True)
                    for index in sorted_indexes_to_remove:
                        if len(self.atom1) > num_neighbours:
                            del self.atom1[index]
                except:
                    pass
                
                
        if self.metal2 != 0:
            for i in range(len(self.index)):
                if self.resname[i] in residue_list and self.atomname[i] in atom_list:
                    if self.distance(self.metal2, i) <= distance_value:
                        self.atom2.append(i)
            # find the neighbour atoms which are from the same residue.
            
            index_to_remove = []
            if len(self.atom2) > num_neighbours:
                for i in range(len(self.atom2)-1):
                    if abs(self.atom2[i] - self.atom2[i+1]) == 1:
                        if self.distance(self.metal2, self.atom2[i]) < self.distance(self.metal2, self.atom2[i+1]):
                            index_to_remove.append(i+1)
                        else:
                            index_to_remove.append(i)
                try:
                    sorted_indexes_to_remove = sorted(index_to_remove, reverse=True)
                    for index in sorted_indexes_to_remove:
                        if len(self.atom2) > num_neighbours:
                            del self.atom2[index]
                except:
                    pass
                
        if self.metal3 != 0:
            for i in range(len(self.index)):
                if self.resname[i] in residue_list and self.atomname[i] in atom_list:
                    if self.distance(self.metal3, i) <= distance_value:
                        self.atom3.append(i)
            # find the neighbour atoms which are from the same residue.
            index_to_remove = []
            if len(self.atom3) > num_neighbours:
                for i in range(len(self.atom3)-1):
                    if abs(self.atom3[i] - self.atom3[i+1]) == 1:
                        if self.distance(self.metal3, self.atom3[i]) < self.distance(self.metal3, self.atom3[i+1]):
                            index_to_remove.append(i+1)
                        else:
                            index_to_remove.append(i)
                try:
                    sorted_indexes_to_remove = sorted(index_to_remove, reverse=True)
                    for index in sorted_indexes_to_remove:
                        if len(self.atom3) > num_neighbours:
                            del self.atom3[index]
                except:
                    pass

        if self.metal4 != 0:
            for i in range(len(self.index)):
                if self.resname[i] in residue_list and self.atomname[i] in atom_list:
                    if self.distance(self.metal4, i) <= distance_value:
                        self.atom4.append(i)
            # find the neighbour atoms which are from the same residue.
            index_to_remove = []
            if len(self.atom4) > num_neighbours:
                for i in range(len(self.atom4)-1):
                    if abs(self.atom4[i] - self.atom4[i+1]) == 1:
                        if self.distance(self.metal4, self.atom4[i]) < self.distance(self.metal4, self.atom4[i+1]):
                            index_to_remove.append(i+1)
                        else:
                            index_to_remove.append(i)
                try:
                    sorted_indexes_to_remove = sorted(index_to_remove, reverse=True)
                    for index in sorted_indexes_to_remove:
                        if len(self.atom4) > num_neighbours:
                            del self.atom4[index]
                except:
                    pass

        if self.metal5 != 0:
            for i in range(len(self.index)):
                if self.resname[i] in residue_list and self.atomname[i] in atom_list:
                    if self.distance(self.metal5, i) <= distance_value:
                        self.atom5.append(i)
            # find the neighbour atoms which are from the same residue.
            index_to_remove = []
            if len(self.atom5) > num_neighbours:
                for i in range(len(self.atom5)-1):
                    if abs(self.atom5[i] - self.atom5[i+1]) == 1:
                        if self.distance(self.metal5, self.atom5[i]) < self.distance(self.metal5, self.atom5[i+1]):
                            index_to_remove.append(i+1)
                        else:
                            index_to_remove.append(i)
                try:
                    sorted_indexes_to_remove = sorted(index_to_remove, reverse=True)
                    for index in sorted_indexes_to_remove:
                        if len(self.atom5) > num_neighbours:
                            del self.atom5[index]
                except:
                    pass

        if self.metal6 != 0:
            for i in range(len(self.index)):
                if self.resname[i] in residue_list and self.atomname[i] in atom_list:
                    if self.distance(self.metal6, i) <= distance_value:
                        self.atom6.append(i)
            # find the neighbour atoms which are from the same residue.
            index_to_remove = []
            if len(self.atom6) > num_neighbours:
                for i in range(len(self.atom6)-1):
                    if abs(self.atom6[i] - self.atom6[i+1]) == 1:
                        if self.distance(self.metal6, self.atom6[i]) < self.distance(self.metal6, self.atom6[i+1]):
                            index_to_remove.append(i+1)
                        else:
                            index_to_remove.append(i)
                try:
                    sorted_indexes_to_remove = sorted(index_to_remove, reverse=True)
                    for index in sorted_indexes_to_remove:
                        if len(self.atom6) > num_neighbours:
                            del self.atom6[index]
                except:
                    pass

        if self.metal7 != 0:
            for i in range(len(self.index)):
                if self.resname[i] in residue_list and self.atomname[i] in atom_list:
                    if self.distance(self.metal7, i) <= distance_value:
                        self.atom7.append(i)                        
            # find the neighbour atoms which are from the same residue.
            index_to_remove = []
            if len(self.atom7) > num_neighbours:
                for i in range(len(self.atom7)-1):
                    if abs(self.atom7[i] - self.atom7[i+1]) == 1:
                        if self.distance(self.metal7, self.atom7[i]) < self.distance(self.metal7, self.atom7[i+1]):
                            index_to_remove.append(i+1)
                        else:
                            index_to_remove.append(i)
                try:
                    sorted_indexes_to_remove = sorted(index_to_remove, reverse=True)
                    for index in sorted_indexes_to_remove:
                        if len(self.atom7) > num_neighbours:
                            del self.atom7[index]
                except:
                    pass

        
    def distance(self, index1, index2):
        distance = math.sqrt((self.x[index2] - self.x[index1])**2 + (self.y[index2] - self.y[index1])**2 + (self.z[index2] - self.z[index1])**2)
        return distance
    
    def calculate_distance(self, point1, point2):
        x1, y1, z1 = point1
        x2, y2, z2 = point2
        distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
        return distance
    
    def calculate_angle(self, point1, point2, point3):
        vector_ab = np.array(point2) - np.array(point1)
        vector_bc = np.array(point2) - np.array(point3)
    
        dot_product = np.dot(vector_ab, vector_bc)
        norm_ab = np.linalg.norm(vector_ab)
        norm_bc = np.linalg.norm(vector_bc)
    
        cos_angle = dot_product / (norm_ab * norm_bc)
        angle_rad = np.arccos(cos_angle)
        angle_deg = np.degrees(angle_rad)
    
        return angle_deg

    def bond_cal(self,  metal1, atom1, metal2, atom2, metal3, atom3, metal4, atom4, metal5, atom5, metal6, atom6, bond_strength):
        for i in [1,2,3,4,5,6]:
            if locals()["metal" + str(i)] !=0:
                metal = locals()["metal" + str(i)]+1
                atom = [i + 1 for i in locals()["atom" + str(i)]]
                metal_point = [self.x[metal-1],self.y[metal-1],self.z[metal-1]]
                print("; please add below to topol.top's distance part")
                for i in atom:
                    locals()["atom" + str(i)] = [self.x[i-1],self.y[i-1],self.z[i-1]]
                    print("%5d%6d%6d%7.2f%9d" % (metal, i, 6, self.calculate_distance(metal_point, locals()["atom" + str(i)]), bond_strength))
        
    def pair_cal(self,  metal1, atom1, metal2, atom2, metal3, atom3, metal4, atom4, metal5, atom5, metal6, atom6):
        for i in [1,2,3,4,5,6]:
            if locals()["metal" + str(i)] !=0:
                metal = locals()["metal" + str(i)]+1
                atom = [i + 1 for i in locals()["atom" + str(i)]]
                target_resid = [self.resid[x-1] for x in atom]
        
                print("; I added the pairs - so the zn will not nonbonded interact with the cyx residues")
                for i in range(len(self.atomname)):
                    if self.atomname[i] == 'CA' and self.resid[i] in target_resid:
                        print("%5d%6d%6d" % (metal,i+1,1))
           
    
    def angle_cal(self,  metal1, atom1, metal2, atom2, metal3, atom3, metal4, atom4, metal5, atom5, metal6, atom6, angle_strength):
        for i in [1,2,3,4,5,6]:
            if locals()["metal" + str(i)] !=0:
                metal = locals()["metal" + str(i)]+1
                atom = [i + 1 for i in locals()["atom" + str(i)]]
                # define how many neighbour atoms
                neighbour = len(atom)
                if neighbour >= 2:
                    metal_point = [self.x[metal-1],self.y[metal-1],self.z[metal-1]]
                    for i in atom:
                        locals()["atom" + str(i)] = [self.x[i-1],self.y[i-1],self.z[i-1]]
                    print("[ angle_restraints ]")
                    if neighbour == 2:
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[0], metal, atom[1],1,self.calculate_angle(locals()["atom" + str(atom[0])],metal_point, locals()["atom" + str(atom[1])]),angle_strength, 1))
                    elif neighbour == 3:
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[0], metal, atom[1],1,self.calculate_angle(locals()["atom" + str(atom[0])],metal_point, locals()["atom" + str(atom[1])]),angle_strength, 1))
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[0], metal, atom[2],1,self.calculate_angle(locals()["atom" + str(atom[0])],metal_point, locals()["atom" + str(atom[2])]),angle_strength, 1))
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[1], metal, atom[2],1,self.calculate_angle(locals()["atom" + str(atom[1])],metal_point, locals()["atom" + str(atom[2])]),angle_strength, 1))
                    elif neighbour == 4:
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[0], metal, atom[1],1,self.calculate_angle(locals()["atom" + str(atom[0])],metal_point, locals()["atom" + str(atom[1])]),angle_strength, 1))
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[1], metal, atom[2],1,self.calculate_angle(locals()["atom" + str(atom[1])],metal_point, locals()["atom" + str(atom[2])]),angle_strength, 1))
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[2], metal, atom[3],1,self.calculate_angle(locals()["atom" + str(atom[2])],metal_point, locals()["atom" + str(atom[3])]),angle_strength, 1))
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[0], metal, atom[3],1,self.calculate_angle(locals()["atom" + str(atom[0])],metal_point, locals()["atom" + str(atom[3])]),angle_strength, 1))
            
            
    def GROwriter(self, gro):
        print(self.head)
        print("%5d"  % (self.total_atom))
        for i in range(len(self.resid)):    
            print("%5d%-3s%7s%5d%8.3f%8.3f%8.3f" %  (self.resid[i], self.resname[i], self.atom[i], self.index[i], self.x[i], self.y[i], self.z[i]))
        print(self.last)

##########################################################################################################################################################################################

class gromerger():
    def __init__(self, receptor_gro, ligand_gro, ligand_itp, receptor_top):
        self.merge_gro_files(receptor_gro, ligand_gro, ligand_itp)
        
    def merge_gro_files(self, receptor_gro, ligand_gro, ligand_itp):
        # Merge the two gro files
        with open(ligand_gro, 'r') as ligand_file, open(receptor_gro, 'r') as receptor_file:
            ligand_lines = ligand_file.readlines()[1:-1]
            receptor_lines = receptor_file.readlines()
            a = int(receptor_lines[1].split()[0])
            b = int(ligand_lines[0].split()[0])
            c = a + b
            receptor_lines[1] = f"{c}\n"
            
            with open('complex.gro', 'w') as complex_file:
                complex_file.writelines(receptor_lines[0:-1])
                complex_file.writelines(ligand_lines[1:])
                complex_file.writelines(receptor_lines[-1])
        
        # Edit the topol.top file         
        with open('include.dat', 'w') as include_file:
            include_file.write("; Include ligand topology\n")
            include_file.write(f"#include \"{ligand_itp}\"\n")
            include_file.write("#ifdef POSRES_LIG\n")
            include_file.write("#include \"posre_lig.itp\"\n")
            include_file.write("#endif\n")
        
        # adding content to the end of topol file
        with open(ligand_itp, 'r') as ligand_itp_file:
            # find ligand molecule name
            nrexcl_line = None
            lines = ligand_itp_file.readlines()
            for i in range(len(lines)):
                if 'nrexcl' in lines[i]:
                    nrexcl_line = lines[i+1].strip()
                    break
            
            if nrexcl_line:
                B = nrexcl_line.split()[0]
                with open('topol.top', 'r') as topol_file:
                    topol_lines = topol_file.readlines()
                # add lig_GMX.itp record to topol.top
                for i, line in enumerate(topol_lines):
                    if 'moleculetype' in line:
                        topol_lines.insert(i, f"#endif\n")
                        topol_lines.insert(i, f"#include \"posre_lig.itp\"\n")
                        topol_lines.insert(i, f"#ifdef POSRES_LIG\n")
                        topol_lines.insert(i, f"#include \"{ligand_itp}\"\n")
                        topol_lines.insert(i, f"; Include ligand topology\n")
                        break
                # add ligand molecule type to topol.top
                for i, line in enumerate(topol_lines):
                    if 'molecules' in line:
                        topol_lines.insert(-1, f"{B}                1\n")
                        break
                
                with open('topol.top', 'w') as topol_file:
                    topol_file.writelines(topol_lines)
    
    
##########################################################################################################################################################################################
class contact_map_detect(): # read uploaded files
    protein = ''
    ligand = ''
    def __init__(self, topol, traj, lig, output, distance):
        contact_map = self.calculate_contact(topol, traj, lig, output, distance)
        self.plot(contact_map, output_name, distance)
        self.csv_writer(contact_map, output_name)
        
    def calculate_contact(self, topol, traj, lig, output, distance):
        # 加载蛋白质和配体的拓扑和轨迹文件
        u = mda.Universe(topol, traj)
        
        # 选择蛋白质和配体
        self.protein = u.select_atoms('protein')
        self.ligand = u.select_atoms('resname ' + lig)  
        
        # 初始化接触图矩阵
        contact_map = np.zeros((len(u.trajectory), len(self.protein.residues)))
    
        # 计算每帧的接触情况
        for ts in u.trajectory:
            y_ticks = [] 
            y_labels= []
            frame_index = ts.frame
            for i, residue in enumerate(self.protein.residues):
                min_dist = np.min(contacts.distance_array(residue.atoms.positions, self.ligand.positions))
                y_ticks.append(i)
                y_labels.append(f'{residue.resid} {residue.resname}')
                if min_dist < distance:
                    contact_map[frame_index, i] = 1
        
        return contact_map
    
    def plot(self, contact_map, output_name, distance):
        resid_list = [i for i in range(self.protein.residues.resids[0], len(self.protein.residues.resids)+self.protein.residues.resids[0], 5)]
        resname_list =[]
        for i, residue in enumerate(self.protein.residues):
            if residue.resid in resid_list:
                resname_list.append(residue.resname)
        i_list = range(len(resid_list))
        # print(resid_list)
        # print(resname_list)
        plt.imshow(contact_map.T, aspect='auto', origin='lower', cmap='Greys')
        plt.xlabel('Time (ns)')
        plt.yticks(ticks=resid_list, labels=['' if (resid_list[i]-self.protein.residues.resids[0]) % 20 != 0 else f'{resid_list[i]} {resname_list[i]}' for i in i_list])
        plt.ylabel('Residue Index')
        plt.title('Protein-Ligand Contact Map')
        plt.suptitle('Distance < ' + str(distance))
        plt.colorbar(label='Contact', ticks=[0, 1])
        # plt.show()
        plt.savefig(output_name)


    def csv_writer(self, contact_map, output_name):
        # 将contact_map数据写入CSV文件
        with open('ligand_contact.csv', 'w', newline='') as f:
            writer = csv.writer(f)
            # 写入标题行，假设每列代表一个残基，每行代表一个时间帧
        #     header = ['Frame'] + [f'Residue_{i}' for i in range(1, len(protein.residues) + 1)]
            header = ['Time (ns)'] + [f'{residue.resid}_{residue.resname}' for i, residue in enumerate(self.protein.residues)]
            writer.writerow(header)
        
            # 写入每一行的数据
            for frame_index, contacts_ in enumerate(contact_map):
                row = [frame_index] + contacts_.tolist()
                writer.writerow(row)
        
        

##########################################################################################################################################################################################
class pep2lig(): 
    pdb     = ''
    pepname = ''
    chain   = 'A'
    resnum  = 1
    atomic_index        = []
    atomic_name         = []
    residue_name        = []
    chain_name          = []
    residue_index       = []
    X_peratom           = []
    Y_peratom           = []
    Z_peratom           = []
    bfactor_per_factor  = []
    temp_factor         = []
    Atomtype_per_atom   = []
    def __init__(self, pdb, pepname):
        self.pdb = pdb
        self.pepname = pepname
        self.converter(pdb, pepname, self.chain, self.resnum)
        self.PDBwriter("modified_" + pdb)
        
    def converter(self, pdb, pepname, chain, resnum):
            with open(pdb, 'r') as infile:
                for line in infile:                                                              # iterate each line in file "f"                     
                        if(line.split()[0] in["ATOM","HETATM"]):       # Judgment Sentence，Used to split each row and then determine whether the first column of the row == ATOM or HETATM
                            self.atomic_index.append(int(line[6:11].strip()))                # The second column is the atomic number
                            self.atomic_name.append(line[12:16].strip())                        # The 3rd column is the atom name C CA CD1 CD2 and so on
                            # self.residue_name.append(line[17:20].strip())                       # Column 4 is the residue name TYR ALA etc.
                            self.residue_name.append(pepname)
                            # self.chain_name.append(line[21].strip())                         # The 5th column is the name of the chain it is on
                            self.chain_name.append(chain)
                            # self.residue_index.append(int(line[22:26].strip()))               # The sixth column is the residue number
                            self.residue_index.append(resnum)
                            self.X_peratom.append(float(line[30:38].strip()))
                            self.Y_peratom.append(float(line[38:46].strip()))
                            self.Z_peratom.append(float(line[46:54].strip()))
                            self.bfactor_per_factor.append(float(line[54:60].strip()) if line[54:60].strip() else 0.0)
                            self.temp_factor.append(float(line[60:66].strip()) if line[60:66].strip() else 0.0 )
                            try:
                                self.Atomtype_per_atom.append(line[76:78].strip())
                            except:
                                self.Atomtype_per_atom.append(" ")              
        
    def PDBwriter(self,filename):
        f = open(filename, "w")                                                             # e.g: f = linesplit[0]+"_PO3.pdb"
        for i in range (0 ,len(self.atomic_index)):                                         # Create a loop, i is a sequence starting from 0, and the number of atoms is the length  
            print("%4s%7d  %-4s%1s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s" %  ("ATOM" ,     # Formatted output, %4s, right-aligned, the output occupies 4 columns in total. If the length is less than 4 columns, the left end will be filled with spaces. If it is greater than 4 columns, the actual length will be output as a string
                                             self.atomic_index[i],                          # %7d, right-aligned, the output occupies a total of 7 columns, if the length is less than 7 columns, the left end is filled with spaces, signed decimal certificate integer
                                             self.atomic_name[i],                           # %-4s, left-aligned, the output occupies a total of 4 columns, if the length is less than 4 columns, the right end is filled with spaces, if it is greater than 4 columns, the actual length is output as a string
                                             self.residue_name[i],                          # %1s, right-aligned, the output occupies a total of 1 column. If it is less than 1 column, it will be filled with spaces from the left end. If it is greater than 1 column, the actual length will be output as a string
                                             "A",                            # %2s, right-aligned, the output occupies 2 columns in total. If it is less than 2 columns, it will be filled with spaces from the left end. If it is greater than 2 columns, the actual length will be output as a string
                                             self.residue_index[i],                         # %4d, right-aligned, the output occupies a total of 4 columns, if the length is less than 4 columns, the left end is filled with spaces, a signed decimal certificate integer
                                             self.X_peratom[i],                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                             self.Y_peratom[i],                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                             self.Z_peratom[i],                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                             self.bfactor_per_factor[i],                    # %6.2f, right-aligned, the output occupies a total of 6 columns, including 2 decimal places, if the width of the value is less than 6, fill in a space at the left end, decimal
                                             self.temp_factor[i],                     # %6.2f, right-aligned, the output occupies a total of 6 columns, including 2 decimal places, if the width of the value is less than 6, fill in a space at the left end, decimal
                                             self.Atomtype_per_atom[i]), file = f )         # %12s, right-aligned, the output occupies a total of 12 columns, if it is less than 12 columns, it will be filled with spaces from the left end
        print("END", file = f)
        f.close()
##########################################################################################################################################################################################
class gmx_dssp():
    def __init__(self, ds_data, ds_traj, ds_output_name, ds_original, ds_color):
        self.plot_figure(ds_data, ds_traj, ds_output_name, ds_original, ds_color)
    
    def read_time_and_residue(self, traj):
        # 载入PDB文件
        u = mda.Universe(traj)
        
        # 提取时间戳和残基列表
        times = [ts.time for ts in u.trajectory]  # 假设时间单位是ns 
        residues = [res.resname + str(res.resid) for res in u.residues]
        return times, residues
 
    def detect_break(self, first_line, residue_list):
        # Find the positions of '=' in the first line
        equal_positions = [pos for pos, char in enumerate(first_line) if char == "="]
        # Sort the positions in reverse order
        equal_positions.sort(reverse=True)
        # Insert 'break' into the residue_list at the corresponding positions from the end
        for pos in equal_positions:
            if pos < len(residue_list):  # Only insert if the position is within the bounds of the residue_list
                residue_list.insert(pos, 'break')
        return residue_list
        
    def read_data(self, data, traj, original, unique_color_bar): 
        times, residues = self.read_time_and_residue(traj)
        # 读取DSSP数据（dssp.dat）
        with open(data, "r") as file:
            dssp_lines = file.readlines()
            first_line = dssp_lines[0]
            residues = self.detect_break(first_line, residues)
    
        # 转换DSSP数据为列表
        dssp_data = [list(line.strip()) for line in dssp_lines]
        
        # 创建DataFrame
        df = pd.DataFrame(data=dssp_data, index=times, columns=residues)
        
        if original == 'false' and unique_color_bar == 'true':
            # times = times.append(times[-1]+1) # for scale the color
            # Calculate the length of each third
            third_length = len(df.columns) // 3
            
            # Create a new row
            new_row = [1] * third_length + [2] * third_length + [3] * (len(df.columns) - 2 * third_length)
            
            # Append the new row at the bottom of the DataFrame
            df.loc[len(df)] = new_row
        elif original != 'false' and unique_color_bar == 'true':
            tenth_length = len(df.columns) // 10
            # Initialize an empty list for the new row
            new_row = []
            # Create new row by appending numbers from 1 to 10, each repeated 'tenth_length' times
            # Note that the last segment fills the remainder of the row if it's not evenly divisible by 10
            for i in range(10):
                if i < 9:
                    new_row += [i] * tenth_length
                else:
                    # The last segment includes any extra columns
                    new_row += [i] * (len(df.columns) - len(new_row))

            # Append the new row at the bottom of the DataFrame
            df.loc[len(df)] = new_row
            pass
        
        return df

    def replace_letters(self, original):
        # 定义转换字典
        if original != 'false':
            structure_values = {
                'H': 9,   # 'alpha-helix',
                'B': 8,   # 'residue in isolated beta-bridge',
                'E': 7,   # 'extended strand that participates in beta-ladder',
                'G': 6,   # '3_10-helix',
                'I': 5,   # 'pi-helix',
                'P': 4,   # 'kappa-helix',  # Assuming kappa-helix is synonymous 
                'S': 3,   # 'bend',
                'T': 2,   # 'hydrogen-bonded turn',
                '=': 1,   # 'break',
                '~': 0    # 'loop'  # No special secondary structure designation
            }
        else:
            structure_values = {
                'H': 2,   # 'alpha-helix',
                'B': 1,   # 'residue in isolated beta-bridge',
                'E': 3,   # 'beta-sheet',
                'G': 1,   # '3_10-helix',
                'I': 1,   # 'pi-helix',
                'P': 1,   # 'kappa-helix',  # Assuming kappa-helix is synonymous 
                'S': 1,   # 'beta-bend',
                'T': 1,   # 'beta-turn',
                '=': 1,   # 'break',
                '~': 1    # 'loop'  # No special secondary structure designation
            }
        return structure_values

    def plot_figure(self, data, traj, outputname, original, color):
        df = self.read_data(data, traj, original,color)
        structure_values = self.replace_letters(original)    
        # 使用字典转换DataFrame中的值
        df.replace(structure_values, inplace=True)


        # Define color scale and color bar settings
        colorscale = [[0, 'red'], [0.5, 'green'], [1.0, 'rgb(0, 0, 255)']]
        colorscale = [[0.00, "gold"],   [0.33, "gold"], [0.33, "mediumturquoise"], [0.66, "mediumturquoise"], [0.66, "lightsalmon"],  [1.00, "lightsalmon"]]
        original_colorscale = [[0.0, "#636EFA"], [0.1, "#636EFA"], [0.1, "#EF553B"], [0.2, "#EF553B"], [0.2, "#00CC96"], [0.3, "#00CC96"], [0.3, "#AB63FA"], [0.4, "#AB63FA"], [0.4, "#FFA15A"], [0.5, "#FFA15A"], [0.5, "#19D3F3"], [0.6, "#19D3F3"], [0.6, "#FF6692"], [0.7, "#FF6692"], [0.7, "#B6E880"], [0.8, "#B6E880"], [0.8, "#FF97FF"], [0.9, "#FF97FF"], [0.9, "#FECB52"], [1.0, "#FECB52"]]
        # colorscale = [[0.00, "red"],   [0.33, "red"], [0.33, "green"], [0.66, "green"], [0.66, "blue"],  [1.00, "blue"]]
        # 将颜色条的刻度设置为固定的值
        if original == 'false':
            colorbar_ticks = [1, 2, 3]  # The actual values in your data
            colorbar_ticktext = ['loop', 'helix', 'beta sheet']  # Descriptive labels
        elif original != 'false':
            colorbar_ticks = [0,1,2,3,4,5,6,7,8,9]
            colorbar_ticktext = ['loop', 'break', 'h-bond turn', 'bend', 'kappa helix', 'pi helix', '3_10 helix', 'strand', 'beta bridge', 'a-helix']
        if original != 'false':
            # 创建热图
            fig = go.Figure(data=go.Heatmap(
                z=df.T.values,
                x=df.index,
                y=df.columns,
                colorscale=original_colorscale,
                colorbar=dict(tickmode = 'array', tickvals=colorbar_ticks, ticktext=colorbar_ticktext),
                hoverongaps=False))
        elif original == 'false':
            # 创建热图
            fig = go.Figure(data=go.Heatmap(
                z=df.T.values,
                x=df.index,
                y=df.columns,
                colorscale=colorscale,
                colorbar=dict(tickmode = 'array', tickvals=colorbar_ticks, ticktext=colorbar_ticktext),
                hoverongaps=False))
           
        fig.update_layout(
            title='Secondary Structure Analysis Over Time',
            yaxis_title='Residue',
            xaxis_title='Time (ns)',
            width=800,
            height=600)

        # 显示热图
        pio.write_image(fig, outputname)

class renumber_MODEL():
    def __init__(self, files):
        # 初始化计数器
        count = 1
        
        # 打开输入文件和输出文件
        with open(files, 'r') as file, open(f'renumbered_{files}', 'w') as output:
            # 遍历文件中的每一行
            for line in file:
                # 使用正则表达式检查行是否包含"MODEL"和后面的数字
                if re.match(r"^\s*MODEL\s+\d+", line):
                    # 替换匹配的行为"MODEL"后面接计数器的值，并将计数器加一
                    new_line = re.sub(r"^\s*MODEL\s+\d+", f"MODEL        {count}", line)
                    count += 1
                    output.write(new_line)
                else:
                    output.write(line)


##########################################################################################################################################################################################
class ff_res_adder():

    def __init__(self, lig_itp, lig_gro, ff_rtp, ff_hdb, ff_bonded, ff_nonbonded, ff_aomtypes, ff_restype, atom_types, isamino_acid, output_name):
        self.add_res_to_ff(lig_itp, lig_gro, ff_rtp, ff_hdb, ff_bonded, ff_nonbonded, ff_aomtypes, ff_restype, atom_types, isamino_acid, output_name)

    def add_res_to_ff(self, lig_itp, lig_gro, ff_rtp, ff_hdb, ff_bonded, ff_nonbonded, ff_aomtypes, ff_restype, atom_types, isamino_acid, output_name):
        # parse the lig_GMX.itp file
        itp_content_dicts    =   self.parse_itp(lig_itp, atom_types)
        # parse the lig_GMX.gro file
        residue_numbers, residue_names, atom_names, atom_numbers, positions = self.parse_gro(lig_gro)
        # replace the digit number to atom_types for bonds, angles, dihedrals, impropers
        bond_records_dict   =   self.generate_bond_records(itp_content_dicts)
        # the notes for residuetypes.dat
        ff_restype_dict     =   self.parse_residuetypes(itp_content_dicts)
        # the content to be added to aminoacids.rtp
        ff_rtp_dict         =   self.parse_rtp(itp_content_dicts, isamino_acid) 
        # the content to be added to aminoacids.hdb          
        ff_hdb_dict         =   self.parse_hdb(itp_content_dicts) 
        # parse the ffbonded.itp file
        ff_bonded_dict      =   self.parse_bonded(ff_bonded) 
        # find the missed bonds records
        missing_bonds_dict  =   self.find_missing_bonds(bond_records_dict, ff_bonded_dict)
        # 打印嵌套字典
        # self.print_nested_dict(bond_records_dict) 
        # self.print_nested_dict(ff_bonded_dict) 
        # Get the list of atom types
        # ff_nonbonded_dict   =   self.parse_nonbonded(ff_nonbonded)          
        ff_aomtypes_dict    =   self.parse_atomtypes(ff_aomtypes, ff_nonbonded) 
        # self.print_nested_dict(ff_aomtypes_dict)  
        # Compare ff_atomtype and user provided atom_types
        missing_atom_type   =   self.find_missing_atomtypes(ff_aomtypes_dict, atom_types)
        # bond_values_dict    =   self.calculate_bond_values(positions, itp_content_dicts)

     
        
         
    def parse_itp(self, lig_itp, atom_types):
        header_atomtypes = ['name', 'bond_type', 'mass', 'charge', 'ptype', 'sigma', 'epsilon', 'Amb_sigma', 'Amb_epsilon']      
        header_moleculetype = ['name', 'nrexcl']        
        header_atoms = ['atom_number', 'atom_type', 'resid', 'resname', 'atom_name', 'cgnr', 'partial_charge', 'mass', 'qtot', 'total_charge']        
        header_bonds = ['a1', 'a2', 'funct', 'bond_length', 'bond_strength', 'a1_name', 'a2_name']        
        header_pairs = ['a1', 'a2', 'funct', 'a1_name', 'a2_name']        
        header_angles = ['a1', 'a2', 'a3', 'funct', 'angle_degree', 'angle_strength', 'a1_name', 'a2_name', 'a3_name']        
        header_dihedrals = ['a1', 'a2', 'a3', 'a4', 'funct', 'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'a1_name', 'a2_name', 'a3_name', 'a4_name']        
        header_impropers = ['a1', 'a2', 'a3', 'a4', 'funct', 'angle_degree', 'angle_strength', 'multiplicity', 'a1_name', 'a2_name', 'a3_name', 'a4_name']    
        header_list = [header_atomtypes, header_moleculetype,  header_atoms ,  header_bonds ,  header_pairs ,  header_angles ,  header_dihedrals ,  header_impropers ]
    
        # Define dictionaries to store different sections
        sections = {
            'atomtypes': {},
            'moleculetype': {},
            'atoms': {},
            'bonds': {},
            'pairs': {},
            'angles': {},
            'dihedrals': {},
            'impropers': {}
        }
    
        # Open and read the file
        with open(lig_itp, 'r') as file:
            lines = file.readlines()

        # Helper variables to track sections and headers
        current_section = None
        headers = []
        data_rows = []
        count = 0

        for line in lines:
            line = line.strip()
            if line.startswith('['):  # Detect new section
                if current_section and headers:
                    # save the datas to target section
                    transposed_data = list(zip(*data_rows))
                    sections[current_section] = {key: tuple(values) for key, values in zip(headers, transposed_data)}          
                # Extract section name from format [ section_name ] and reset the previous section
                section_name = line.split()[1]
                current_section = section_name
                # headers = []  # Reset headers for the new section
                headers = header_list[count]
                count += 1 
                data_rows = []
            # define whether it's improper dihedral part
            elif current_section == 'dihedrals' and line.startswith(';') and 'propers' in line:
                current_section = 'impropers'
            # read contents
            elif line and not line.startswith(';'):  # Valid data lines
                if current_section and headers:  # We already have headers, parse data
                    a = line.split()
                    # remove the values = ; or -
                    a = [item for item in a if item not in [';', '-']]
                    # remove the "-" for dihedrals or improper dihedrals
                    if current_section in ['dihedrals', 'impropers']:
                        a = a[:-4] + [item.replace('-', '') for item in a[-4:]]
                    # a = [item for item in a if item != ';']
                    data_rows.append(a)
        # 最后一个节的数据保存
        if current_section and headers:
            transposed_data = list(zip(*data_rows))
            sections[current_section] = {key: tuple(values) for key, values in zip(headers, transposed_data)}
        # change the atom_type to user's input
        if len(atom_types) >= len(sections['atoms']['atom_type']):
            sections['atoms']['atom_type'] = tuple(atom_types)
        else:
            print("Error: atom_types list does not have enough elements to replace all atom_types.")

        return sections
    
        
    def parse_gro(self, lig_gro):
        
        
        # 定义空列表以存储每列的数据
        residue_numbers = []
        residue_names = []
        atom_names = []
        atom_numbers = []
        positions    = []
        
        # 读取文件并提取所需行
        with open(lig_gro, 'r') as file:
            lines = file.readlines()
            # 获取第3行到倒数第二行的数据行
            data_lines = lines[2:-1]
        
        # 解析每行数据
        for line in data_lines:
            # 假设每个元素之间由多个空格分隔，这是一个典型的固定列宽格式
            parts = line.split()
            if len(parts) >= 7:  # 确保每行数据都是完整的
                residue_numbers.append(parts[0])
                residue_names.append(parts[1])
                atom_names.append(parts[2])
                atom_numbers.append(parts[3])
                positions.append([float(parts[4]), float(parts[5]), float(parts[6])])


        return residue_numbers, residue_names, atom_names, atom_numbers, positions 
    def parse_rtp(self, itp_content_dict, isamino_acids):  
        print("######################## Please add below to aminoacids.rtp! ########################")
        print(f"[ {itp_content_dict['atoms']['resname'][0]} ]")
        print(" [ atoms ]")
        for atomname, atomtype, charge, indices in zip(itp_content_dict['atoms']['atom_name'], itp_content_dict['atoms']['atom_type'], itp_content_dict['atoms']['partial_charge'], itp_content_dict['atoms']['atom_number']):
            print("%6s    %-12s%8.5f%6d" % (atomname, atomtype, float(charge), int(indices)))
        print(" [ bonds ]")
        for atom1, atom2 in zip(itp_content_dict['bonds']['a1_name'], itp_content_dict['bonds']['a2_name']):
            print("%6s%6s" % (atom1,atom2))
        if isamino_acids == "true":
            print("%6s%6s" % ("-C","N"))
        print(" [ impropers ]")
        if isamino_acids == "true":
            print("%6s%6s%6s%6s" % ("-C","CA","N","H"))
            print("%6s%6s%6s%6s" % ("CA","+N","C", "O"))
        for atom1,atom2,atom3,atom4 in zip(itp_content_dict['impropers']['a1_name'], itp_content_dict['impropers']['a2_name'], itp_content_dict['impropers']['a3_name'], itp_content_dict['impropers']['a4_name']):
            print("%6s%6s%6s%6s" % (atom1,atom2,atom3,atom4))
        


         
    def parse_hdb(self, itp_content_dict):  
        hdb_bond_types = """
1. one planar hydrogen, e.g. rings or peptide bond
2. one single hydrogen, e.g. hydroxyl
3. two planar hydrogens, e.g. ethylene -C=CH2, or amide -C(=O)NH2
4. two or three tetrahedral hydrogens, e.g. -CH3
5. one tetrahedral hydrogen, e.g. C3* CH*
6. two tetrahedral hydrogens, e.g. C-CH2 *-C*
"""
        def find_connections(atom, bonds):
            connections = {'H': [], 'non_H': []}
            for a1, a2 in zip(bonds['a1_name'], bonds['a2_name']):
                if a1 == atom:
                    (connections['H'] if 'H' in a2 else connections['non_H']).append(a2)
                if a2 == atom:
                    (connections['H'] if 'H' in a1 else connections['non_H']).append(a1)
            return connections
        def define_bond_type(atoms, counts):
            bond_type = "5"
            if atoms == "N":
                bond_type = str(1)
            elif atoms == "CA" and counts == 2:
                bond_type = str(6)
            elif atoms == "CA" and counts == 1:
                bond_type = str(5)
            elif atoms == "CB" and counts == 2:
                bond_type = str(6)
            elif atoms == "CB" and counts == 1:
                bond_type = str(5)
            elif counts == 1:
                bond_type = str(5)
            elif counts == 2:
                bond_type = str(6)
            elif counts == 3:
                bond_type = str(4)
            return bond_type
        # Main processing
        output = []
        for atom in itp_content_dict['atoms']['atom_name']:
            if 'H' not in atom:
                conns = find_connections(atom, itp_content_dict['bonds'])
                if conns['H']:
                    # Only consider atoms directly bonded to H
                    count_hydrogens = len(conns['H'])
                    bond_type = define_bond_type(atom, count_hydrogens)
                    hydrogens = "H" + str(atom[1:] if len(atom) > 1 else "")
                    others = '       '.join(conns['non_H'])
                    # output.append(f"{count_hydrogens}\t6\t{hydrogens}\t{atom}\t{', '.join(conns['non_H'])}")
                    output.append(f"{count_hydrogens}\t{bond_type}\t{hydrogens}\t{atom}\t{others}")
        # Print output as per the format discussed
        print("######################## Please add below to aminoacids.hdb! ########################")
        print("CR1\t" + str(len(output)))
        for line in output:
            print(line)
        print("\n","##Rules of the Hydrogen adding!##")        
        print(hdb_bond_types)
        
               
    def parse_residuetypes(self, itp_content_dict):  
        print("######################## Please add below to residuetypes.dat! ########################")
        print("%-8s%-7s" % (itp_content_dict['atoms']['resname'][0],"Protein"))

    def generate_bond_records(self, itp_content_dict):
        # 建立atom_number到atom_type的映射
        atom_type_mapping = {num: typ for num, typ in zip(itp_content_dict['atoms']['atom_number'], itp_content_dict['atoms']['atom_type'])}
        
        # 替换bonds字典中的a1和a2为对应的atom_type
        updated_a1 = [atom_type_mapping[num] for num in itp_content_dict['bonds']['a1']]
        updated_a2 = [atom_type_mapping[num] for num in itp_content_dict['bonds']['a2']]
        
        # 更新bonds字典
        itp_content_dict['bonds']['a1'] = tuple(updated_a1)
        itp_content_dict['bonds']['a2'] = tuple(updated_a2)

        # 替换angles字典中的a1, a2, a3
        updated_a1 = [atom_type_mapping[num] for num in itp_content_dict['angles']['a1']]
        updated_a2 = [atom_type_mapping[num] for num in itp_content_dict['angles']['a2']]
        updated_a3 = [atom_type_mapping[num] for num in itp_content_dict['angles']['a3']]

        # 更新angles字典
        itp_content_dict['angles']['a1'] = tuple(updated_a1)
        itp_content_dict['angles']['a2'] = tuple(updated_a2)
        itp_content_dict['angles']['a3'] = tuple(updated_a3)

        # 替换dihedrals字典中的a1, a2, a3, a4
        updated_a1 = [atom_type_mapping[num] for num in itp_content_dict['dihedrals']['a1']]
        updated_a2 = [atom_type_mapping[num] for num in itp_content_dict['dihedrals']['a2']]
        updated_a3 = [atom_type_mapping[num] for num in itp_content_dict['dihedrals']['a3']]
        updated_a4 = [atom_type_mapping[num] for num in itp_content_dict['dihedrals']['a4']]

        # 更新dihedrals字典
        itp_content_dict['dihedrals']['a1'] = tuple(updated_a1)
        itp_content_dict['dihedrals']['a2'] = tuple(updated_a2)
        itp_content_dict['dihedrals']['a3'] = tuple(updated_a3)
        itp_content_dict['dihedrals']['a4'] = tuple(updated_a4)

        # 替换impropers字典中的a1, a2, a3, a4
        updated_a1 = [atom_type_mapping[num] for num in itp_content_dict['impropers']['a1']]
        updated_a2 = [atom_type_mapping[num] for num in itp_content_dict['impropers']['a2']]
        updated_a3 = [atom_type_mapping[num] for num in itp_content_dict['impropers']['a3']]
        updated_a4 = [atom_type_mapping[num] for num in itp_content_dict['impropers']['a4']]

        # 更新impropers字典
        itp_content_dict['impropers']['a1'] = tuple(updated_a1)
        itp_content_dict['impropers']['a2'] = tuple(updated_a2)
        itp_content_dict['impropers']['a3'] = tuple(updated_a3)
        itp_content_dict['impropers']['a4'] = tuple(updated_a4)

        return itp_content_dict

    def parse_bonded(self, ff_bonded): 
        header_bondtypes = ['a1', 'a2', 'funct', 'bond_length', 'bond_strength']      
        header_constrainttypes = ['a1', 'a2', 'funct', 'bond_length']        
        header_angletypes = ['a1', 'a2', 'a3', 'funct', 'angle_degree', 'angle_strength']        
        header_dihedral_f1 = ['a1', 'a2', 'a3', 'a4', 'funct', 'angle_degree', 'angle_strength', 'multiplicity']        # acpype for non-periodic improper dihedrals, to fix the dihedral
        header_dihedral_f3 = ['a1', 'a2', 'a3', 'a4', 'funct', 'c0', 'c1', 'c2', 'c3', 'c4', 'c5']                      # acpype for proper dihedrals
        header_dihedral_f4 = ['a1', 'a2', 'a3', 'a4', 'funct', 'angle_degree', 'angle_strength', 'multiplicity']        # Gromacs for periodic improper dihedral, have some flexibility to spin
        header_dihedral_f9 = ['a1', 'a2', 'a3', 'a4', 'funct', 'angle_degree', 'angle_strength', 'multiplicity']        # Gromacs for proper dihedral (multiple), multiple is the number of energy peaks or valleys
        header_list = [header_bondtypes, header_constrainttypes,  header_angletypes ,  header_dihedral_f1 ,  header_dihedral_f3 ,  header_dihedral_f4 ,  header_dihedral_f9]
    
        # Define dictionaries to store different sections
        sections = {
            'bondtypes': {},
            'constrainttypes': {},
            'angletypes': {},
            'dihedral_f1': {},
            'dihedral_f3': {},
            'dihedral_f4': {},
            'dihedral_f9': {}
        }
    
        # Open and read the file
        with open(ff_bonded, 'r') as file:
            lines = file.readlines()

        # Helper variables to track sections and headers
        current_section = None
        headers = []
        data_rows = []
        count = 0

        for line in lines:
            line = line.strip()
            if line.startswith('['):  # Detect new section
                # step-5. write the data to dictionary!
                if current_section and headers:
                    # save the datas to target section
                    transposed_data = list(zip(*data_rows))
                    sections[current_section] = {key: tuple(values) for key, values in zip(headers, transposed_data)}          
                # step-1. Extract section name from format [ section_name ] and reset the previous section
                section_name = line.split()[1]
                current_section = section_name
                # step-2. Assign headers for the columns
                if current_section == "bondtypes":
                    count = 0
                elif current_section == "constrainttypes":
                    count = 1
                elif current_section == "angletypes":
                    count = 2
                elif current_section == "dihedraltypes":
                    count = 4
                elif current_section == "dihedraltypes":
                    count = 4
                elif current_section == "dihedraltypes":
                    count = 4
                elif current_section == "dihedraltypes":
                    count = 4                    
                headers = header_list[count]

                data_rows = []
            # step-3. Define the type of dihedrals, also write the 1st row to data_rows
            elif current_section == 'dihedraltypes' and not line.startswith((';', '[', '#')):
                a = line.split()
                # remove the values = ; or -
                a = [item for item in a if item not in [';', '-']]
                data_rows.append(a)                
                if str(line.split()[4]) == '1':
                    current_section = 'dihedral_f1'
                    count = 3
                if str(line.split()[4]) == '3':
                    current_section = 'dihedral_f3'
                    count = 4
                if str(line.split()[4]) == '4':
                    current_section = 'dihedral_f4'
                    count = 5
                if str(line.split()[4]) == '9':
                    current_section = 'dihedral_f9' 
                    count = 6
                headers = header_list[count]
            # step-4. Read contents
            elif line and not line.startswith((';', '[', '#')):  # Valid data lines
                if current_section and headers:  # We already have headers, parse data
                    a = line.split()
                    # remove the values = ; or -
                    a = [item for item in a if item not in [';', '-']]
                    data_rows.append(a)
        # step-6. 最后一个节的数据保存
        if current_section and headers:
            transposed_data = list(zip(*data_rows))
            sections[current_section] = {key: tuple(values) for key, values in zip(headers, transposed_data)}

        return sections
    
    def print_nested_dict(self, dictionary, indent=0):
        for key, value in dictionary.items():
            if isinstance(value, dict):
                print(f"{' ' * indent}{key}:")
                self.print_nested_dict(value, indent + 4)
            else:
                print(f"{' ' * indent}{key}: {value}")
       
    # Used in find_missing_bonds method
    def matches_with_x(self, ff_entry, lig_entry):
        # ff_entry 和 lig_entry 都是包含四个元素的列表，例如 dihedral ['C', 'N', 'CA', 'X']
        for ff_atom, lig_atom in zip(ff_entry, lig_entry):
            if ff_atom != 'X' and ff_atom != lig_atom:
                return False
        return True
    
    def find_missing_bonds(self, lig_dict, ffbonded_dict):
        # 提取 ligand 的键合信息
        lig_bonds_list = [list(pair) for pair in zip(lig_dict['bonds']['a1'], lig_dict['bonds']['a2'])]  
        # lig_angles
        lig_angles_list = [list(pair) for pair in zip(lig_dict['angles']['a1'], lig_dict['angles']['a2'], lig_dict['angles']['a3'])]  
        # lig_dihedrals
        lig_dihedrals_list = [list(pair) for pair in zip(lig_dict['dihedrals']['a1'], lig_dict['dihedrals']['a2'], lig_dict['dihedrals']['a3'], lig_dict['dihedrals']['a4'])]   
        # lig_impropers
        lig_impropers_list = [list(pair) for pair in zip(lig_dict['impropers']['a1'], lig_dict['impropers']['a2'], lig_dict['impropers']['a3'], lig_dict['impropers']['a4'])]          

        # 提取 forcefield 中已有的键合信息
        ff_bonds_list = [list(pair) for pair in zip(ffbonded_dict['bondtypes']['a1'], ffbonded_dict['bondtypes']['a2'])]  
        # ff_angles
        ff_angles_list = [list(pair) for pair in zip(ffbonded_dict['angletypes']['a1'], ffbonded_dict['angletypes']['a2'], ffbonded_dict['angletypes']['a3'])]  
        # ff_dihedrals
        ff_dihedrals_list = []
        for dihedral_dict in [ffbonded_dict['dihedral_f3'], ffbonded_dict['dihedral_f9']]:
            if dihedral_dict:
                dihedrals_list = [list(pair) for pair in zip(dihedral_dict['a1'], dihedral_dict['a2'], dihedral_dict['a3'], dihedral_dict['a4'])]
                ff_dihedrals_list.extend(dihedrals_list)
        ff_impropers_list = []
        for dihedral_dict in [ffbonded_dict['dihedral_f1'], ffbonded_dict['dihedral_f4']]:
            if dihedral_dict:
                dihedrals_list = [list(pair) for pair in zip(dihedral_dict['a1'], dihedral_dict['a2'], dihedral_dict['a3'], dihedral_dict['a4'])]
                ff_impropers_list.extend(dihedrals_list)

        # 找出 forcefield 中缺少的键合信息
        missing_bonds = [bond for bond in lig_bonds_list if bond not in ff_bonds_list and bond[::-1] not in ff_bonds_list]
        missing_angles = [angle for angle in lig_angles_list if angle not in ff_angles_list and angle[::-1] not in ff_angles_list]
        # missing_dihedrals = [dihedral for dihedral in lig_dihedrals_list if dihedral not in ff_dihedrals_list and dihedral[::-1] not in ff_dihedrals_list]
        missing_dihedrals = []
        for lig_dihedral in lig_dihedrals_list:
            if not any(self.matches_with_x(ff_dihedral, lig_dihedral) or self.matches_with_x(ff_dihedral[::-1], lig_dihedral) for ff_dihedral in ff_dihedrals_list):
                missing_dihedrals.append(lig_dihedral)
        # missing_impropers = [dihedral for dihedral in lig_impropers_list if dihedral not in ff_impropers_list and dihedral[::-1] not in ff_impropers_list]
        missing_impropers = []
        for lig_improper in lig_impropers_list:
            if not any(self.matches_with_x(ff_dihedral, lig_improper) or self.matches_with_x(ff_dihedral[::-1], lig_improper) for ff_dihedral in ff_impropers_list):
                missing_impropers.append(lig_improper)
        # 删除重复
        missing_bonds_dedu = []
        missing_angles_dedu =[]
        missing_dihedrals_dedu = []
        missing_impropers_dedu = []
        [missing_bonds_dedu.append(bonds) for bonds in missing_bonds if bonds not in missing_bonds_dedu]
        [missing_angles_dedu.append(bonds) for bonds in missing_angles if bonds not in missing_angles_dedu]
        [missing_dihedrals_dedu.append(bonds) for bonds in missing_dihedrals if bonds not in missing_dihedrals_dedu]
        [missing_impropers_dedu.append(bonds) for bonds in missing_impropers if bonds not in missing_bonds_dedu]
        
        print("######################## Please add below to ffbonded.itp! ########################") 
        print(f"There are in total {len(missing_bonds_dedu)} missed bonds: {missing_bonds_dedu}")
        print(f"There are in total {len(missing_angles_dedu)} missed angles: {missing_angles_dedu}")
        print(f"There are in total {len(missing_dihedrals_dedu)} missed proper dihedrals: {missing_dihedrals_dedu}")
        print(f"There are in total {len(missing_impropers_dedu)} missed improper dihedrals: {missing_impropers_dedu}")
        print("########### Add below to [ bondtypes ] segment ###########")
        count = 0
        for bond in missing_bonds_dedu:
            for a1, a2, funct, bond_length, bond_strength in zip(lig_dict['bonds']['a1'], lig_dict['bonds']['a2'], lig_dict['bonds']['funct'], lig_dict['bonds']['bond_length'], lig_dict['bonds']['bond_strength']):
                if bond[0] == a1 and bond[1] == a2:
                    print("   %-3s%-11s%1d%11.5f%11.1f" % (a1, a2, int(funct), float(bond_length), float(bond_strength)))
                    break
        print("########### Add below to [ angletypes ] segment ###########")
        for angle in missing_angles_dedu:
            for a1, a2, a3, funct, angle_degree, angle_strength in zip(lig_dict['angles']['a1'], lig_dict['angles']['a2'], lig_dict['angles']['a3'], lig_dict['angles']['funct'], lig_dict['angles']['angle_degree'], lig_dict['angles']['angle_strength']):
                if angle[0] == a1 and angle[1] == a2 and angle[2] == a3:
                    print("%-4s%-4s%-13s%1d%10.3f%11.3f" % (a1, a2, a3, int(funct), float(angle_degree), float(angle_strength)))
                    break
        print("########### Add below to [ dihedraltypes ] function = 1 segment ###########")
        for dihedral in missing_impropers_dedu:
            for a1, a2, a3, a4, funct, angle_degree, angle_strength, multiplicity in zip(lig_dict['impropers']['a1'], lig_dict['impropers']['a2'], lig_dict['impropers']['a3'], lig_dict['impropers']['a4'], lig_dict['impropers']['funct'], lig_dict['impropers']['angle_degree'], lig_dict['impropers']['angle_strength'], lig_dict['impropers']['multiplicity']):
                if dihedral[0] == a1 and dihedral[1] == a2 and dihedral[2] == a3 and dihedral[3] == a4:
                    print("%-4s%-4s%-4s%-9s%1d%12.2f%12.5f%6d" % (a1, a2, a3, a4, int(funct), float(angle_degree), float(angle_strength), int(multiplicity)))
                    break
        print("########### Add below to [ dihedraltypes ] function = 3 segment ###########")
        for dihedral in missing_dihedrals_dedu:
            for a1, a2, a3, a4, funct, c0, c1, c2, c3, c4, c5 in zip(lig_dict['dihedrals']['a1'], lig_dict['dihedrals']['a2'], lig_dict['dihedrals']['a3'], lig_dict['dihedrals']['a4'], lig_dict['dihedrals']['funct'], lig_dict['dihedrals']['C0'], lig_dict['dihedrals']['C1'], lig_dict['dihedrals']['C2'], lig_dict['dihedrals']['C3'], lig_dict['dihedrals']['C4'], lig_dict['dihedrals']['C5']):
                if dihedral[0] == a1 and dihedral[1] == a2 and dihedral[2] == a3 and dihedral[3] == a4:
                    print("%-4s%-4s%-4s%-9s%1d%12.6f%11.6f%11.6f%11.6f%11.6f%11.6f" % (a1, a2, a3, a4, int(funct), float(c0), float(c1), float(c2), float(c3), float(c4), float(c5)))
                    break

    # def parse_nonbonded(self, ff_nonbonded):  
        
    def parse_atomtypes(self, ff_aomtypes, ff_nonbonded):  # check if the user provided atom_types has any new atom type, if yes, then ask user the new atomtype's at.num(atomtypes), mass(atomtypes), sigma(ffnonbonded) and epsilon(ffnonbonded), give a guess value based on the database
        # Dictionary to store combined atom data
        atom_data = {}
        
        # Read data from atomtypes.atp
        with open(ff_aomtypes, 'r') as file:
            for line in file:
                if line.strip() and not line.startswith(';'):
                    parts = line.split()
                    atom_type = parts[0]
                    mass_atp = float(parts[1])
                    # Initialize or update the dictionary for this atom type
                    if atom_type not in atom_data:
                        atom_data[atom_type] = {'mass_atp': mass_atp}
                    else:
                        atom_data[atom_type]['mass_atp'] = mass_atp
        
        # Read data from ffnonbonded.itp
        with open(ff_nonbonded, 'r') as file:
            for line in file:
                if line.strip() and not line.startswith(';') and not line.startswith('['):
                    parts = line.split()
                    atom_type = parts[0]
                    at_num = int(parts[1])
                    mass_ffn = float(parts[2])
                    charge = float(parts[3])
                    ptype = parts[4]
                    sigma = float(parts[5])
                    epsilon = float(parts[6])
                    # Initialize or update the dictionary for this atom type
                    if atom_type not in atom_data:
                        atom_data[atom_type] = {
                            'at_num': at_num, 'mass_ffn': mass_ffn, 'charge': charge,
                            'ptype': ptype, 'sigma': sigma, 'epsilon': epsilon
                        }
                    else:
                        atom_data[atom_type].update({
                            'at_num': at_num, 'mass_ffn': mass_ffn, 'charge': charge,
                            'ptype': ptype, 'sigma': sigma, 'epsilon': epsilon
                        })
        
        return atom_data

    def find_missing_atomtypes(self, ff_atomtype, atom_types):
        # Check for missing atom types from the user-provided list
        missing_atom_types = []
        missing_atom_types = [atype for atype in atom_types if atype not in ff_atomtype]
        if missing_atom_types:
            missing_atom_type = list(set(missing_atom_types))
            print(f"The missed atom types are: {missing_atom_type}")
        return missing_atom_type
        
    
    
    
    # def calculate_bond_values(self, lig_gro):
        
##########################################################################################################################################################################################

if multi_files[0] != 0:
    x = plotly_go(multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, violin, x_low, x_high, y_low, y_high, error_bar, replica_number, transparency)
elif mr_file != 0:
    x = mr(mr_file, mr_num_neighbours, mr_distance_value, mr_atom_list, mr_metal_list, mr_residue_list, mr_bond_strength, mr_angle_strength)
elif gm_receptor_gro != 0:
    x = gromerger(gm_receptor_gro, gm_ligand_gro, gm_ligand_itp, gm_receptor_top )
elif cm_topol != 0:
    x = contact_map_detect(cm_topol, cm_traj, cm_lig, cm_output, cm_distance)
elif pl_pdb != 0:
    x = pep2lig(pl_pdb, pl_name)
elif ds_dssp != 0:
    x = gmx_dssp(ds_dssp, ds_traj, ds_output, ds_original, ds_color)
elif rn_file !=0:
    renumber_MODEL(rn_file)
elif ar_lig_itp !=0:
    x = ff_res_adder(ar_lig_itp, ar_lig_gro, ar_ff_rtp, ar_ff_hdb, ar_ff_bonded, ar_ff_nonbonded, ar_ff_aomtypes, ar_ff_restype, ar_atom_types, ar_isamino_acid, ar_output_name)


