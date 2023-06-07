#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 16:22:53 2023

@author: Shangze Xu
"""

import re
import argparse

parser = argparse.ArgumentParser(description='Version 1.0 \n'
                                 'Optimised SASA  \n'
                                 'Usage: ./sasa_cool -z NIST.pdb -x sasa1.xvg sasa2.xvg sasa3.xvg -o nist_sasa_good.txt \n'
                                 'Usage: ./sasa_cool -z NIST.pdb -x sasa1.xvg sasa2.xvg sasa3.xvg -o nist_sasa_good.txt -e 0 \n'
                                 'Usage: ./sasa_cool -z NIST.pdb -x sasa1.xvg sasa2.xvg sasa3.xvg -o nist_sasa_good.txt -co 10 \n',
                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-z', '--pdb', nargs='+', required=True, default=0, help='input gro file')
parser.add_argument('-o', '--output', default="sasa_better.txt", help='input gro file')
parser.add_argument('-x', '--sasa', nargs='+', required=True, default=0, help='e.g: sasa1.xvg sasa2.xvg sasa3.xvg')
parser.add_argument('-co', '--cutoff', default=0, help='e.g: cutoff of the number of residues, default is 0')
parser.add_argument('-r', '--arginine', default=1, help='charges for arginine, default value = +1')
parser.add_argument('-k', '--lysine', default=1, help='charges for lysine, default value = +1')
parser.add_argument('-hi', '--histidine', default=0, help='charges for histidine, default value = 0')
parser.add_argument('-d', '--asparticacid', default=-1, help='charges for asparticacid, default value = -1')
parser.add_argument('-e', '--glutamicacid', default=-1, help='charges for glutamicacid, default value = -1')
parser.add_argument('-s', '--serine', default=0, help='charges for serine, default value = 0')
parser.add_argument('-t', '--threonine', default=0, help='charges for threonine, default value = 0')
parser.add_argument('-n', '--asparagine', default=0, help='charges for asparagine, default value = 0')
parser.add_argument('-q', '--glutamine', default=0, help='charges for glutamine, default value = 0')
parser.add_argument('-c', '--cysteine', default=0, help='charges for cysteine, default value = 0')
parser.add_argument('-u', '--selenocysteine', default=0, help='charges for selenocysteine, default value = 0')
parser.add_argument('-g', '--glycine', default=0, help='charges for glycine, default value = 0')
parser.add_argument('-p', '--proline', default=0, help='charges for proline, default value = 0')
parser.add_argument('-a', '--alanine', default=0, help='charges for alanine, default value = 0')
parser.add_argument('-v', '--valine', default=0, help='charges for valine, default value = 0')
parser.add_argument('-i', '--isoleucine', default=0, help='charges for isoleucine, default value = 0')
parser.add_argument('-l', '--leucine', default=0, help='charges for leucine, default value = 0')
parser.add_argument('-m', '--methionine', default=0, help='charges for methionine, default value = 0')
parser.add_argument('-f', '--phenylalanine', default=0, help='charges for phenylalanine, default value = 0')
parser.add_argument('-y', '--tyrosine', default=0, help='charges for tyrosine, default value = 0')
parser.add_argument('-w', '--tryptophan', default=0, help='charges for tryptophan, default value = 0')


args = parser.parse_args()
file = args.pdb[0]
sasa_list = [str(x) for x in args.sasa]
output_name = args.output
cutoff = int(args.cutoff)
# charge_dic = {'arginine' : args.arginine,
#               'lysine'   : args.lysine,
#               'histidine': args.lysine,
#               'asparticacid' : args.asparticacid,
#               'glutamicacid' : args.glutamicacid,
#               'serine' : args.serine,
#               'threonine' : args.threonine,
#               'asparagine' : args.asparagine,
#               'glutamine' : args.glutamine,
#               'cysteine' : args.cysteine,
#               'selenocysteine' : args.selenocysteine,
#               'glycine' : args.glycine,
#               'proline' : args.proline,
#               'alanine' : args.alanine,
#               'valine' : args.valine,
#               'isoleucine' : args.isoleucine,
#               'leucine' : args.leucine,
#               'methionine' : args.methionine,
#               'phenylalanine' : args.phenylalanine,
#               'tyrosine' : args.tyrosine,
#               'tryptophan' : args.tryptophan,
#     }

charge_dic = {'ARG' : int(args.arginine),
              'LYS' : int(args.lysine),
              'HIS': int(args.histidine),
              'ASP' : int(args.asparticacid),
              'GLU' : int(args.glutamicacid),
              'SER' : int(args.serine),
              'THR' : int(args.threonine),
              'ASN' : int(args.asparagine),
              'GLN' : int(args.glutamine),
              'CYS' : int(args.cysteine),
              'SEC' : int(args.selenocysteine),
              'GLY' : int(args.glycine),
              'PRO' : int(args.proline),
              'ALA' : int(args.alanine),
              'VAL' : int(args.valine),
              'ILE' : int(args.isoleucine),
              'LEU' : int(args.leucine),
              'MET' : int(args.methionine),
              'PHE' : int(args.phenylalanine),
              'TYR' : int(args.tyrosine),
              'TRP' : int(args.tryptophan),
    }


class sasa_better:
    atomic_index = []
    atomic_name = []
    residue_name = []
    chain_name = []
    residue_index = []
    
    residue_list = []
    sasa_values = []
    mean_sasa_values = []
    nor_mean_sasa_values = []
    charge_list = []
    produc = []
        
        
    def __init__(self, pdb, sasa_list, charge_dic, output_name, cutoff):
        self.PDBreader(pdb, cutoff)
        self.SASAreader(sasa_list, cutoff)
        self.residue_charges(charge_dic)
        self.save_file(output_name)
        
    def PDBreader(self, pdb, cutoff):                                                   # first method to be called in __main__, used for creating object and charactors.
        self.atomic_index.clear()
        self.atomic_name.clear()
        self.residue_name.clear()
        self.chain_name.clear()
        self.residue_index.clear()
        # self.X_peratom.clear()
        # self.Y_peratom.clear()
        # self.Z_peratom.clear()
        # self.bfactor_per_factor.clear()
        # self.charge_per_factor.clear()
        # self.Atomtype_per_atom.clear()
        f = open(pdb, "r")                                                     # "filename" = $PATH/crankpep_docking_results/PLDAYL_corrected_top_1.pdb
        for line in f:                                                              # iterate each line in file "f"                     
            if(line.split()[0] == "ATOM" or line.split()[0] == "HETATM"):       # Judgment Sentence锛孶sed to split each row and then determine whether the first column of the row == ATOM or HETATM
                self.atomic_index.append(float(line.split()[1]))                # The second column is the atomic number
                self.atomic_name.append(line.split()[2])                        # The 3rd column is the atom name C CA CD1 CD2 and so on
                self.residue_name.append(line.split()[3])                       # Column 4 is the residue name TYR ALA etc.       
                pattern1 = r'([A-Z]+)(\s{1})([A-Za-z]{1})(\s+)(\d+)'
                pattern2 = r'([A-Z]+)(\s{1})([A-Za-z]{1})(\d+)'
                matches1 = re.findall(pattern1, line)  
                matches2 = re.findall(pattern2, line)
                # print(matches1, matches2)
                if matches1:
                    position, position, chain, position, resid = matches1[0]
                    self.chain_name.append(chain)
                    self.residue_index.append(resid)
                    pass
                elif matches2:
                    position, position, chain, resid = matches2[0]
                    self.chain_name.append(chain)
                    self.residue_index.append(resid)
                    pass
        
        
        # generating the residue list              
        for i in range(len(self.atomic_index)-1):
            if self.residue_index[i+1] != self.residue_index[i]:
                self.residue_list.append(self.residue_name[i])
        self.residue_list.append(self.residue_name[-1])
                
        if cutoff == 0:
            pass
        elif cutoff !=0:
            self.residue_list = self.residue_list[:cutoff]
        
                
    def SASAreader(self, sasa_list, cutoff):
        for i in sasa_list:
            locals()[i] = []           
            with open(i, 'r') as file:
                for line in file:
                    if not line.startswith('#') and not line.startswith('@'):
                        values = line.split()
                        if len(values) >= 2:
                            locals()[i].append(float(values[1]))
            
            # add the values to sasa_values
            self.sasa_values.append(locals()[i]) 
        # calculate the mean value for sasa replicas    
        self.mean_sasa_values = [sum(sublist) / len(sublist) for sublist in zip(*self.sasa_values)]
        # get the max value
        max_sasa = max(self.mean_sasa_values)
        # normalise
        for i in self.mean_sasa_values:
            self.nor_mean_sasa_values.append(i/max_sasa)
        
        if cutoff == 0:
            pass
        elif cutoff != 0:
            self.mean_sasa_values = self.mean_sasa_values[:cutoff]
            self.nor_mean_sasa_values = self.nor_mean_sasa_values[:cutoff]
        
    def residue_charges(self, charge_dic):
        for i in range(len(self.residue_list)):
            self.charge_list.append(charge_dic[self.residue_list[i]])
        
        # product
        self.produc = [x * y for x, y in zip(self.nor_mean_sasa_values, self.charge_list)]
            
    def save_file(self, output_name):
        with open(output_name, 'w') as w:
            w.write("%-15s %15s %15s %15s %15s %15s" %  ("Residue Index", "Residue", "Mean value", "Normalised Mean", "Net Charge", "Product"))
            w.write("\n")
            a = 1
            for i in range(len(self.residue_list)):
                w.write("%-15d %15s %15.2f %15.2f %15d %15.2f" %  (a, self.residue_list[i], self.mean_sasa_values[i], self.nor_mean_sasa_values[i], self.charge_list[i], self.produc[i]))
                w.write("\n")
                a += 1
            w.write("sum of the product = " + str(round(sum(self.produc), 2)))
                       
x = sasa_better(file, sasa_list, charge_dic, output_name, cutoff)