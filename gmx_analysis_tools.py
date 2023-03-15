#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 10:15:35 2023

@author: Shang

RMSD start at line 19
RMSF start at line 18
"""

import getopt
import sys
import csv

args=sys.argv[1:]  
file1 = ''
file2 = ''
file3 = ''
renumber = 'fault'

try:
    opts, args = getopt.getopt(args,"f:s:t:r:h",["file1=",
                                             "file2=",
                                             "file3=",
                                             "renumber="
                                             "help"])
except getopt.GetoptError:
    print('Usage: ./gmx_analysis_tools -f <file1> -s <file2> -t <file3> -r <true/fault>  \nCurrently it works for rmsd, rmsf, sasa')
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print('Usage: ./gmx_analysis_tools -f <file1> -s <file2> -t <file3> -r <true/fault> \nCurrently it works for rmsd, rmsf, sasa')
        sys.exit()
    elif opt in ("-f", "--file1"):
        file1 = str(arg)
    elif opt in ("-s", "--file2"):
        file2 = str(arg)
    elif opt in ("-t", "--file3"):
        file3 = str(arg)
    elif opt in ("-r", "--renumber"):
        renumber = str(arg)

class merger():
    flag = ''
    sasa_flag = ''
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
    
    def __init__(self,file1, file2, file3, renumber):
        self.flag_recognizer(file1, file2, file3)
        if self.flag == 'rmsd':
            self.rmsd_averager(file1,file2,file3)
        elif self.flag == 'rmsf':
            self.rmsf_averager(file1, file2,file3)
        elif self.flag == 'sasa' and self.sasa_flag != '-or':
            self.sasa_averager(file1, file2, file3)
        elif self.flag == 'sasa' and self.sasa_flag == '-or':
            self.sasa_res_averager(file1, file2, file3,renumber)
        
    def flag_recognizer(self,file1, file2, file3):                                                   # first method to be called in __main__, used for creating object and charactors.
        #self.atomic_index.clear()                                                   # cleaveage the information of previous object before put new record into these charactors
        with open(file1, 'r') as f:
            lines = f.readlines()                                                 # "filename" = $PATH/crankpep_docking_results/PLDAYL_corrected_top_1.pdb
            if len(lines)  >= 3:
                self.flag = lines[2].split()[5]
                if self.flag == 'rms,' :
                    self.flag = 'rmsd'
                elif self.flag == 'rmsf,':
                    self.flag = 'rmsf'
                elif self.flag == 'sasa,':
                    self.flag = 'sasa'
            if len(lines) >= 9 and '-or' in lines[8]:
                self.sasa_flag = '-or'
                
                
                    
    def rmsd_averager(self, file1, file2, file3):
        a = 18
        with open(file1, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
        with open(file2, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time2.append(float(lines[i].split()[0]))
                    self.values2.append(float(lines[i].split()[1]))
        with open(file3, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time3.append(float(lines[i].split()[0]))
                    self.values3.append(float(lines[i].split()[1]))
                    
        for i in [self.values1, self.values2, self.values3]:
            self.max_value.append(max(i))
            self.average_value.append(sum(i)/len(i))
        with open('rmsd_ave_max.txt', 'w') as w:
            w.write('    No.         average_value   max_value')
            for i in range(3):
               w.write(f"\n replica{i}    {self.average_value[i]:-8.3f}{self.max_value[i]:16.3f}")
        
        data = [["Replica_No.", "Average_value", "Max_value"]]
        for i in range(len(self.max_value)):
            data.append([i, self.average_value[i], self.max_value[i]])
        with open('rmsd_ave_max.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(data)
        

                    
    def rmsf_averager(self, file1, file2, file3):
        a = 17
        with open(file1, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
        with open(file2, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time2.append(float(lines[i].split()[0]))
                    self.values2.append(float(lines[i].split()[1]))
        with open(file3, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time3.append(float(lines[i].split()[0]))
                    self.values3.append(float(lines[i].split()[1]))

        with open('rmsf_merged.txt', 'w') as w:
            w.write("%4s %12s %12s %12s" %  ("resid", "replica1", "replica2", "replica3"))
            for i in range(len(self.values1)):
                w.write('\n')
                w.write("%4d %12.3f %12.3f %12.3f" %  (self.time1[i], self.values1[i], self.values2[i], self.values3[i]))
                
        data = [["Resid", "Replica1", "Replica2", "Replica3"]]
        for i in range(len(self.time1)):
            data.append([self.time1[i], self.values1[i], self.values2[i], self.values3[i]])
        with open('rmsf_merged.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(data)
           
            
    def sasa_averager(self, file1, file2, file3):
        a = 24
        with open(file1, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
        with open(file2, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time2.append(float(lines[i].split()[0]))
                    self.values2.append(float(lines[i].split()[1]))
        with open(file3, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time3.append(float(lines[i].split()[0]))
                    self.values3.append(float(lines[i].split()[1]))

        with open('sasa_merged.txt', 'w') as w:
            w.write("%4s %12s %12s %12s" %  ("time", "replica1", "replica2", "replica3"))
            for i in range(len(self.values1)):
                w.write('\n')
                w.write("%4d %12.3f %12.3f %12.3f" %  (self.time1[i], self.values1[i], self.values2[i], self.values3[i]))
                
        data = [["Time", "Replica1", "Replica2", "Replica3"]]
        for i in range(len(self.time1)):
            data.append([self.time1[i], self.values1[i], self.values2[i], self.values3[i]])
        with open('sasa_merged.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(data)
            
    def sasa_res_averager(self, file1, file2, file3,renumber):
        a = 25
        with open(file1, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
                    self.sd1.append(float(lines[i].split()[2]))
        with open(file2, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time2.append(float(lines[i].split()[0]))
                    self.values2.append(float(lines[i].split()[1]))
                    self.sd2.append(float(lines[i].split()[2]))
        with open(file3, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time3.append(float(lines[i].split()[0]))
                    self.values3.append(float(lines[i].split()[1]))
                    self.sd3.append(float(lines[i].split()[2]))
        if renumber == 'true':
            discontinuous_positions = []
            for i in range(1, len(self.time1)):
                if self.time1[i] != self.time1[i-1]+1:
                    discontinuous_positions.append(i)
            for i in range(len(self.time1)):
                self.time1[i] = i+1

         
        data = [['# discontinuous_positions are here: ' + str(discontinuous_positions)], ["ResID", "Replica1", "Replica2", "Replica3"]]
        for i in range(len(self.time1)):
            data.append([self.time1[i], self.values1[i], self.values2[i], self.values3[i]])
        with open('sasa_res_average.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(data)
            
        data = [['# discontinuous_positions are here: ' + str(discontinuous_positions)], ["ResID", "Replica1", "Replica2", "Replica3"]]
        for i in range(len(self.time1)):
            data.append([self.time1[i], self.sd1[i], self.sd2[i], self.sd3[i]])
        with open('sasa_res_sd.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(data)

                

                    
x = merger(file1,file2,file3,renumber)
