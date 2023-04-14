#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:15:48 2023

@author: dozeduck
"""
import getopt
import sys
import csv
import plotly.graph_objs as go

args=sys.argv[1:]  
file1 = ''
file2 = ''
file3 = ''
renumber = 'fault'
ave = 'fault'
output_name = 'plotly'
font_size = 46
font_family = 'Arial'
font_color = 'black'

try:
    opts, args = getopt.getopt(args,"f:s:t:r:a:o:h",["file1=",
                                             "file2=",
                                             "file3=",
                                             "renumber=",
                                             "average=",
                                             "output_name="
                                             "help"])
except getopt.GetoptError:
    print('Usage: ./gmx_analysis_tools -f <file1> -s <file2> -t <file3> -o <output_name> -r <true/fault> -a <true/fault> \
          \nCurrently it works for rmsd, rmsf, sasa \
          \n-r true: default is fault, renumber the residues, mainly for the uncontinues residue numbers while work on rmsf per residue\
          \n-a true: default is fault, output the average score for each column at the bottom line, mainly for rmsd and radius of gyration calculation')
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print('Usage: ./gmx_analysis_tools -f <file1> -s <file2> -t <file3> -r <true/fault> -a <true/fault> \
          \nCurrently it works for rmsd, rmsf, sasa \
          \n-r true: default is fault, renumber the residues, mainly for the uncontinues residue numbers while work on rmsf per residue\
          \n-a true: default is fault, output the average score for each column at the bottom line, mainly for rmsd and radius of gyration calculation')
        sys.exit()
    elif opt in ("-f", "--file1"):
        file1 = str(arg)
    elif opt in ("-s", "--file2"):
        file2 = str(arg)
    elif opt in ("-t", "--file3"):
        file3 = str(arg)
    elif opt in ("-r", "--renumber"):
        renumber = str(arg)
    elif opt in ("-a", "--average"):
        ave = str(arg)
    elif opt in ("-o", "--output_name"):
        output_name = str(arg)
        
        
class plotly_go():
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
    
    def __init__(self,file1, file2, file3, output_name, renumber, ave):
        self.flag_recognizer(file1, file2, file3)
        if self.flag == 'rmsd':
            self.rmsd_averager(file1,file2,file3, ave, output_name)
        elif self.flag == 'rmsf':
            self.rmsf_averager(file1, file2,file3, output_name)
        elif self.flag == 'sasa' and self.sasa_flag != '-or':
            self.sasa_averager(file1, file2, file3, output_name)
        elif self.flag == 'sasa' and self.sasa_flag == '-or':
            self.sasa_res_averager(file1, file2, file3,renumber,ave, output_name)
        elif self.flag == 'gyrate':
            self.gyrate_averager(file1, file2,file3, ave, output_name)
        elif self.flag == 'dipoles':
            self.dipoles_averager(file1, file2,file3, ave, output_name)
        
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
                elif self.flag == 'gyrate,':
                    self.flag = 'gyrate'
                elif self.flag == 'dipoles,':
                    self.flag = 'dipoles'
            if len(lines) >= 9 and '-or' in lines[8]:
                self.sasa_flag = '-or'

                
                
                    
    def rmsd_averager(self, file1, file2, file3, ave, output_name):
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
        
        # calculate the average and max value for each replica
        for i in [self.values1, self.values2, self.values3]:
            self.max_value.append(max(i))
            self.average_value.append(sum(i)/len(i))
        
        # define time unit is ps or ns
        with open(file1, 'r') as f: # define time unit
            lines = f.readlines() 
            if len(lines) >= 15 and '(ns)"' in lines[14]:
                timeflag = 'ns'
            elif len(lines) >= 15 and '(ps)"' in lines[14]:
                timeflag = 'ps'
        if timeflag == 'ps':
            divisor = 1000
            self.time1 = [x / divisor for x in self.time1]
            

        # Line chart
        plot_title = 'RMSD'
        x_name = 'Time (ns)'
        y_name = 'RMSD (nm)'
        self.plotly_dawang(plot_title, x_name, y_name, self.time1, self.values1, self.values2, self.values3)     
        #trace1 = go.Scatter(x=self.time1, y=self.values1, name='Replica 1')
        #trace2 = go.Scatter(x=self.time1, y=self.values2, name='Replica 2')
        #trace3 = go.Scatter(x=self.time1, y=self.values3, name='Replica 3')
        #data = [trace1, trace2, trace3]
        #layout = go.Layout(title='RMSD',
        #                   xaxis=dict(title='Time (ns)', titlefont=dict(size=46, color='black', family='Arial'), zeroline=False, autorange=True,
        #                              showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
        #                   yaxis=dict(title='RMSD (nm)', titlefont=dict(size=46, color='black', family='Arial'), zeroline=False, autorange=True,
        #                              showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
        #                   legend=dict(x=1, y=1, orientation='v', font=dict(size=30)),
        #                   plot_bgcolor='rgba(255, 255, 255, 0.1)',
        #                   paper_bgcolor='rgba(255, 255, 255, 0.2)')
        #fig = go.Figure(data=data, layout=layout)
        #fig.show()
        
        

                    
    def rmsf_averager(self, file1, file2, file3, output_name):
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

        # define axis is atom or residues
        with open(file1, 'r') as f: # define time unit
            lines = f.readlines() 
            if len(lines) >= 15 and '"Residue"' in lines[14]:
                x_name = 'Residue Number'
            elif len(lines) >= 15 and '"Atom"' in lines[14]:
                x_name = 'Atom Number'
            
        plot_title = 'RMSF'
        y_name = 'RMSF (nm)'
        self.plotly_dawang(plot_title, x_name, y_name, self.time1, self.values1, self.values2, self.values3)            

        # Line chart
        #trace1 = go.Scatter(x=self.time1, y=self.values1, name='Replica 1')
        #trace2 = go.Scatter(x=self.time1, y=self.values2, name='Replica 2')
        #trace3 = go.Scatter(x=self.time1, y=self.values3, name='Replica 3')
        #data = [trace1, trace2, trace3]
        #layout = go.Layout(title='RMSF',
        #                   xaxis=dict(title=xaxis_name, titlefont=dict(size=46, color='black', family='Arial'), zeroline=False, autorange=True,
        #                              showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
        #                   yaxis=dict(title='RMSF (nm)', titlefont=dict(size=46, color='black', family='Arial'), zeroline=False, autorange=True,
        #                              showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
        #                   legend=dict(x=1, y=1, orientation='v', font=dict(size=30)),
        #                   plot_bgcolor='rgba(255, 255, 255, 0.1)',
        #                   paper_bgcolor='rgba(255, 255, 255, 0.2)')
        #fig = go.Figure(data=data, layout=layout)
        #fig.show()



           
            
    def sasa_averager(self, file1, file2, file3, output_name):
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


            
    def sasa_res_averager(self, file1, file2, file3,renumber,ave, output_name):
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


        
            
    def gyrate_averager(self, file1, file2, file3, ave, output_name):
        a = 27
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
        

            
    def dipoles_averager(self, file1, file2, file3, ave, output_name):
        a = 27
        with open(file1, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[4]))
        with open(file2, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time2.append(float(lines[i].split()[0]))
                    self.values2.append(float(lines[i].split()[4]))
        with open(file3, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time3.append(float(lines[i].split()[0]))
                    self.values3.append(float(lines[i].split()[4]))
                    
        for i in [self.values1, self.values2, self.values3]:
            self.max_value.append(max(i))
            self.average_value.append(sum(i)/len(i))


    def plotly_dawang(self, plot_title, x_name, y_name, time1, values1, values2, values3):
            trace1 = go.Scatter(x=time1, y=values1, name='Replica 1')
            trace2 = go.Scatter(x=time1, y=values2, name='Replica 2')
            trace3 = go.Scatter(x=time1, y=values3, name='Replica 3')
            data = [trace1, trace2, trace3]
            layout = go.Layout(title=plot_title,
                               xaxis=dict(title=x_name, titlefont=dict(size=46, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               yaxis=dict(title=y_name, titlefont=dict(size=46, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               legend=dict(x=1, y=1, orientation='v', font=dict(size=30)),
                               plot_bgcolor='rgba(255, 255, 255, 0.1)',
                               paper_bgcolor='rgba(255, 255, 255, 0.2)')
            fig = go.Figure(data=data, layout=layout)
            fig.show()        
            
                

                    
x = plotly_go(file1,file2,file3, output_name, renumber, ave)