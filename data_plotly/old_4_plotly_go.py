#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:15:48 2023

@author: dozeduck
"""
import getopt
import sys
# import plotly
import plotly.graph_objs as go
import plotly.io as pio

args=sys.argv[1:]  
file1 = ''
file2 = ''
file3 = ''
renumber = 'fault'
ave = 'fault'
output_name = 'plotly.png'
font_size = 40
font_family = 'Arial'
font_color = 'black'
xaxis_name=0
yaxis_name=0

try:
    opts, args = getopt.getopt(args,"f:s:t:r:a:o:x:y:h",["file1=",
                                             "file2=",
                                             "file3=",
                                             "renumber=",
                                             "average=",
                                             "output_name=",
                                             "xaxis_name=",
                                             "yaxis_name=",
                                             "help"])
except getopt.GetoptError:
    print('Usage: ./plotly_go -f <file1> -s <file2> -t <file3> -o <output_name> -r <true/fault> -a <2/3> -x <xaxis_name> -y <yaxis_name> \
          \nCurrently it works for rmsd, rmsf, sasa_time, sasa_residue, gyration, dipole movement, distance  \
          \nIt can read one, two or three same type files \
          \n-o output_name.png, suitable file format for output: png, jpg, jpeg, webp, svg, pdf, eps, json \
          \n-r true: default is fault, renumber the residues, mainly for the duplicated residue numbers while work on rmsf and sasa per residue\
          \n-a 2/3: default is fault, output the average score for each replica\
          \nE.g: \
          \n ./plotly_go -f rmsd1.xvg -o rmsd1.png \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -o rmsd12.png \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -t rmsd3.xvg -o rmsd123.png \
          \n ./plotly_go -f rmsf1.xvg -s rmsf2.xvg -t rmsf3.xvg -o rmsf123.png -r true \
          \n ./plotly_go -f rmsf1.xvg -s rmsf2.xvg -t rmsf3.xvg -o rmsf123.png -a 3 \
          \n ./plotly_go -f rmsf1.xvg -s rmsf2.xvg -o rmsf12.png -a 2 \
          \n ./plotly_go -f sasa1.xvg -s sasa2.xvg -t sasa3.xvg -o sasa123.png -r true -x "Time (ns)" -y "SASA (nm<sup>2</sup>)"')
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print('Usage: ./plotly_go -f <file1> -s <file2> -t <file3> -o <output_name> -r <true/fault> -a <2/3> -x <xaxis_name> -y <yaxis_name> \
          \nCurrently it works for rmsd, rmsf, sasa_time, sasa_residue, gyration, dipole movement, distance  \
          \nIt can read one, two or three same type files \
          \n-o output_name.png, suitable file format for output: png, jpg, jpeg, webp, svg, pdf, eps, json \
          \n-r true: default is fault, renumber the residues, mainly for the duplicated residue numbers while work on rmsf and sasa per residue\
          \n-a 2/3: default is fault, output the average score for each replica \
          \nE.g: \
          \n ./plotly_go -f rmsd1.xvg -o rmsd1.png \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -o rmsd12.png \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -o rmsd12.png -a 2 \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -t rmsd3.xvg -o rmsd123.png \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -t rmsd3.xvg -o rmsd123.png -a 3\
          \n ./plotly_go -f sasa1.xvg -s sasa2.xvg -t sasa3.xvg -o sasa123.png -r true -x "Time (ns)" -y "SASA (nm<sup>2</sup>)"')
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
        ave = int(arg)
    elif opt in ("-o", "--output_name"):
        output_name = str(arg)
    elif opt in ("-x", "--xaxis_name"):
        xaxis_name = str(arg)
    elif opt in ("-y", "--yaxis_name"):
        yaxis_name = str(arg)
        
        
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
    
    def __init__(self,file1, file2, file3, output_name, renumber, ave, xaxis_name, yaxis_name):
        self.flag_recognizer(file1, file2, file3)
        if self.flag == 'rmsd':
            self.rmsd_averager(file1,file2,file3, ave, output_name, xaxis_name, yaxis_name)
        elif self.flag == 'rmsf':
            self.rmsf_averager(file1, file2,file3, renumber, ave, output_name, xaxis_name, yaxis_name)
        elif self.flag == 'sasa' and self.sasa_flag != '-or':
            self.sasa_averager(file1, file2, file3, ave, output_name, xaxis_name, yaxis_name)
        elif self.flag == 'sasa' and self.sasa_flag == '-or':
            self.sasa_res_averager(file1, file2, file3,renumber,ave, output_name, xaxis_name, yaxis_name)
        elif self.flag == 'gyrate':
            self.gyrate_averager(file1, file2,file3, ave, output_name, xaxis_name, yaxis_name)
        elif self.flag == 'dipoles':
            self.dipoles_averager(file1, file2,file3, ave, output_name, xaxis_name, yaxis_name)
        elif self.flag == 'distance':
            self.distance_averager(file1, file2,file3, ave, output_name, xaxis_name, yaxis_name)
        
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
                elif self.flag == 'distance,':
                    self.flag = 'distance'

            if len(lines) >= 9 and '-or' in lines[8]:
                self.sasa_flag = '-or'

                
                
                    
    def rmsd_averager(self, file1, file2, file3, ave, output_name, xaxis_name, yaxis_name):
        a = 18
        with open(file1, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
        try:
            with open(file2, 'r') as f:
                lines = f.readlines() 
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time2.append(float(lines[i].split()[0]))
                        self.values2.append(float(lines[i].split()[1]))
        except:
            self.values2 = ''
            pass
        try:
            with open(file3, 'r') as f:
                lines = f.readlines() 
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time3.append(float(lines[i].split()[0]))
                        self.values3.append(float(lines[i].split()[1]))
        except:
            self.values3 = ''
            pass
        
        # calculate the average and max value for each replica
        try:
            for i in [self.values1, self.values2, self.values3]:
                self.max_value.append(max(i))
                self.average_value.append(sum(i)/len(i))
        except:
            pass
        
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
        if xaxis_name == 0:
            x_name = 'Time (ns)'
        else:
            x_name = xaxis_name
        if yaxis_name == 0:
            y_name = 'RMSD (nm)'
        else:
            y_name = yaxis_name
        self.plotly_noSD(output_name, plot_title, x_name, y_name, self.time1, self.values1, self.values2, self.values3, ave)     

        

                    
    def rmsf_averager(self, file1, file2, file3, renumber, ave, output_name, xaxis_name, yaxis_name):
        a = 17
        with open(file1, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
        try:
            with open(file2, 'r') as f:
                lines = f.readlines() 
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time2.append(float(lines[i].split()[0]))
                        self.values2.append(float(lines[i].split()[1]))
        except:
            self.values2 = ''
            pass
        try:
            with open(file3, 'r') as f:
                lines = f.readlines() 
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time3.append(float(lines[i].split()[0]))
                        self.values3.append(float(lines[i].split()[1]))
        except:
            self.values3 = ''
            pass

        if renumber == 'true':
            discontinuous_positions = []
            for i in range(1, len(self.time1)):
                if self.time1[i] != self.time1[i-1]+1:
                    discontinuous_positions.append(i)
            for i in range(len(self.time1)):
                self.time1[i] = i+1

        # define axis is atom or residues
        if xaxis_name == 0:
            with open(file1, 'r') as f: # define time unit
                lines = f.readlines() 
                if len(lines) >= 15 and '"Residue"' in lines[14]:
                    x_name = 'Residue Number'
                elif len(lines) >= 15 and '"Atom"' in lines[14]:
                    x_name = 'Atom Number'
        else:
            x_name = xaxis_name
            
        plot_title = 'RMS fluctuation'
        if yaxis_name == 0:
            y_name = 'RMSF (nm)'
        else:
            y_name = yaxis_name
        self.plotly_noSD(output_name, plot_title, x_name, y_name, self.time1, self.values1, self.values2, self.values3, ave)            




           
            
    def sasa_averager(self, file1, file2, file3, ave, output_name, xaxis_name, yaxis_name):
        a = 24
        with open(file1, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
        try:
            with open(file2, 'r') as f:
                lines = f.readlines() 
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time2.append(float(lines[i].split()[0]))
                        self.values2.append(float(lines[i].split()[1]))
        except:
            self.values2 = ''
            pass
        try:
            with open(file3, 'r') as f:
                lines = f.readlines() 
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time3.append(float(lines[i].split()[0]))
                        self.values3.append(float(lines[i].split()[1]))
        except:
            self.values3 = ''
            pass
        
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
        plot_title = 'Solvent Accessible Surface'
        if xaxis_name == 0:
            x_name = 'Time (ns)'
        else:
            x_name = xaxis_name
        if yaxis_name == 0:
            y_name = 'SASA Area (nm<sup>2</sup>)'
        else:
            y_name = yaxis_name
        self.plotly_noSD(output_name, plot_title, x_name, y_name, self.time1, self.values1, self.values2, self.values3, ave)  


            
    def sasa_res_averager(self, file1, file2, file3,renumber,ave, output_name, xaxis_name, yaxis_name):
        a = 25
        with open(file1, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
                    self.sd1.append(float(lines[i].split()[2]))
        try:
            with open(file2, 'r') as f:
                lines = f.readlines() 
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time2.append(float(lines[i].split()[0]))
                        self.values2.append(float(lines[i].split()[1]))
                        self.sd2.append(float(lines[i].split()[2]))
        except:
            self.values2 = ''
            self.sd2 = ''
            pass
        try:
            with open(file3, 'r') as f:
                lines = f.readlines() 
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time3.append(float(lines[i].split()[0]))
                        self.values3.append(float(lines[i].split()[1]))
                        self.sd3.append(float(lines[i].split()[2]))
        except:
            self.values3 = ''
            self.sd3 = ''
            pass
        
        if renumber == 'true':
            discontinuous_positions = []
            for i in range(1, len(self.time1)):
                if self.time1[i] != self.time1[i-1]+1:
                    discontinuous_positions.append(i)
            for i in range(len(self.time1)):
                self.time1[i] = i+1
            

        # Line chart
        plot_title_noSD = 'Area per residue over the trajectory (noSD)'
        plot_title_SD = 'Area per residue over the trajectory (SD)'
        if xaxis_name == 0:
            x_name = 'Residue Number'
        else:
            x_name = xaxis_name
        if yaxis_name == 0:
            y_name = 'SASA Area (nm<sup>2</sup>)'
        else:
            y_name = yaxis_name
        self.plotly_noSD(output_name, plot_title_noSD, x_name, y_name, self.time1, self.values1, self.values2, self.values3, ave)
        self.plotly_SD(output_name, plot_title_SD, x_name, y_name, self.time1, self.values1, self.values2, self.values3, self.sd1, self.sd2, self.sd3, ave)  


        
            
    def gyrate_averager(self, file1, file2, file3, ave, output_name, xaxis_name, yaxis_name):
        a = 27
        with open(file1, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
        try:
            with open(file2, 'r') as f:
                lines = f.readlines() 
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time2.append(float(lines[i].split()[0]))
                        self.values2.append(float(lines[i].split()[1]))
        except:
            self.values2 = ''
            pass
        try:
            with open(file3, 'r') as f:
                lines = f.readlines() 
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time3.append(float(lines[i].split()[0]))
                        self.values3.append(float(lines[i].split()[1]))
        except:
            self.values3 = ''
            pass
                    
        try:
            for i in [self.values1, self.values2, self.values3]:
                self.max_value.append(max(i))
                self.average_value.append(sum(i)/len(i))
        except:
            pass
        
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
        plot_title = 'Radius of Gyration'
        if xaxis_name == 0:
            x_name = 'Time (ns)'
        else:
            x_name = xaxis_name
        if yaxis_name == 0:
            y_name = 'Rg (nm)'
        else:
            y_name = yaxis_name
        self.plotly_noSD(output_name, plot_title, x_name, y_name, self.time1, self.values1, self.values2, self.values3, ave) 
        

            
    def dipoles_averager(self, file1, file2, file3, ave, output_name, xaxis_name, yaxis_name):
        a = 27
        with open(file1, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[4]))
        try:
            with open(file2, 'r') as f:
                lines = f.readlines() 
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time2.append(float(lines[i].split()[0]))
                        self.values2.append(float(lines[i].split()[4]))
        except:
            self.values2 = ''
            pass
        try:
            with open(file3, 'r') as f:
                lines = f.readlines() 
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time3.append(float(lines[i].split()[0]))
                        self.values3.append(float(lines[i].split()[4]))
        except:
            self.values3 = ''
            pass
                    
        try:
            for i in [self.values1, self.values2, self.values3]:
                self.max_value.append(max(i))
                self.average_value.append(sum(i)/len(i))
        except:
            pass
            
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
        plot_title = 'Total dipole moment of the simulation box vs. time'
        if xaxis_name == 0:
            x_name = 'Time (ns)'
        else:
            x_name = xaxis_name
        if yaxis_name == 0:
            y_name = 'Total Dipole Moment (Debye)'
        else:
            y_name = yaxis_name
        
        
        
        self.plotly_noSD(output_name, plot_title, x_name, y_name, self.time1, self.values1, self.values2, self.values3, ave) 


    def distance_averager(self, file1, file2, file3, ave, output_name, xaxis_name, yaxis_name):
        a = 24
        with open(file1, 'r') as f:
            lines = f.readlines() 
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
        try:
            with open(file2, 'r') as f:
                lines = f.readlines() 
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time2.append(float(lines[i].split()[0]))
                        self.values2.append(float(lines[i].split()[1]))
        except:
            self.values2 = ''
            pass
        try:
            with open(file3, 'r') as f:
                lines = f.readlines() 
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time3.append(float(lines[i].split()[0]))
                        self.values3.append(float(lines[i].split()[1]))
        except:
            self.values3 = ''
            pass
        
        # calculate the average and max value for each replica
        try:
            for i in [self.values1, self.values2, self.values3]:
                self.max_value.append(max(i))
                self.average_value.append(sum(i)/len(i))
        except:
            pass
        
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
        plot_title = 'Distance'
        if xaxis_name == 0:
            x_name = 'Time (ns)'
        else:
            x_name = xaxis_name
        if yaxis_name == 0:
            y_name = 'Distance (nm)'
        else:
            y_name = yaxis_name
        self.plotly_noSD(output_name, plot_title, x_name, y_name, self.time1, self.values1, self.values2, self.values3, ave) 


    def plotly_noSD(self, output_file_name, plot_title, x_name, y_name, time1, values1, values2, values3, ave):
            trace1 = go.Scatter(x=time1, y=values1, name='Replica 1')
            try:
                trace2 = go.Scatter(x=time1, y=values2, name='Replica 2')
            except:
                trace2 = 0
                pass
            try:
                trace3 = go.Scatter(x=time1, y=values3, name='Replica 3')
            except:
                trace3 = 0
                pass
            # data = [trace1, trace2, trace3]
            data = []
            for i in trace1, trace2, trace3:
                if i != 0:
                    data.append(i)

            layout = go.Layout(title=plot_title, title_x=0.5, title_y=1, font=dict(size=24),
                               xaxis=dict(title=x_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               yaxis=dict(title=y_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               legend=dict(x=1, y=1, orientation='v', font=dict(size=30)), showlegend=True,
                               plot_bgcolor='rgba(255, 255, 255, 0.1)',
                               paper_bgcolor='rgba(255, 255, 255, 0.2)')
            fig = go.Figure(data=data, layout=layout)
            # fig.show()
            pio.write_image(fig, output_file_name)
            try:
                if ave == 3:
                    average_value = [(values1[i] + values2[i] + values3[i]) / 3 for i in range(len(values1))] 
                elif ave == 2:
                    average_value = [(values1[i] + values2[i]) / 2 for i in range(len(values1))] 
                trace_ave = trace1 = go.Scatter(x=time1, y=average_value, name='Mean Value')
                data_ave = [trace_ave]
                fig_ave = go.Figure(data=data_ave, layout=layout)
                pio.write_image(fig_ave, 'mean_'+output_file_name)
            except:
                pass
            
            
    def plotly_SD(self, output_file_name, plot_title, x_name, y_name, time1, values1, values2, values3, sd1, sd2, sd3, ave):
            trace1 = go.Scatter(x=time1, y=values1, name='Replica 1')
            error_y1 = dict(type='data', array=sd1, visible=True, thickness=1)
            trace1.update(error_y=error_y1)
            try:
                trace2 = go.Scatter(x=time1, y=values2, name='Replica 2')
                error_y2 = dict(type='data', array=sd2, visible=True, thickness=1)
                trace2.update(error_y=error_y2)
            except:
                trace2 = 0
                pass
            try:
                trace3 = go.Scatter(x=time1, y=values3, name='Replica 3')
                error_y3 = dict(type='data', array=sd3, visible=True, thickness=1)
                trace3.update(error_y=error_y3)
            except:
                trace3 = 0
                pass
            
            # data = [trace1, trace2, trace3]
            data = []
            for i in trace1, trace2, trace3:
                if i != 0:
                    data.append(i)
            layout = go.Layout(title=plot_title, title_x=0.5, title_y=1, font=dict(size=24),
                               xaxis=dict(title=x_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               yaxis=dict(title=y_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               legend=dict(x=1, y=1, orientation='v', font=dict(size=30)),
                               plot_bgcolor='rgba(255, 255, 255, 0.1)',
                               paper_bgcolor='rgba(255, 255, 255, 0.2)')
            fig = go.Figure(data=data, layout=layout)
            # fig.show()
            pio.write_image(fig, 'SD_' + output_file_name)
            try:
                if ave == 3:
                    average_value = [(values1[i] + values2[i] + values3[i]) / 3 for i in range(len(values1))]
                    average_sd = [(sd1[i] + sd2[i] + sd3[i]) / 3 for i in range(len(sd1))]
                elif ave == 2:
                    average_value = [(values1[i] + values2[i]) / 2 for i in range(len(values1))]
                    average_sd = [(sd1[i] + sd2[i]) / 2 for i in range(len(sd1))]
                trace_ave = go.Scatter(x=time1, y=average_value, name='Mean Value')
                error_y_ave = dict(type='data', array=average_sd, visible=True, thickness=1)
                trace_ave.update(error_y=error_y_ave)
                data_ave = [trace_ave]
                fig_ave = go.Figure(data=data_ave, layout=layout)
                pio.write_image(fig_ave, 'SD_mean_'+output_file_name)
            except:
                pass
            
                

                    
x = plotly_go(file1,file2,file3, output_name, renumber, ave, xaxis_name, yaxis_name)