#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:15:48 2023

@author: dozeduck
"""
import getopt
import sys
import re
# import plotly
import plotly.graph_objs as go
import plotly.io as pio
# for PCA
import numpy as np
from rpy2.robjects import r
# for histogram dist plot
import plotly.figure_factory as ff

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
xaxis_name = 0
yaxis_name = 0
rdf_cutoff = 0
multi_files = 0
plot_name = ''
pca = 0
rscript = 0
nbin = 500
size = 400
move_average = 0
mean_value = 'false' # flag = -d
histogram = 'false'

try:
    opts, args_lala = getopt.getopt(args,"f:s:t:r:a:o:x:y:c:m:p:l:j:n:z:b:d:i:h",["file1=",
                                             "file2=",
                                             "file3=",
                                             "renumber=",
                                             "average=",
                                             "output_name=",
                                             "xaxis_name=",
                                             "yaxis_name=",
                                             "rdf_cutoff=",
                                             "multi_files=",
                                             "plot_name=",
                                             "x_file=",
                                             "pca",
                                             "nbin=",
                                             "size=",
                                             "move_average=",
                                             "mean_value=",
                                             "histogram=",
                                             "help"])
except getopt.GetoptError:
    print('Version: 1.7 \
          \nCurrently it works for rmsd, rmsf, sasa_time, sasa_residue, gyration, dipole movement, rdf, distance, PCA  \
          \nIt can read one, two or three same type files \
          \n-o output_name.png, suitable file format for output: png, jpg, jpeg, webp, svg, pdf, eps, json \
          \n-r true: default is fault, renumber the residues, mainly for the duplicated residue numbers while work on rmsf and sasa per residue \
          \n-c number: rdf distance cutoff value \
          \n-a 2 or 3 or true: default is fault, output the average score for each replica \
          \n-p plot_name: the title shoed in plot \
          \n-j whether reading cluster_PCA0.xvg cluster_PCA1.xvg ..., default is 0, if you want generate PCA contour heat map, then type 1 \
          \n-n represent "nbin", mainly used for pca and/or histogram ploting, default value=360, the larger the smoother, however, if there is white line in your PCA_Density plot, please change the value to get a cleaner plot \
          \n-z represent "size", mainly used for pca ploting, default value=500, the larger the higher resolution, if there is white line in your PCA Density plot, please change the value to get a cleaner plot \
          \n-b the window size of moving average analysis for SASA \
          \n-d whether user want to draw a straight line to represent the mean value for each trace. default value = false; user can change it to true \
          \n-i whether user want to generate histogram plot, default value = false \
          \nUsage: \
          \n ./plotly_go -f <file1> -s <file2> -t <file3> -o <output_name> -r <true/fault> -a <2/3/true> -x <xaxis_name> -y <yaxis_name> -c <rdf_cutoff> -p <plot_title> \
          \n ./plotly_go -m <file1> <file2> <file3>  -o <output_name> -r <true/fault> -a <2/3/true> -x <xaxis_name> -y <yaxis_name> -c <rdf_cutoff> -p <plot_title>  \
          \n ./plotly_go -f rmsd1.xvg -o rmsd1.png \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -o rmsd12.png \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -t rmsd3.xvg -o rmsd123.png \
          \n ./plotly_go -f rmsf1.xvg -s rmsf2.xvg -t rmsf3.xvg -o rmsf123.png -r true \
          \n ./plotly_go -f rmsf1.xvg -s rmsf2.xvg -t rmsf3.xvg -o rmsf123.png -a 3 \
          \n ./plotly_go -f rmsf1.xvg -s rmsf2.xvg -o rmsf12.png -a 2 \
          \n ./plotly_go -f sasa1.xvg -s sasa2.xvg -t sasa3.xvg -o sasa123.png -r true -x "Time (ns)" -y "SASA (nm<sup>2</sup>)" \
          \n ./plotly_go -f rdf1.xvg -s rdf2.xvg -t rdf3.xvg -o rdf123.png -c 4.7 -o rdf123.png \
          \n ./plotly_go -m sasa1.xvg sasa2.xvg sasa3.xvg sasa4.xvg sasa5.xvg sasa6.xvg sasa7.xvg sasa8.xvg sasa9.xvg -o multy_sasa_9.png -a true -b 10 -d true\
          \n ./plotly_go -m resarea1.xvg resarea2.xvg resarea3.xvg -o test_resare123.png -a true -r true \
          \n ./plotly_go -m 2dproj1.xvg -n 1000 -z 500 -o pca_myname.png \
          \n ./plotly_go -m dihedral-ave.xvg -o dihedral-ave.png -n 1 -i true -p "Dihedral F38@C-T39@N-T39@CA-T39@C')
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print('Version: 1.7 \
          \nCurrently it works for rmsd, rmsf, sasa_time, sasa_residue, gyration, dipole movement, rdf, distance, PCA  \
          \nIt can read one, two or three same type files \
          \n-o output_name.png, suitable file format for output: png, jpg, jpeg, webp, svg, pdf, eps, json \
          \n-r true: default is fault, renumber the residues, mainly for the duplicated residue numbers while work on rmsf and sasa per residue\
          \n-a 2 or 3 or true: default is fault, output the average score for each replica \
          \n-c number: rdf distance cutoff value \
          \n-p plot_name: the title shoed in plot \
          \n-j whether reading cluster_PCA0.xvg cluster_PCA1.xvg ..., default is 0, if you want generate PCA contour heat map, then type 1 \
          \n-n represent "nbin", mainly used for pca ploting, default value=1000, the larger the smoother, however, if there is white line in your PCA_Density plot, please change the value(-n 500) to get a cleaner plot \
          \n-z represent "size", mainly used for pca ploting, default value=500, the larger the higher resolution, if there is white line in your PCA Density plot, please change the value to get a cleaner plot \
          \n-b the window size of moving average analysis for SASA \
          \n-d whether user want to draw a straight line to represent the mean value for each trace. default value = false; user can change it to true \
          \n-i whether user want to generate histogram plot, default value = false \
          \nUsage: \
          \n ./plotly_go -f <file1> -s <file2> -t <file3> -o <output_name> -r <true/fault> -a <2/3/true> -x <xaxis_name> -y <yaxis_name> -c <rdf_cutoff> -p <plot_title> \
          \n ./plotly_go -m <file1> <file2> <file3>  -o <output_name> -r <true/fault> -a <2/3/true> -x <xaxis_name> -y <yaxis_name> -c <rdf_cutoff> -p <plot_title> \
          \n ./plotly_go -f rmsd1.xvg -o rmsd1.png \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -o rmsd12.png \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -o rmsd12.png -a 2 \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -t rmsd3.xvg -o rmsd123.png \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -t rmsd3.xvg -o rmsd123.png -a 3\
          \n ./plotly_go -f sasa1.xvg -s sasa2.xvg -t sasa3.xvg -o sasa123.png -r true -x "Time (ns)" -y "SASA (nm<sup>2</sup>)" \
          \n ./plotly_go -f rdf1.xvg -s rdf2.xvg -t rdf3.xvg -o rdf123.png -c 4.7 -o rdf123.png \
          \n ./plotly_go -m sasa1.xvg sasa2.xvg sasa3.xvg sasa4.xvg sasa5.xvg sasa6.xvg sasa7.xvg sasa8.xvg sasa9.xvg -o multy_sasa_9.png -a true -b 10 -d true\
          \n ./plotly_go -m resarea1.xvg resarea2.xvg resarea3.xvg -o test_resare123.png -a true -r true \
          \n ./plotly_go -m 2dproj1.xvg -n 1000 -z 500 -o pca_myname.png \
          \n ./plotly_go -m dihedral-ave.xvg -o dihedral-ave.png -n 1 -i true -p "Dihedral F38@C-T39@N-T39@CA-T39@C')
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
        ave = arg
    elif opt in ("-o", "--output_name"):
        output_name = str(arg)
    elif opt in ("-x", "--xaxis_name"):
        xaxis_name = str(arg)
    elif opt in ("-y", "--yaxis_name"):
        yaxis_name = str(arg)
    elif opt in ("-c", "--rdf_cutoff"):
        rdf_cutoff = round(float(arg),1)
    elif opt in ("-l", "--multi_files"):
        multi_files = arg.split(',')
    elif opt in ("-p", "--plot_name"):
        plot_name = str(arg)
    # how to recognize the space seperated command line file inputs
    elif opt in ("-j", "--pca"):
        pca = int(arg)
    elif opt in ("-n", "--nbin"):
        nbin = int(arg)
    elif opt in ("-z", "--size"):
        size = int(arg)
    elif opt in ("-b", "--move_average"):
        move_average = int(arg)
    elif opt in ("-d", "--mean_value"):
        mean_value = str(arg)
    elif opt in ("-i", "--histogram"):
        histogram = str(arg)
    elif opt in ("-m", "--multi_files"):
        value = 1
        multi_files = []
        index = args.index(opt)
        for number in range(index, len(args)):
            if number+1 < len(args):
                if not args[number+1].startswith("-") and value <= len(args):
                    multi_files.append(args[number+1])
                    value += 1
                else:
                    value = value + 100000
        args.remove('-m')
        for value in multi_files:
            args.remove(value)
        try:
            opts, args_lala = getopt.getopt(args,"f:s:t:r:a:o:x:y:c:m:p:l:j:n:z:b:d:i:h",["file1=",
                                                     "file2=",
                                                     "file3=",
                                                     "renumber=",
                                                     "average=",
                                                     "output_name=",
                                                     "xaxis_name=",
                                                     "yaxis_name=",
                                                     "rdf_cutoff=",
                                                     "multi_files=",
                                                     "plot_name=",
                                                     "x_file=",
                                                     "pca",
                                                     "nbin=",
                                                     "size=",
                                                     "move_average=",
                                                     "mean_value=",
                                                     "histogram=",
                                                     "help"])
        except getopt.GetoptError:
            sys.exit(2)

        for opt, arg in opts:
            if opt in ("-f", "--file1"):
                file1 = str(arg)
            elif opt in ("-s", "--file2"):
                file2 = str(arg)
            elif opt in ("-t", "--file3"):
                file3 = str(arg)
            elif opt in ("-r", "--renumber"):
                renumber = str(arg)
            elif opt in ("-a", "--average"):
                ave = arg
            elif opt in ("-o", "--output_name"):
                output_name = str(arg)
            elif opt in ("-x", "--xaxis_name"):
                xaxis_name = str(arg)
            elif opt in ("-y", "--yaxis_name"):
                yaxis_name = str(arg)
            elif opt in ("-c", "--rdf_cutoff"):
                rdf_cutoff = round(float(arg),1)
            elif opt in ("-l", "--multi_files"):
                multi_files = arg.split(',')
            elif opt in ("-p", "--plot_name"):
                plot_name = str(arg)
            elif opt in ("-j", "--pca"):
                pca = int(arg)
            elif opt in ("-n", "--nbin"):
                nbin = int(arg)
            elif opt in ("-z", "--size"):
                size = int(arg)
            elif opt in ("-b", "--move_average"):
                move_average = int(arg)
            elif opt in ("-d", "--mean_value"):
                mean_value = str(arg)
            elif opt in ("-i", "--histogram"):
                histogram = str(arg)





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

    def __init__(self,file1, file2, file3, output_name, renumber, ave, xaxis_name, yaxis_name, rdf_cutoff, multi_files, plot_name, pca, nbin, size, move_average, mean_value, histogram):
        # if multi_files == 0:
        #     self.flag_recognizer(file1, file2, file3)
        #     if self.flag == 'rmsd':
        #         self.rmsd_averager(file1,file2,file3, ave, output_name, xaxis_name, yaxis_name)
        #     elif self.flag == 'rmsf':
        #         self.rmsf_averager(file1, file2,file3, renumber, ave, output_name, xaxis_name, yaxis_name)
        #     elif self.flag == 'sasa' and self.sasa_flag != '-or':
        #         self.sasa_averager(file1, file2, file3, ave, output_name, xaxis_name, yaxis_name)
        #     elif self.flag == 'sasa' and self.sasa_flag == '-or':
        #         self.sasa_res_averager(file1, file2, file3,renumber,ave, output_name, xaxis_name, yaxis_name)
        #     elif self.flag == 'gyrate':
        #         self.gyrate_averager(file1, file2,file3, ave, output_name, xaxis_name, yaxis_name)
        #     elif self.flag == 'dipoles':
        #         self.dipoles_averager(file1, file2,file3, ave, output_name, xaxis_name, yaxis_name)
        #     elif self.flag == 'distance':
        #         self.distance_averager(file1, file2,file3, ave, output_name, xaxis_name, yaxis_name)
        #     elif self.flag == 'rdf':
        #         self.rdf_averager(file1,file2,file3, ave, output_name, xaxis_name, yaxis_name, rdf_cutoff)



        if len(multi_files) >=1:
            # print(multi_files)
            file1 = multi_files[0]
            self.flag_recognizer(file1, file2, file3)
            if self.pca_flag != 1 and self.flag != 'pca':
                self.plotly_multy(self.flag, multi_files, xaxis_name, yaxis_name, renumber, ave, output_name, rdf_cutoff, plot_name, move_average, mean_value, histogram, nbin)
            elif self.pca_flag == 1:
                self.plotly_pca(multi_files, xaxis_name, yaxis_name, renumber, ave, output_name, rdf_cutoff, plot_name, nbin, size, move_average)
            elif self.flag == 'pca':
                self.plotly_pca(multi_files, xaxis_name, yaxis_name, renumber, ave, output_name, rdf_cutoff, plot_name, nbin, size, move_average)





    def flag_recognizer(self,file1, file2, file3):                                                   # first method to be called in __main__, used for creating object and charactors.
        #self.atomic_index.clear()                                                   # cleaveage the information of previous object before put new record into these charactors
        with open(file1, 'r') as f:
            lines = f.readlines()                                                 # "filename" = $PATH/crankpep_docking_results/PLDAYL_corrected_top_1.pdb
            if len(lines)  >= 3:
                try:
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
                    elif self.flag == 'rdf,':
                        self.flag = 'rdf'
                    elif self.flag == 'convergence':
                        self.flag = 'convergence'
                    elif self.flag == 'anaeig,':
                        self.flag = 'pca'
                    elif self.flag == 'angle,':
                        self.flag = 'angle'
                except:
                    pass

            if len(lines) >= 9 and '-or' in lines[8]:
                self.sasa_flag = '-or'

            if 'pca' in str(file1).lower() or '2dproj' in str(file1):
                self.pca_flag = 1
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


    def plotly_multy(self, flag, multi_files, xaxis_name, yaxis_name, renumber, ave, output_file_name, rdf_cutoff, plot_name, move_average, mean_value, histogram, nbin):
        Plotly = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
        a=1
        data = []
        histogram_data = []
        group_labels = []

        ################## for regular expression substance ##################
        regex = r"\[|\]|'"

        ################## define plot title, x axis name and y axis name ##################
        if plot_name == '':
            with open(multi_files[0], "r") as f:
                plot_title = re.sub(regex, "", str(re.findall('"([^"]*)"', f.readlines()[13])))
        else:
            plot_title = str(plot_name)
        if xaxis_name == 0:
            with open(multi_files[0], "r") as f:
                x_name = re.sub(regex, "", str(re.findall('"([^"]*)"', f.readlines()[14])))
        else:
            x_name = xaxis_name
            pass

        if yaxis_name == 0 and plot_title not in ['Solvent Accessible Surface', 'Area per residue over the trajectory']:
            with open(multi_files[0], "r") as f:
                y_name = re.sub(regex, "", str(re.findall('"([^"]*)"', f.readlines()[15])))
        elif yaxis_name == 0 and plot_title in ['Solvent Accessible Surface', 'Area per residue over the trajectory']:
            y_name = 'Area (nm<sup>2</sup>)'
        elif yaxis_name == 0 and self.flag == 'dihedral_distribution' and histogram == 'true':
            yname = 'Probability'
        else:
            y_name = yaxis_name
            pass
        ################## reading the datas!! ##################
        for i in multi_files:
            # create empty list
            locals()["x_" + str(a)] = []
            locals()["y_" + str(a)] = []
            locals()["sd_" + str(a)] = []
            # grab datas from input files
            with open(i, "r") as f:
                lines = f.readlines()
                for num in range(len(lines)):
                    if lines[num].startswith("#") or lines[num].startswith("@"):
                        pass
                    else:
                        if x_name == 'Time (ps)':       # set time value into ns
                            locals()["x_" + str(a)].append(float(lines[num].split()[0])/1000)
                        else:
                            locals()["x_" + str(a)].append(float(lines[num].split()[0]))
                        locals()["y_" + str(a)].append(float(lines[num].split()[1]))
                        try:
                            locals()["sd_" + str(a)].append(float(lines[num].split()[2]))
                        except:
                            pass
            if x_name == 'Residue' and renumber == 'true':
                for k in range(len(locals()["x_" + str(a)])):
                    locals()["x_" + str(a)][k] = k+1
            # find out if the residues indexs are consited        
            if x_name == 'Residue':
                seq1, seq2 = self.consist(locals()["x_" + str(a)])
                if seq2 == []:
                    ################## define traces ##################
                    locals()["trace" + str(a)] = go.Scatter(x=locals()["x_" + str(a)], y=locals()["y_" + str(a)], line=dict(color=Plotly[a-1]), name=str(i).split('.')[0])
                    data.append(locals()["trace" + str(a)])
                    # print(Plotly[a-1])
                else:
                    locals()["trace1" + str(a)] = go.Scatter(x=seq1, y=locals()["y_" + str(a)][:len(seq1)], line=dict(color=Plotly[a-1]), name=str(i).split('.')[0])
                    locals()["trace2" + str(a)] = go.Scatter(x=seq2, y=locals()["y_" + str(a)][-len(seq2):], line=dict(color=Plotly[a-1]), showlegend=False)
                    data.append(locals()["trace1" + str(a)])
                    data.append(locals()["trace2" + str(a)])
            else:
                locals()["trace" + str(a)] = go.Scatter(x=locals()["x_" + str(a)], y=locals()["y_" + str(a)], line=dict(color=Plotly[a-1]), name=str(i).split('.')[0])
                data.append(locals()["trace" + str(a)])
            
            # add mean value straight line
            if mean_value == "true":
                mean_value_number = np.mean(locals()["y_" + str(a)])
                my_list = [mean_value_number] * len(locals()["x_" + str(a)])
                locals()["trace_mean_value" + str(a)] = go.Scatter(x=locals()["x_" + str(a)], y=my_list, line=dict(color=Plotly[a-1], dash='dash'), name='mean_value', showlegend=False)
                data.append(locals()["trace_mean_value" + str(a)])
            else:
                pass
            ######## histogram distribution plot ########
            if histogram == 'true':
            
                # Group data together
                histogram_data.append(locals()["y_" + str(a)])
                # Legend names
                group_labels.append(str(i).split('.')[0])
                

            
            a += 1
        ################## test if time unit is ns, if not then change it from ps to ns ##################
        if x_name == 'Time (ps)':
            x_name = 'Time (ns)'


        ################## plot the datas ##################
        # if mean_value == "true":
        #     layout = go.Layout(title=plot_title, title_x=0.5, title_y=1, font=dict(size=24),
        #                        xaxis=dict(title=x_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
        #                                   showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
        #                        yaxis=dict(title=y_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
        #                                   showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
        #                        legend=dict(x=1, y=1, orientation='v', font=dict(size=30)), showlegend=True,
        #                        plot_bgcolor='rgba(255, 255, 255, 0.1)',
        #                        paper_bgcolor='rgba(255, 255, 255, 0.2)',
        #                        width=800, height=600,
        #                        shapes=[{'type': 'line',
        #                                 'x0': 0,
        #                                 'y0': mean_value,
        #                                 'x1': max(data),
        #                                 'y1': mean_value,
        #                                 'line': {'color': 'black', 'width': 2, 'dash': 'dashdot'}}])
        # else:
        layout = go.Layout(title=plot_title, title_x=0.5, title_y=1, font=dict(size=24),
                           xaxis=dict(title=x_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                      showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                           yaxis=dict(title=y_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                      showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                           legend=dict(x=1, y=1, orientation='v', font=dict(size=30)), showlegend=True,
                           plot_bgcolor='rgba(255, 255, 255, 0.1)',
                           paper_bgcolor='rgba(255, 255, 255, 0.2)',
                           width=800, height=600)
            
        fig = go.Figure(data=data, layout=layout)
        pio.write_image(fig, output_file_name)
        #################### histogram distribution plotting #####################
        if histogram == 'true':
            # give an output file name
            histogram_output_filename = "histogram_" + output_file_name
            # generate figures
            fig_hist = ff.create_distplot(histogram_data, group_labels, bin_size=nbin, colors=Plotly, curve_type='normal')
            # Create distplot with curve_type set to 'normal'
            # fig = ff.create_distplot([x1, x2], group_labels, bin_size=.5, curve_type='normal', colors=colors)
            # Add title
            # fig_hist.update(layout_title_text="lalala")
            fig_hist.update_layout(title_text=plot_title, titlefont=dict(size=24),
                                   title_x=0.5, xaxis_range=[-180, 180], 
                                   title_font=dict(size=24),  # Set title font size
                                   xaxis=dict(
                                        title="Degree",  # x-axis label
                                        title_font=dict(size=20),  # x-axis label font size
                                        range=[-180, 180],  # Set the range for the x-axis
                                        tickvals=[-180,-90, 0, 90, 180],  # Set ticks at -180 and +180
                                        tickfont=dict(size=24)
                                    ),
                                   yaxis=dict(
                                        title="Probability",  # y-axis label
                                        title_font=dict(size=24),  # y-axis label font size
                                        tickfont=dict(size=24)
                                    ),
                                   legend=dict(font=dict(size=24)),  # Set legend font size
                                   plot_bgcolor='rgba(255, 255, 255, 0.1)', 
                                   paper_bgcolor='rgba(255, 255, 255, 0.2)',
                                   width=800, height=600)
            pio.write_image(fig_hist, histogram_output_filename)
            
            
        ################## if user ask for average the inputs ##################
        if ave == 'true':
            number = len(multi_files)
            # average_value = [(locals()["x_" + str(a)][i] + values2[i] + values3[i]) / 3 for i in range(len(values1))]
            # average_sd = [(sd1[i] + sd2[i] + sd3[i]) / 3 for i in range(len(sd1))]
            average_value = locals()["y_" + str(1)]
            for a in range(1, number):
                average_value = [x + y for x, y in zip(average_value, locals()["y_" + str(a+1)])]
            average_value = [x/number for x in average_value]

            while len(locals()["x_" + str(a)]) != len(average_value):
                if len(locals()["x_" + str(a)]) > len(average_value):
                    locals()["x_" + str(a)].pop()
                else:
                    average_value.pop()
            data = []
            # trace_ave = go.Scatter(x=locals()["x_" + str(a)], y=average_value, name='Average Values')
            # data.append(trace_ave)
            
            # find out if the residues indexs are consited        
            if x_name == 'Residue':
                seq1, seq2 = self.consist(locals()["x_" + str(a)])
                if seq2 == []:
                    ################## define traces ##################
                    trace_ave = go.Scatter(x=locals()["x_" + str(a)], y=average_value, name='Average Values')
                    data.append(trace_ave)
                else:
                    trace_ave1 = go.Scatter(x=seq1, y=average_value[:len(seq1)], line=dict(color=Plotly[0]), name='Average Values')
                    trace_ave2 = go.Scatter(x=seq2, y=average_value[-len(seq2):], line=dict(color=Plotly[0]), showlegend=False)
                    data.append(trace_ave1)
                    data.append(trace_ave2)
            else:
                trace_ave = go.Scatter(x=locals()["x_" + str(a)], y=average_value, name='Average Values')
                data.append(trace_ave)
            if mean_value == "true":
                mean_value_number = np.mean(average_value)
                my_list = [mean_value_number] * len(locals()["x_" + str(a)])
                trace_mean_value = go.Scatter(x=locals()["x_" + str(a)], y=my_list, line=dict(color='black', dash='dash'), name='mean_value', showlegend=False)
                data.append(locals()["trace_mean_value" + str(a)])
            else:
                pass

            
            layout = go.Layout(title=plot_title, title_x=0.5, title_y=1, font=dict(size=24),
                               xaxis=dict(title=x_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               yaxis=dict(title=y_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               legend=dict(x=1, y=1, orientation='v', font=dict(size=30)), showlegend=True,
                               plot_bgcolor='rgba(255, 255, 255, 0.1)',
                               paper_bgcolor='rgba(255, 255, 255, 0.2)',
                               width=800, height=600)
            fig = go.Figure(data=data, layout=layout)
            pio.write_image(fig, "Average_" + output_file_name)
            ######### output the average.xvg file ########
            with open(output_file_name[:-4]+"_average.xvg", 'w') as f:
                with open(multi_files[0], "r") as a:
                    lines = a.readlines()
                    for num in range(len(lines)):
                        if lines[num].startswith("#") or lines[num].startswith("@"):
                            f.write(lines[num])
                        else:
                            pass
                for num in range(len(average_value)):
                    average_line = "{}     {}\n".format(locals()["x_" + str(1)][num], average_value[num])
                    f.write(average_line)
                    
        ################## if user ask for movement average ##################
        if move_average != 0:
            ma_plot_title = plot_title + "<br><sup>Moving Average (window size=" + str(move_average) + ")</sup>"
            ma_output_file_name = "MovingAverage_" + output_file_name
            data = []
            # 窗口大小
            window_size = move_average   
            # count the number of input files
            number = len(multi_files)
            a = 1
            for i in multi_files:
                try:
                    # 计算移动平均线
                    moving_average = np.convolve(locals()["y_" + str(a)], np.ones(window_size) / window_size, mode='valid')
                    locals()["moving_average_" + str(a)] =  moving_average
                    # print(len(locals()["moving_average_" + str(a)]))  # 打印属性值
                    # print(locals()["moving_average_" + str(a)]) # 打印属性值

                    # set the traces
                    locals()["average_move_trace_" + str(a)] = go.Scatter(x=locals()["x_" + str(a)][window_size-1:], y=locals()["moving_average_" + str(a)], line=dict(color=Plotly[a-1]), name=str(i).split('.')[0])
                    data.append(locals()["average_move_trace_" + str(a)])
                    if mean_value == "true":
                        mean_value_number = np.mean(locals()["moving_average_" + str(a)])
                        my_list = [mean_value_number] * len(locals()["x_" + str(a)])
                        trace_mean_value = go.Scatter(x=locals()["x_" + str(a)][window_size-1:], y=my_list, line=dict(color=Plotly[a-1], dash='dash'), name='mean_value', showlegend=False)
                        data.append(locals()["trace_mean_value" + str(a)])
                    else:
                        pass
                    a += 1
                except IndexError:
                    pass    
            layout = go.Layout(title=ma_plot_title, title_x=0.5, title_y=0.95, font=dict(size=24),
                               xaxis=dict(title=x_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               yaxis=dict(title=y_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               legend=dict(x=1, y=1, orientation='v', font=dict(size=30)), showlegend=True,
                               plot_bgcolor='rgba(255, 255, 255, 0.1)',
                               paper_bgcolor='rgba(255, 255, 255, 0.2)',
                               width=800, height=600)
            fig = go.Figure(data=data, layout=layout)
            pio.write_image(fig, ma_output_file_name)
            
                          

    def plotly_pca(self, multi_files, xaxis_name, yaxis_name, renumber, ave, output_name, rdf_cutoff, plot_name, nbin, size, move_average):
        x_points = []
        y_points = []

        ################## define plot title, x axis name and y axis name ##################
        if plot_name == '':
            plot_title = 'PCA 2D projection of trajectory'
        else:
            plot_title = str(plot_name)
            pass

        if xaxis_name == 0:
            x_name = 'projection on eigenvector 1 (nm)'
        else:
            x_name = xaxis_name
            pass

        if yaxis_name == 0 and plot_title not in ['Solvent Accessible Surface', 'Area per residue over the trajectory']:
            y_name = 'projection on eigenvector 2 (nm)'
        else:
            y_name = yaxis_name
            pass
        ################## reading the datas!! ##################
        for i in multi_files:
            x_points=[]
            y_points=[]
            labels = []
            # create empty list
            with open(i, "r") as f:
                lines = f.readlines()
                for num in range(len(lines)):
                    if lines[num].startswith("#") or lines[num].startswith("@"):
                        pass
                    else:
                        x_points.append(float(lines[num].split()[0]))
                        y_points.append(float(lines[num].split()[1]))
                        labels.append(num)

            # # 创建散点图轨迹
            # scatter_trace = go.Scatter(
            #     x=x_points,
            #     y=y_points,
            #     mode='markers',
            #     marker=dict(
            #         size=5,
            #         color='black'
            #     )
            # )

            # layout = go.Layout(title=plot_title, title_x=0.5, title_y=1, font=dict(size=20),
            #                     xaxis=dict(title=x_name, titlefont=dict(size=20, color='black', family='Arial'), zeroline=False, autorange=True,
            #                               showgrid=False, gridwidth=1, gridcolor='rgba(0,0,0,0.1)', tickfont=dict(size=20)),
            #                     yaxis=dict(title=y_name, titlefont=dict(size=20, color='black', family='Arial'), zeroline=False, autorange=True,
            #                               showgrid=False, gridwidth=1, gridcolor='rgba(0,0,0,0.1)', tickfont=dict(size=20)),
            #                     legend=dict(x=1, y=1, orientation='v', font=dict(size=30)), showlegend=False,
            #                     plot_bgcolor='rgba(255, 255, 255, 0.1)',
            #                     paper_bgcolor='rgba(255, 255, 255, 0.2)')

            # # 创建图形对象
            # fig = go.Figure(data=scatter_trace, layout=layout)

            # # 显示图形
            # # fig.show()
            # if output_name == 'plotly.png':
            #     pio.write_image(fig, "PCA_Scatter_"+i.split('.')[0]+".png")
            # else:
            #     pio.write_image(fig, "Scatter_" + output_name)

            # 创建一个 Scatter 对象
            scatter = go.Scatter(
                x=x_points,
                y=y_points,
                mode='markers',
                marker=dict(
                    color=labels,  # 设置颜色为标签的数值
                    colorscale='rainbow',  # 颜色映射，你可以根据需要选择不同的颜色映射
                    colorbar=dict(title='Label Range'),  # 添加颜色条
                ),
            )

            # 创建数据列表
            data = [scatter]

            # 创建布局
            layout = go.Layout(
                title='PCA plot with Color Bar for frame order', title_x=0.5, title_y=1, font=dict(size=24),
                xaxis=dict(title='PC1 (nm)', titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                           showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                yaxis=dict(title='PC2 (nm)', titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                           showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                plot_bgcolor='rgba(255, 255, 255, 0.1)',
                paper_bgcolor='rgba(255, 255, 255, 0.2)',
                width=800, height=600
            )

            # 创建 Figure 对象
            fig = go.Figure(data=data, layout=layout)
            if output_name == 'plotly.png':
                pio.write_image(fig, "PCA_Scatter_"+i.split('.')[0]+".png")
            else:
                pio.write_image(fig, "Scatter_" + output_name)

####################################################################################################################################
        # 用R创建PCA图
        x_diff = np.max(x_points) - np.min(x_points)
        y_diff = np.max(y_points) - np.min(y_points)
        width= (x_diff / y_diff) * size
        height = 1 * size
        # 读取数据
        for i in multi_files:
            r('data <- read.table(%r, skip = 17, header = FALSE)' % (i))

            # 设置变量名
            r('names(data) <- c("PC1", "PC2")')

            # 导入所需的R包
            r('library(RColorBrewer)')

            # 执行其他指令
            r('zBot <- 1.52')
            r('zTop <- 3.42')
            r('zW <- 0.83')
            r('buylrd <- rev(brewer.pal(11,"RdYlBu"))')

            # 保存为PNG文件
            if output_name == 'plotly.png':
                # r('png(file=%r, height=600, width=450)' % ("PCA_Density_" + i.split('.')[0] + ".png"))
                # r('png(file=%r)' % ("PCA_Density_" + i.split('.')[0] + ".png"))
                r('png(file=%r, height=%r, width=%r)' % ("PCA_Density_" + i.split('.')[0] + ".png", height, width))
                r('smoothScatter(data$PC2 ~ data$PC1, nbin=%r, colramp = colorRampPalette(c(buylrd)),nrpoints=Inf, pch="", cex=.7,col="black",main=%r, xlab=%r, ylab=%r,transformation = function(x) x^.45)' % (nbin, plot_title, x_name, y_name))
                r('dev.off()')
            else:
                # r('png(file=%r, height=600, width=450)' % (output_name))
                # r('png(file=%r)' % (output_name))
                r('png(file=%r, height=%r, width=%r)' % ("Density_" + output_name, height, width))
                r('smoothScatter(data$PC2 ~ data$PC1, nbin=%r, colramp = colorRampPalette(c(buylrd)),nrpoints=Inf, pch="", cex=.7,col="black",main=%r, xlab=%r, ylab=%r,transformation = function(x) x^.45)' % (nbin, plot_title, x_name, y_name))
                r('dev.off()')

            # # 保存为PDF文件
            # r('pdf(file="PCA.pdf", height=1600, width=1600,paper = "a4")')
            # r('smoothScatter(data$PC2 ~ data$PC1, nbin=1000, colramp = colorRampPalette(c(buylrd)),nrpoints=Inf, pch="", cex=.7,col="black",main="Shangze MD simulation PCA analysis", xlab="PCA1", ylab="PCA2",transformation = function(x) x^.5)')
            # r('dev.off()')








x = plotly_go(file1,file2,file3, output_name, renumber, ave, xaxis_name, yaxis_name, rdf_cutoff, multi_files, plot_name, pca, nbin, size, move_average, mean_value, histogram)
