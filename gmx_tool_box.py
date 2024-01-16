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
# from rpy2.robjects import r
# for histogram dist plot
import plotly.figure_factory as ff
import argparse
# for bool values
import ast
# metal restraints adding
import math
          
parser = argparse.ArgumentParser(description='Version: 1.0  \n'
                                 'Usage: python gmx_tool_box.py can generate figures based on *.xvg files; adding restraints for metals; merge lig_GMX.gro to rec.gro and modify topol.top automaticlly \n' \
                                  './gmx_tool_box -m <file1> <file2> <file3>  -o <output_name> -r <true/fault> -a <2/3/true> -x <xaxis_name> -y <yaxis_name> -c <rdf_cutoff> -p <plot_title> \n ' \
                                  './gmx_tool_box -m rmsd1.xvg rmsd2.xvg rmsd3.xvg -o rmsd123.png -a true \n '\
                                  './gmx_tool_box -m rmsf1.xvg rmsf2.xvg rmsf3.xvg -rn true -o rmsd123.png -a true \n '\
                                  './gmx_tool_box -m sasa1.xvg sasa2.xvg sasa3.xvg -o sasa123.png -rn true -xn "Time (ns)" -yn "SASA (nm<sup>2</sup>)" \n ' \
                                  './gmx_tool_box -m rdf1.xvg rdf2.xvg rdf3.xvg -rc 4.7 -o rdf123.png \n ' \
                                  './gmx_tool_box -m sasa1.xvg sasa2.xvg sasa3.xvg sasa4.xvg sasa5.xvg sasa6.xvg sasa7.xvg sasa8.xvg sasa9.xvg -o multy_sasa_9.png -a true -ma 10 -ml true \n '\
                                  './gmx_tool_box -m resarea1.xvg resarea2.xvg resarea3.xvg -o test_resare123.png -a true -r true \n ' \
                                  './gmx_tool_box -m 2dproj1.xvg -pn 1000 -pz 500 -o pca_myname.png \n ' \
                                  './gmx_tool_box -m dihedral-ave.xvg -o dihedral-ave.png -pn 1 -hg true -t "Dihedral F38@C-T39@N-T39@CA-T39@C" \n ' \
                                  '##################### Adding restraints to metal atoms ##################### \n' \
                                  './plotly_go -mtgro rec.gro -mtml MG MN ZN CA -mtrl HIS GLU ASP ASN CYS LYS TYR -mtal ND1 OE1 OE2 OD1 OD2 ND2 SG NZ OH -mtd 0.4 -mtn 3 -mtbr 200000 -mtar 10000 \n ' \
                                  '##################### Merge receptor.gro and ligand_GMX.gro and edit topol.top ##################### \n' \
                                  './plotly_go -gmgro rec.gro -gmlgro lig_GMX.gro -gmlitp lig_GMX.itp -gmtop topol.top',
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
parser.add_argument('-tf', '--title_font', default=24, help='User indicate the size of plot title font, default = 24')
parser.add_argument('-ls', '--legend_show', default='True', help='Whether user want to show the legend, default = True')
parser.add_argument('-lf', '--legend_font', default=30, help='Define the legend font size, defaul = 30')
parser.add_argument('-fm', '--font_family', default='Arial', help='Font family default = Arial')
parser.add_argument('-fc', '--font_color', default='black', help='Font color default = black')
parser.add_argument('-gs', '--grid_show', default='True', help='Whether show the grid line, defaul = True')
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
xaxis_size      = int(args.xaxis_size)
yaxis_size      = int(args.yaxis_size)
xy_font         = int(args.xy_font)
title_font      = int(args.title_font)
legend_show     = ast.literal_eval(args.legend_show)
legend_font     = int(args.legend_font)
font_family     = str(args.font_family)
font_color      = str(args.font_color)
grid_show       = ast.literal_eval(args.grid_show)
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

    def __init__(self, multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show):

        if len(multi_files) >=1:
            # print(multi_files)
            file1 = multi_files[0]
            self.flag_recognizer(file1)
            if self.pca_flag != 1 and self.flag != 'pca':
                self.plotly_multy(multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, self.flag)
            elif self.pca_flag == 1:
                self.plotly_pca(multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, self.flag)
            elif self.flag == 'pca':
                self.plotly_pca(multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, self.flag)





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
            'angle,': 'angle'
        }
        
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

    def define_trace(self, x_data, y_data, file_name, colour, flag=0, labels=0):
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
        else:
            trace = go.Scatter(x=x_data, y=y_data, line=dict(color=colour), name=str(file_name).split('.')[0])
        return trace

    def setup_layout(self, plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show, flag=0):
        # 设置布局
        if flag == 'pca':
            legend_show = False
        layout = go.Layout(
            title=plot_title, title_x=0.5, title_y=1, font=dict(size=title_font, color=font_color),
            xaxis=dict(title=x_name, titlefont=dict(size=xy_font, color=font_color, family=font_family), zeroline=False, autorange=True,
                       showgrid=grid_show, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
            yaxis=dict(title=y_name, titlefont=dict(size=xy_font, color=font_color, family=font_family), zeroline=False, autorange=True,
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


    def plotly_multy(self, multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, flag):
        Plotly = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
        data, histogram_data, group_labels = [], [], []

        # 读取plot_title, x_name, y_name
        plot_title, x_name, y_name = self.extract_plot_details(multi_files, plot_name, xaxis_name, yaxis_name, flag, histogram)

        # 读取数据并创建迹线
        for i, file in enumerate(multi_files):
            x_data, y_data, _ = self.read_data(file, xaxis_name, renumber)
            trace = self.define_trace(x_data, y_data, file, Plotly[i % len(Plotly)])
            data.append(trace)

            # 添加直方图数据
            if histogram == 'true':
                histogram_data.append(y_data)
                group_labels.append(str(file).split('.')[0])

        # 设置布局
        layout = self.setup_layout(plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show)

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


        # # 调用上述方法
        # for i, file in enumerate(multi_files):
        #     x_data, y_data, sd_data = self.read_data(file, xaxis_name, renumber)
        #     trace = self.define_trace(x_data, y_data, file, Plotly[i % len(Plotly)])
        #     data.append(trace)

        # layout = self.setup_layout(plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show)
        # self.plot_graph(data, layout, output_name)


         
            
    def plotly_pca(self, multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, flag):
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
        layout = self.setup_layout(plot_title, title_font, 'PC1 (nm)', 'PC2 (nm)', xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show, flag)

        # 使用 plot_graph 绘制图形
        self.plot_graph(data, layout, "Scatter_" + output_name)

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





if multi_files[0] != 0:
    x = plotly_go(multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show)
elif mr_file != 0:
    x = mr(mr_file, mr_num_neighbours, mr_distance_value, mr_atom_list, mr_metal_list, mr_residue_list, mr_bond_strength, mr_angle_strength)
elif gm_receptor_gro != 0:
    x = gromerger(gm_receptor_gro, gm_ligand_gro, gm_ligand_itp, gm_receptor_top )


