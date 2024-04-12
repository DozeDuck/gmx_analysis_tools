import MDAnalysis as mda
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
# parser大法
import argparse

# Define command line arguments
parser = argparse.ArgumentParser(description='dssp ploting. \n'
                                 'Usage: python plot_dssp.py -f dssp.dat -t md_noPBC_dt1000.pdb -o dssp.png',
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-f', '--dat', required=True, default='dssp.dat', help='input dssp.dat file')
parser.add_argument('-t', '--traj', required=True, default='md_noPBC_dt1000.pdb', help='the traj  file')
parser.add_argument('-o', '--output', required=False, default='dssp.png', help='output name')

args = parser.parse_args()

# define the parameters
dssp = args.dat
traj = args.traj
output_name = args.output
#####################################################################################################

# 载入PDB文件
u = mda.Universe(traj)

# 提取时间戳和残基列表
# times = [ts.time / 1000 for ts in u.trajectory]  # 假设时间单位是ps，转换为ns 
times = [ts.time for ts in u.trajectory]  # 假设时间单位是ns 
residues = [res.resname + str(res.resid) for res in u.residues]

# 读取DSSP数据（dssp.dat）
with open(dssp, "r") as file:
    dssp_lines = file.readlines()

# 转换DSSP数据为列表
dssp_data = [list(line.strip()) for line in dssp_lines]

# 创建DataFrame
df = pd.DataFrame(data=dssp_data, index=times, columns=residues)

# 定义转换字典
structure_values_original = {
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

# 使用字典转换DataFrame中的值
# df.replace(structure_values, inplace=True)
df.replace(structure_values_original, inplace=True)


#######################################################################
# times = times.append(times[-1]+1) # for scale the color
# # Calculate the length of each third
# third_length = len(df.columns) // 3
# 
# # Create a new row
# new_row = [1] * third_length + [2] * third_length + [3] * (len(df.columns) - 2 * third_length)
# 
# # Append the new row at the bottom of the DataFrame
# df.loc[len(df)] = new_row
#######################################################################

# Define color scale and color bar settings
colorscale = [[0, 'red'], [0.5, 'green'], [1.0, 'rgb(0, 0, 255)']]
colorscale = [[0, 'gold'], [0.5, 'mediumturquoise'], [1, 'lightsalmon']]
colorscale = [[0.00, "gold"],   [0.33, "gold"], [0.33, "mediumturquoise"], [0.66, "mediumturquoise"], [0.66, "lightsalmon"],  [1.00, "lightsalmon"]]
# colorscale = [[0.00, "red"],   [0.33, "red"], [0.33, "green"], [0.66, "green"], [0.66, "blue"],  [1.00, "blue"]]
# 将颜色条的刻度设置为固定的值
colorbar_tickvals = [1, 2, 3]
colorbar_ticktext = ['red', 'green', 'blue']
#
colorbar_ticks = [1, 2, 3]  # The actual values in your data
colorbar_ticktext = ['loop', 'helix', 'beta sheet']  # Descriptive labels

# 创建热图
fig = go.Figure(data=go.Heatmap(
    z=df.T.values,
    x=df.index,
    y=df.columns,
    # colorscale=colorscale,
    # colorbar=dict(title = 'type', titleside = 'top', tickmode = 'array', tickvals=colorbar_tickvals, ticktext=colorbar_ticktext),
    hoverongaps=False))

fig.update_layout(
    title='Secondary Structure Analysis Over Time',
    yaxis_title='Residue',
    xaxis_title='Time (ns)',
    width=800,
    height=600)

# 显示热图
pio.write_image(fig, output_name)
