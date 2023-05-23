import argparse
import numpy as np
import re
import math

# Define command line arguments
parser = argparse.ArgumentParser(description='Metal atoms restraints adding. \n'
                                 'Usage: python mr_adding_multi.py -f rec.gro -m1 888 -a1 358 832 404 -m2 889 -a2 157 524 128',
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-f', '--gro', nargs='+', required=True, default=0, help='input gro file')
parser.add_argument('-m1', '--metal1', required=True, default=0, help='the indexs number for your metal atoms, showed in rec.gro file')
parser.add_argument('-a1', '--atom1', nargs='+', required=True, default=0, help='O/N atom index numbers for each metal atoms, for calculating the b-length, angle')
parser.add_argument('-m2', '--metal2', default=0, help='the indexs number for your metal atoms, showed in rec.gro file')
parser.add_argument('-a2', '--atom2', nargs='+',  default=[0], help='O/N atom index numbers for each metal atoms, for calculating the b-length, angle')
parser.add_argument('-m3', '--metal3', default=0, help='the indexs number for your metal atoms, showed in rec.gro file')
parser.add_argument('-a3', '--atom3', nargs='+',  default=[0], help='O/N atom index numbers for each metal atoms, for calculating the b-length, angle')
parser.add_argument('-m4', '--metal4', default=0, help='the indexs number for your metal atoms, showed in rec.gro file')
parser.add_argument('-a4', '--atom4', nargs='+',  default=[0], help='O/N atom index numbers for each metal atoms, for calculating the b-length, angle')
parser.add_argument('-m5', '--metal5', default=0, help='the indexs number for your metal atoms, showed in rec.gro file')
parser.add_argument('-a5', '--atom5', nargs='+', default=[0], help='O/N atom index numbers for each metal atoms, for calculating the b-length, angle')
parser.add_argument('-m6', '--metal6', default=0, help='the indexs number for your metal atoms, showed in rec.gro file')
parser.add_argument('-a6', '--atom6', nargs='+', default=[0], help='O/N atom index numbers for each metal atoms, for calculating the b-length, angle')
# parser.add_argument('-t', '--topology_file', default='topol.top', help='the path for topol.top')
parser.add_argument('-l', '--bond', default=200000, help='give the path for your em.mdp file, or this script will use the default_em.mdp')
parser.add_argument('-g', '--angle', default=10000, help='give the path for your nvt.mdp file, or this script will use the default_nvt.mdp')

args = parser.parse_args()

# define the parameters
file = args.gro[0]
metal1 = int(args.metal1)
atom1 = [int(x) for x in args.atom1]
metal2 = int(args.metal2)
atom2 = [int(x) for x in args.atom2]
metal3 = int(args.metal3)
atom3 = [int(x) for x in args.atom3]
metal4 = int(args.metal4)
atom4 = [int(x) for x in args.atom4]
metal5 = int(args.metal5)
atom5 = [int(x) for x in args.atom5]
metal6 = int(args.metal6)
atom6 = [int(x) for x in args.atom6]
bond_strength = args.bond
angle_strength = args.angle


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
    last = ''
    
    
    def __init__(self, gro, metal1, atom1, metal2, atom2, metal3, atom3, metal4, atom4, metal5, atom5, metal6, atom6, bond_strength, angle_strength):
        self.GROreader(gro)
        self.bond_cal( metal1, atom1, metal2, atom2, metal3, atom3, metal4, atom4, metal5, atom5, metal6, atom6,bond_strength)
        self.pair_cal( metal1, atom1, metal2, atom2, metal3, atom3, metal4, atom4, metal5, atom5, metal6, atom6)
        self.angle_cal( metal1, atom1, metal2, atom2, metal3, atom3, metal4, atom4, metal5, atom5, metal6, atom6, angle_strength)

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
                self.resname.append(match.group(2))
            self.atomname.append(line.split()[1])                        # The 3rd column is the atom name C CA CD1 CD2 and so on
            self.index.append(int(line.split()[2]))                   # Column 4 is the residue name TYR ALA etc.
            self.x.append(float(line.split()[3]))                         # The 5th column is the name of the chain it is on
            self.y.append(float(line.split()[4]))               # The sixth column is the residue number
            self.z.append(float(line.split()[5]))                   # Column 7 is the x-coordinate of the atom
            
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
                metal = locals()["metal" + str(i)]
                atom = locals()["atom" + str(i)]
                metal_point = [self.x[metal-1],self.y[metal-1],self.z[metal-1]]
                print("; please add below to topol.top's distance part")
                for i in atom:
                    locals()["atom" + str(i)] = [self.x[i-1],self.y[i-1],self.z[i-1]]
                    print("%5d%6d%6d%7.2f%9d" % (metal, i, 6, self.calculate_distance(metal_point, locals()["atom" + str(i)]), bond_strength))
        
    def pair_cal(self,  metal1, atom1, metal2, atom2, metal3, atom3, metal4, atom4, metal5, atom5, metal6, atom6):
        for i in [1,2,3,4,5,6]:
            if locals()["metal" + str(i)] !=0:
                metal = locals()["metal" + str(i)]
                atom = locals()["atom" + str(i)]
                target_resid = [self.resid[x-1] for x in atom]
        
                print("; I added the pairs - so the zn will not nonbonded interact with the cyx residues")
                for i in range(len(self.atomname)):
                    if self.atomname[i] == 'CA' and self.resid[i] in target_resid:
                        print("%5d%6d%6d" % (metal,i+1,1))
           
    
    def angle_cal(self,  metal1, atom1, metal2, atom2, metal3, atom3, metal4, atom4, metal5, atom5, metal6, atom6, angle_strength):
        for i in [1,2,3,4,5,6]:
            if locals()["metal" + str(i)] !=0:
                metal = locals()["metal" + str(i)]
                atom = locals()["atom" + str(i)]
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
            
x = mr(file, metal1, atom1, metal2, atom2, metal3, atom3, metal4, atom4, metal5, atom5, metal6, atom6, bond_strength, angle_strength)
