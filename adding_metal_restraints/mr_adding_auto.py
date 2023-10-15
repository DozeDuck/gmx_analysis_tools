import argparse
import numpy as np
import re
import math

# Define command line arguments
parser = argparse.ArgumentParser(description='Version 1.0 Metal atoms restraints adding. \n'
                                 'Usage: python mr_adding_multi.py -f rec.gro -m1 888 -a1 358 832 404 -m2 889 -a2 157 524 128',
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-f', '--gro', nargs='+', required=True, default=0, help='input gro file')
parser.add_argument('-ml', '--metallist', default=["MG","MN","ZN","CA"], help='the list of Metals,default value is MG,MN,ZN,CA')
parser.add_argument('-rl', '--reslist', default=["HIS","GLU","ASP","ASN", "CYS", "LYS", "TYR"], help='the list of residues,default value is "HIS","GLU","ASP","ASN", "CYS", "LYS", "TYR"')
parser.add_argument('-al', '--atomlist', default=["ND1","OE1","OE2", "OD1","OD2","ND2","SG", "NZ", "OH"], help='the list of atom name,default value is "ND1","OE1","OE2", "OD1","OD2","ND2","SG", "NZ", "OH"')
parser.add_argument('-d', '--distance', default=0.4, help='the distance default value is 4 angstrom')
parser.add_argument('-n', '--neighbours', default=3, help='the number of neighbours default value is 3 neighbours')

# parser.add_argument('-m1', '--metal1', required=True, default=0, help='the indexs number for your metal atoms, showed in rec.gro file')
# parser.add_argument('-a1', '--atom1', nargs='+', required=True, default=0, help='O/N atom index numbers for each metal atoms, for calculating the b-length, angle')
# parser.add_argument('-m2', '--metal2', default=0, help='the indexs number for your metal atoms, showed in rec.gro file')
# parser.add_argument('-a2', '--atom2', nargs='+',  default=[0], help='O/N atom index numbers for each metal atoms, for calculating the b-length, angle')
# parser.add_argument('-m3', '--metal3', default=0, help='the indexs number for your metal atoms, showed in rec.gro file')
# parser.add_argument('-a3', '--atom3', nargs='+',  default=[0], help='O/N atom index numbers for each metal atoms, for calculating the b-length, angle')
# parser.add_argument('-m4', '--metal4', default=0, help='the indexs number for your metal atoms, showed in rec.gro file')
# parser.add_argument('-a4', '--atom4', nargs='+',  default=[0], help='O/N atom index numbers for each metal atoms, for calculating the b-length, angle')
# parser.add_argument('-m5', '--metal5', default=0, help='the indexs number for your metal atoms, showed in rec.gro file')
# parser.add_argument('-a5', '--atom5', nargs='+', default=[0], help='O/N atom index numbers for each metal atoms, for calculating the b-length, angle')
# parser.add_argument('-m6', '--metal6', default=0, help='the indexs number for your metal atoms, showed in rec.gro file')
# parser.add_argument('-a6', '--atom6', nargs='+', default=[0], help='O/N atom index numbers for each metal atoms, for calculating the b-length, angle')
# parser.add_argument('-t', '--topology_file', default='topol.top', help='the path for topol.top')
parser.add_argument('-l', '--bond', default=200000, help='give the path for your em.mdp file, or this script will use the default_em.mdp')
parser.add_argument('-g', '--angle', default=10000, help='give the path for your nvt.mdp file, or this script will use the default_nvt.mdp')

args = parser.parse_args()

# define the parameters
file = args.gro[0]
metal_list = args.metallist
residue_list = args.reslist
atom_list = args.atomlist
distance_value = args.distance
num_neighbours = int(args.neighbours)
# metal1 = int(args.metal1)
# atom1 = [int(x) for x in args.atom1]
# metal2 = int(args.metal2)
# atom2 = [int(x) for x in args.atom2]
# metal3 = int(args.metal3)
# atom3 = [int(x) for x in args.atom3]
# metal4 = int(args.metal4)
# atom4 = [int(x) for x in args.atom4]
# metal5 = int(args.metal5)
# atom5 = [int(x) for x in args.atom5]
# metal6 = int(args.metal6)
# atom6 = [int(x) for x in args.atom6]
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
        # try:
        #     if self.metals[0]:
        #         self.metal1 = self.metals[0]-1
        #     if self.metals[1]:
        #         self.metal2 = self.metals[1]-1
        #     if self.metals[2]:
        #         self.metal3 = self.metals[2]-1
        #     if self.metals[3]:
        #         self.metal4 = self.metals[3]-1
        #     if self.metals[4]:
        #         self.metal5 = self.metals[4]-1
        #     if self.metals[5]:
        #         self.metal6 = self.metals[5]-1
        # except:
        #     pass
        
        # print out the metals index for double check
        # formatted_metals = ', '.join(map(str,  self.metals))
        metals_name = [self.atomname[i-1] for i in self.metals]
        sentence = "The metals atom index are: {}".format(list(zip(metals_name, self.metals)))
        print(sentence)
        print(self.metal1,self.metal2,self.metal3)
            
        
        
        
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
            
x = mr(file, num_neighbours, distance_value, atom_list, metal_list, residue_list, bond_strength, angle_strength)
