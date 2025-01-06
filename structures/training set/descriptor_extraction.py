import re #Import RegEx
import os
import math
import ast
import numpy as np
import pandas as pd
from pathlib import Path
from atomic_numbers import element

# Assign directory names
directory = os.getcwd()
spe_directory = directory + '\\spe'
opt_directory = directory + '\\opt'
steric_map_dir = directory + '\\steric_maps_nocoord_training'
print('Directory = ', directory)


#Input excel table containing metal, atom_a and atom_b name
df = pd.read_excel('input.xls')   

xyz_match = ['X           Y           Z']
LUMO_match = ['Alpha virt. eigenvalues']
HOMO_match = ['Alpha  occ. eigenvalues']
nbo_match = ['Natural Population Analysis']
spe_match = ['SCF Done']


#List for constructing dataframe
df_name = []
df_metal = []
df_atom_a = []
df_atom_b = []
df_angle = []
df_metal_atoma_distance = []
df_metal_atomb_distance = []
df_metal_nbo = []
df_a_nbo = []
df_b_nbo = []
df_HOMO_LUMO = []
df_category = []
df_ligand = []
df_AB_NBO = []
df_MAMB_distance = []
df_spe = []
df_natoms = []

for atomic_number, symbol in element.items():
    if symbol == 'Pd':
        pd_atom_number = int(atomic_number)


#Define functions
        
def getNAtoms(line):
    if 'NAtoms' in line:
        number_list = re.findall('-?\d*\.?\d+',line)            # get NAtoms value
        natoms = int(number_list[0])
        print(natoms)

def load_and_clean_data(filename):
    df = pd.read_csv(filename)
    df = df.drop('Unnamed: 0', axis=1)
    return df

def merge_dataframes(df1, df2, merge_columns):
    merged_df = pd.merge(df1, df2, how="outer", on=merge_columns)
    return merged_df

def save_dataframe(df, filename):
    df.to_csv(filename)




# Aims
# Extract descriptor from output file of spe.
# Extract geometric descriptors (bond angles, bond lengths) and electronic descriptors (nbo charge, HOMO-LUMO gap)
# Input needed: No. for the donor atoms (A and B) as well as metal atom (Pd)

# Bond distance: given cartesian coordinates, find distance between atoms
# Bond angle: given three atoms, find angle


for subdir,dirs,files in os.walk(spe_directory):                  # Loop over each directory, subdirectory and files
    for file in files:                                      # Loop over each file
        if any([file.endswith('.out')]):                    # If file is a .out file
            filename = os.path.join(subdir, file)       # Return path to file
            name = Path(filename).stem           # Extract filename from the end of path and return as a string
            print(name)
         
            for x in df['Filename']:                      # For each element in the "Filename" column of the input file
                if x == name:
                    match_file = df[df['Filename'] == name] 
                    print(match_file)
            
                    metal = match_file.iloc[0,1]                    # Extract metal, atom_a and atom_b numbers from input files
                    atom_a = match_file.iloc[0,2]
                    atom_b = match_file.iloc[0,3]
            
            
#             print(name, 'Metal: ',metal, 'Atom_A: ', atom_a, 'Atom_B: ', atom_b)

            mylines = []

            with open (filename, 'rt') as myfile:       # Open .out for reading text
                # myfile = myfile.read()                # Read the entire file to a string
                for myline in myfile:                    # For each line, stored as myline,
                    mylines.append(myline)               # add its contents to mylines list.


                # Find XYZ Coordinates

                for line in mylines:
                    if 'NAtoms' in line:
                        number_list = re.findall('-?\d*\.?\d+',line)            # get NAtoms value
                        natoms = int(number_list[0])
                        print(natoms)

                for line in mylines:
                    for phrase in xyz_match:                                # iterate through each phrases
                        if phrase in line:                                          # check if phrase is in line

                                                                                    # For loop for generating the XYZ coordinates
                            count = 0
                            xyz = []
                            line_number = mylines.index(line) + 2
                            while count < natoms:
                                count = count + 1
                                xyz.append(mylines[line_number])
                                line_number = line_number + 1
                            
#                             # Convert into xyz file 
                            

#                             with open('xyz.txt', 'w') as filehandle:                   # Save xyz coordinates (list) into a txt file
#                                 for listitem in xyz:
#                                     filehandle.write(listitem)

                x_coord = []
                y_coord = []
                z_coord = []
                atom_symbol =[]                        
                
                # Generate XYZ file in .txt form, then find xyz coordinates for metal, atom_a and atom_b
                for line in xyz:
                    number_list = re.findall('-?\d*\.?\d+',line)
                    atom_number = int(number_list[0])
                    element_number = int(number_list[1])
                    # print('ELement number:', element_number)
                    atom_x = float(number_list[3])
                    atom_y = float(number_list[4])
                    atom_z = float(number_list[5])
                    
                    # Make xyz coord into .txt file for percent buried volume
                    x_coord.append(atom_x)
                    y_coord.append(atom_y)
                    z_coord.append(atom_z)
                    atom_symbol.append(element[str(element_number)])
                    name_xyz = name + '_xyz.txt'            
                    xyz_df = pd.DataFrame([atom_symbol,x_coord,y_coord,z_coord])
                    xyz_df = xyz_df.transpose()
#                     xyz_df.to_csv(name_xyz, header=False, index=False, sep = " ")

                    # Identify metal, atom_a and atom_b numbers and obtain their respective x,y,z coords
                    if atom_number == metal:
                        metal_x = atom_x
                        metal_y = atom_y
                        metal_z = atom_z
                        if element_number != pd_atom_number:
                            print('Pd not assigned',filename) # Checking if metal atom was assigned correctly
                    elif atom_number == atom_a:
                        a_x = atom_x
                        a_y = atom_y
                        a_z = atom_z
                    elif atom_number == atom_b:
                        b_x = atom_x
                        b_y = atom_y
                        b_z = atom_z
                        
                
                    

                #Extract metal atom distance
                
                metal_atoma_distance = math.sqrt((metal_x - a_x) ** 2 + (metal_y - a_y) ** 2 + (metal_z - a_z) ** 2)    
                metal_atomb_distance = math.sqrt((metal_x - b_x) ** 2 + (metal_y - b_y) ** 2 + (metal_z - b_z) ** 2)

                metal_atoma_vector = [a_x - metal_x, a_y - metal_y, a_z - metal_z]
                metal_atomb_vector = [b_x - metal_x, b_y - metal_y, b_z - metal_z]

                dot_pdt = np.dot(metal_atoma_vector,metal_atomb_vector)
                    
                #Extract A-M-B angle by using arcos rule
                angle = int(math.degrees(math.acos (dot_pdt / ( metal_atoma_distance * metal_atomb_distance ))))         

#                 print('A-M-B Angle: ', angle)
#                 print('Metal-AtomA Distance: ', metal_atoma_distance)
#                 print('Metal-AtomB Distance: ', metal_atomb_distance)



                # Find HOMO - LUMO energy gap

                LUMO = []
                for phrase in LUMO_match:
                    for line in mylines:
                        if phrase in line:
                            LUMO.append(line)
                number_list = re.findall('-?\d*\.?\d+',LUMO[0])            # get value from line containing LUMO
                LUMO_energy = float(number_list[0])                         # get first value from the line, which represents the LUMO

                HOMO = []
                for phrase in HOMO_match:
                    for line in mylines:
                        if phrase in line:
                            HOMO.append(line)
                number_list = re.findall('-?\d*\.?\d+',HOMO[-1])            # get value from line containing LUMO
                HOMO_energy = float(number_list[-1])                        # get lsat value from line, which represents HOMO

                HOMO_LUMO = LUMO_energy - HOMO_energy                       # HOMO LUMO gap, in A.U.
#                 print('HOMO-LUMO Gap: ', HOMO_LUMO)

                # Find NBO charge

                for phrase in nbo_match:                                       # iterate through each phrases
                    for line in mylines:
                        if phrase in line:                                          # check if phrase is in line

                            count = 0
                            nbo_xyz = []
                            line_number = mylines.index(line) + 6
                            while count < natoms:
                                count = count + 1
                                nbo_xyz.append(mylines[line_number])
                                line_number = line_number + 1

                            with open('nbo_xyz.txt', 'w') as filehandle:                   # Save nbo_xyz  into a txt file
                                for listitem in nbo_xyz:
                                    filehandle.write(listitem)

                                                                                            # Extract NBO charge for metal, atom A and atom B
                for line in nbo_xyz:
                    number_list = re.findall('-?\d*\.?\d+',line)
                    atom_number = int(number_list[0])

                    if atom_number == metal:
                        metal_nbo = float(number_list[1])

                    elif atom_number == atom_a:
                        a_nbo = float(number_list[1])

                    elif atom_number == atom_b:
                        b_nbo = float(number_list[1])
                        
                # Find SPE 
                
                for phrase in spe_match:
                    for line in mylines:
                        if phrase in line:
                            spe = line.replace("SCF Done:  E(RM06) =  " , "")
                            l = len(spe)
                            spe = float(spe[:l-27])
                            df_spe.append(spe)
                            
                # Naming category
                
                if 'nocoord' in name:
                    category = 'nocoord'
                    df_category.append('nocoord')
                    name = name.replace('-spe-nocoord',"")
                elif 'cocoord' in name:
                    category = 'cocoord'
                    df_category.append('cocoord')
                    name = name.replace('-spe-cocoord',"")
                elif 'etcoord' in name:
                    category = 'etcoord'
                    df_category.append('etcoord') 
                    name = name.replace('-spe-etcoord',"")
                else:
                    print('error, no category assigned')
                 
                print(name, category, natoms)
                df_name.append(name)
                df_metal.append(metal)
                df_atom_a.append(atom_a)
                df_atom_b.append(atom_b)
                df_angle.append(angle)
                df_metal_atoma_distance.append(metal_atoma_distance)
                df_metal_atomb_distance.append(metal_atomb_distance)
                df_metal_nbo.append(metal_nbo)
                df_a_nbo.append(a_nbo)
                df_b_nbo.append(b_nbo)
                df_HOMO_LUMO.append(HOMO_LUMO)
                df_AB_NBO.append(a_nbo-b_nbo)
                df_MAMB_distance.append(metal_atoma_distance-metal_atomb_distance)
                df_natoms.append(natoms)
#                 df_ligand.append(name[:-12])
                
                
data = {'Name': df_name,
        'Category' : df_category,
       'Metal': df_metal, 
       'Atom_A' : df_atom_a, 
       'Atom_B' : df_atom_b, 
       'A-M-B Angle' : df_angle, 
       'M-A Distance' : df_metal_atoma_distance, 
       'M-B Distance' : df_metal_atomb_distance, 
       'Metal_NBO' : df_metal_nbo, 
       'A_NBO' : df_a_nbo, 
       'B_NBO' : df_b_nbo, 
       'HOMO_LUMO' : df_HOMO_LUMO,
       'A-B_NBO' : df_AB_NBO,
       'MA_MB_Difference' : df_MAMB_distance,
       'SPE': df_spe}       
                
    




descriptor_data = pd.DataFrame(data)
descriptor_data = descriptor_data.sort_values('Category')
print(descriptor_data)
descriptor_data.to_csv("output_nbo.csv")      # Generate csv file with data table of nbo descriptors
checking_data = {'Name': df_name,
                 'Category' : df_category,
                 'NAtoms': df_natoms} 

checking_df = pd.DataFrame(checking_data)
checking_df.to_csv("checking_df.csv") 

# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Aims
# Extract gibbs free energy from output opt files.

df_filename = []
df_gibbs = []
df_category = []
df_ligand = []

match_list = ['Thermal correction to Gibbs Free Energy']

rootdir = os.getcwd()                                       # Get current working directory

count = 0

for subdir,dirs,files in os.walk(opt_directory):                  # Loop over each directory, subdirectory and files
    for file in files:                                      # Loop over each file
        if any([file.endswith('.out')]):                    # If file is a .out file
            filename = os.path.join(subdir, file)       # Return path to file

            mylines = []

            with open (filename, 'rt') as myfile:       # Open .out for reading text
                # contents = myfile.read()                # Read the entire file to a string
                for myline in myfile:                    # For each line, stored as myline,
                    mylines.append(myline)               # add its contents to mylines.
                    
                for line in mylines:
                    if 'NAtoms' in line:
                        number_list = re.findall('-?\d*\.?\d+',line)            # get NAtoms value
                        natoms = int(number_list[0])


                for phrase in match_list:                                       # iterate through each phrases
                    for line in mylines:
                        if phrase in line:                                          # check if phrase is in line
                            name = Path(filename).stem
                            line = line.replace("Thermal correction to Gibbs Free Energy=         " , "")
                            l = len(line)
#                             line = line[:l-27]
                            gfe_corr = float(line)
                            df_gibbs.append(gfe_corr)
                            
                            if 'nocoord' in name:
                                category = 'nocoord'
                                df_category.append('nocoord')
                                name = name.replace('-opt-nocoord',"")
                            elif 'cocoord' in name:
                                category = 'cocoord'
                                df_category.append('cocoord')
                                name = name.replace('-opt-cocoord',"")
                            elif 'etcoord' in name:
                                category = 'etcoord'
                                df_category.append('etcoord') 
                                name = name.replace('-opt-etcoord',"")
                            else:
                                print('error, no category assigned')
                            df_filename.append(name)
                            print(name, category, natoms)

                            if name in checking_df['Name'].tolist():
                                to_check_df = checking_df[(checking_df['Name'] == name) & (checking_df['Category'] == category)]

                                if len(to_check_df) == 0:
                                    print('opt name not in spe list')
                                
                                elif natoms != to_check_df['NAtoms'].tolist()[0]:
                                    print('opt and spe file doesnt match')
                                    
                            else:
                                print('opt name not in spe list')
 

                            
delta_g_data = {'Name': df_filename, 'Category': df_category, 'Gibbs Free Energy Correction' : df_gibbs}      

descriptor_data = pd.DataFrame(delta_g_data)
descriptor_data.to_csv("output_deltaG.csv")                              
                            
                            # number_list = re.findall('-?\d*\.?\d+',line)            # get value
                            # for x in number_list:                                   # print value
                            #     print(phrase,': ', float(x))

# Merge nbo and gibbs free energy dataset
# Calculate delta G difference between etcoord and nocoord

# Clean up dataframes
df_opt = load_and_clean_data('output_deltaG.csv')
df_nbo = load_and_clean_data('output_nbo.csv')


outer_merge = merge_dataframes(df_nbo, df_opt, ["Name", "Category"])    #Merge nbo and opt dataframes
save_dataframe(outer_merge, "output_merge.csv")

# Seperate dataframe into two dataframes:  df_nocoord and df_etcoord
category_unique = outer_merge["Category"].unique()
for x in list(category_unique):
    df = outer_merge[outer_merge.Category == x]
    df.to_csv("df_" + x + ".csv") 


df_nocoord = load_and_clean_data('df_nocoord.csv')
df_etcoord = load_and_clean_data('df_etcoord.csv')
df_nocoord.rename(columns = {'SPE':'SPE_nocoord', 
                             'Gibbs Free Energy Correction':'Gibbs_Free_Energy_Correction_nocoord'}, inplace = True)
df_etcoord.rename(columns = {'SPE':'SPE_etcoord', 
                             'Gibbs Free Energy Correction':'Gibbs_Free_Energy_Correction_etcoord'}, inplace = True)
df_etcoord = df_etcoord[['Name','SPE_etcoord','Gibbs_Free_Energy_Correction_etcoord']]

#Re-merge the two dataframes horizontally
new_merge = pd.merge(df_nocoord,df_etcoord, how="outer", on=["Name"])

# Add new column with the difference in Gibbs free energy between the etcoord and nocoord intermediates
new_merge['GFE_etcoord-nocoord'] = (627.5*(new_merge['Gibbs_Free_Energy_Correction_etcoord'] + new_merge['SPE_etcoord'] -
                                   (new_merge['Gibbs_Free_Energy_Correction_nocoord'] + new_merge['SPE_nocoord']) + 78.51996755))

# Remove columns which are not descriptors
new_merge = new_merge.drop(['Category',
                            'Metal',
                            'Atom_A',
                            'Atom_B',
                            'SPE_nocoord',
                            'SPE_etcoord',
                            'Gibbs_Free_Energy_Correction_nocoord',
                            'Gibbs_Free_Energy_Correction_etcoord'],axis=1)

# Merge with percent buried volume
df_buriedvol = pd.read_csv('percent_buried_volume.csv')
new_merge_2 = pd.merge(new_merge,df_buriedvol, how="outer", on=["Name"])

print(new_merge_2)

# Adding new steric descriptors from steric map raw data (%Filled_x, %Filled_y and weighted %filled)

# Import list of assigned nonalternating examples for measuring coefficient later
nonalt_list = pd.read_csv('nonalt_list.csv',header=None)
nonalt_list = nonalt_list.values.tolist()
nonalt_list = [val for sublist in nonalt_list for val in sublist]
    

#Set up new dataframe for steric descriptors

data_dict = {'Name':[],
        'Axis':[],
        'Boundary':[],
        'Threshold':[],
        'Unfilled':[],
        'Filled':[],
        'Percent_Filled':[],
        'Nonalt':[]
       }

new_df = pd.DataFrame(data_dict)
# new_df = pd.DataFrame(columns=['Name', 'Axis', 'Boundary', 'Threshold', 'Unfilled', 'Filled', 'Percent_Filled', 'Nonalt'])

# Parse through each steric map csv file to obtain %Filled for each axis, boundary and threshold        
for filename in os.listdir(steric_map_dir):
    csv_name = filename.replace(".csv","")
    name = csv_name.replace("-nocoord-map","")   #obtain catalyst name as a variable
    filepath = os.path.join(steric_map_dir,filename)
    print(name,filepath)
    dataset = pd.read_csv(filepath,header=None)
    dataset.columns = ["X-coord", "Y-coord", "Z"]

    # Set up dataset lists
    axis_list = ["X-coord",'Y-coord']
    min_list = [-1.00, -2.00, -3.00]
    max_list = [1.00, 2.00, 3.00]

    # Define upper and lower boundaries within the dataset
    for axis in axis_list:
        for min_boundary, max_boundary in zip(min_list,max_list):
            df_min = dataset[dataset[axis] >= min_boundary]
            df_minmax = df_min[df_min[axis] <= max_boundary]

    # Identify %Filled for each Z threshold based on voxel value

            threshold_list = list(np.arange(0,2,0.25))  # Define Z thresholds with 0.25 intervals
            boundary = str(min_boundary)+'-'+str(max_boundary)

            for threshold in threshold_list:
                df_minmax.loc[df_minmax['Z'] <= threshold, 'Filled_Unfilled'] = 'Unfilled'
                df_minmax.loc[df_minmax['Z'] > threshold, 'Filled_Unfilled'] = 'Filled'

                if name in nonalt_list:
                    nonalt = 1
                else:
                    nonalt = 0 

                if df_minmax.Filled_Unfilled.value_counts().Unfilled == len(df_minmax.index):  # If all voxel is unfilled
                    filled = 0
                    unfilled = df_minmax.Filled_Unfilled.value_counts().Unfilled
                    filled_percent = 0

                else:
                    unfilled, filled = df_minmax['Filled_Unfilled'].value_counts()  # If some voxel is filled
                    filled_percent = 100 * filled/(unfilled+filled)

                new_row = {'Name':name,
                           'Axis':axis,
                           'Boundary':boundary,
                           'Threshold':threshold, 
                           'Unfilled':unfilled, 
                           'Filled':filled, 
                           'Percent_Filled':filled_percent,
                           'Nonalt':nonalt}
                print('New Row:', new_row)
                # new_df = pd.concat([new_df, new_row], ignore_index=True)
                new_df = new_df._append(new_row, ignore_index=True)

    #Add inversed dataset
    for axis in axis_list:
        for min_boundary, max_boundary in zip(min_list,max_list):
            df_overmax = dataset[dataset[axis] >= max_boundary]                                #Inversed boundaries
            df_lowermin = dataset[dataset[axis] <= min_boundary]
            df_combined = pd.concat([df_overmax,df_lowermin])


           # Identify %Filled for each Z threshold

            threshold_list = list(np.arange(0,2,0.25))
            boundary = '-3.45 - '+str(min_boundary)+' and '+str(max_boundary) + ' - 3.55'

            for threshold in threshold_list:
                df_combined.loc[df_combined['Z'] <= threshold, 'Filled_Unfilled'] = 'Unfilled'
                df_combined.loc[df_combined['Z'] > threshold, 'Filled_Unfilled'] = 'Filled'


                if name in nonalt_list:
                    nonalt = 1
                else:
                    nonalt = 0 


                if df_combined.Filled_Unfilled.value_counts().Unfilled == len(df_combined.index):
                    filled = 0
                    unfilled = df_combined.Filled_Unfilled.value_counts().Unfilled
                    filled_percent = 0

                else:
                    unfilled, filled = df_combined['Filled_Unfilled'].value_counts()
                    filled_percent = 100 * filled/(unfilled+filled)

                new_row = {'Name':name,
                           'Axis':axis,
                           'Boundary':boundary,
                           'Threshold':threshold, 
                           'Unfilled':unfilled, 
                           'Filled':filled, 
                           'Percent_Filled':filled_percent,
                           'Nonalt':nonalt}
                new_df = new_df._append(new_row, ignore_index=True)

              
    
new_df = new_df.sort_values(by =['Axis', 'Boundary', 'Threshold'])
                

new_df.to_csv('steric_descriptors-training.csv')



# Obtaining coefficient for each axis/boundary/threshold combination

group_df = new_df.groupby(['Axis', 'Boundary','Threshold'])
# print(new_df.groupby(['Axis', 'Boundary','Threshold']).size())
# new_df.to_csv('group_df_steric_descriptors-training.csv')

coefficient_list = []
group_name_list = []
axis_list = []
boundary_list = []
threshold_list = []

x_coefficient = 0
y_coefficient = 0

for group_name, group_data in group_df:
#     print(group_df.get_group(key)) 
    axis_group_list = group_data["Axis"].values.tolist()
    axis_array = np.array(axis_group_list)
    axis_value = pd.unique(axis_array)
    axis_list.append(' '.join(axis_value))
    
    boundary_group_list = group_data["Boundary"].values.tolist()
#     print(boundary_group_list)
    boundary_array = np.array(boundary_group_list)
    boundary_value = pd.unique(boundary_array)
    boundary_list.append(' '.join(boundary_value))
    
    threshold_group_list = group_data["Threshold"].values.tolist()
#     print(threshold_group_list)
    threshold_array = np.array(threshold_group_list)
    threshold_value = pd.unique(threshold_array)
    threshold_value = str(threshold_value)[1:-1]
    threshold_list.append(threshold_value)

    if group_name[0] == 'X-coord':
        column_1 = group_data["Percent_Filled"]
        column_2 = group_data["Nonalt"]        
        coefficient = (column_1.corr(column_2))**2
        print('Name: ', group_name , 'Coefficient = %.2f' % coefficient)
        coefficient_list.append(coefficient)
        group_name_list.append(group_name)
        
        # Iterate through all combinations of boundaries/z-theshold until highest coefficient is obtained for the x-axis
        if coefficient > x_coefficient:
            x_coefficient = coefficient
            x_group_name = group_name
        
        
        
    elif group_name[0] == 'Y-coord':
        column_1 = group_data["Percent_Filled"]
        column_2 = group_data["Nonalt"]        
        coefficient = (column_1.corr(column_2))**2
        print('Name: ', group_name ,'Coefficient = %.2f' % coefficient)
        coefficient_list.append(coefficient)
        group_name_list.append(group_name)
        
        # Iterate through all combinations of boundaries/z-theshold until highest coefficient is obtained for the y-axis
        if coefficient > y_coefficient:
            y_coefficient = coefficient
            y_group_name = group_name

# Obtain combination of Axis, Boundary with highest coefficient for X & Y coordinates
print(x_group_name, 'Coefficient = %.2f' % x_coefficient)
print(y_group_name, 'Coefficient = %.2f' % y_coefficient)

#Store group combination into a txt file
with open('x_max_steric.txt','w') as file:
    file.write(str(x_group_name))
with open('y_max_steric.txt','w') as file:
    file.write(str(y_group_name))
    

#Obtain dataframe group with highest coefficient
x_max_steric = group_df.get_group(tuple(x_group_name))
print("x-max-steric: ", x_max_steric)
y_max_steric = group_df.get_group(tuple(y_group_name))
print("y-max-steric: ", y_max_steric)


#Obtain dataframe with group name & coefficient
coefficient_tuples = list(zip(axis_list, boundary_list, threshold_list, coefficient_list))
coefficient_df = pd.DataFrame(coefficient_tuples, columns=['Axis', 'Boundary','Threshold','Coefficient'])
coefficient_df.to_csv('coefficient_df.csv')

 
# printing correlation list
print("The original correlation list : " + str(coefficient_list))
 
# declaring magnitude of repetition with the number of training samples
K = len(new_merge_2)     
 
new_coefficient_list=[]
# repeat elements K times
for i in coefficient_list:
    new_coefficient_list.extend([i]*K)
    
# print("The list after adding elements : " + str(len(new_coefficient_list)))
new_df["Coefficient"] = new_coefficient_list

# Multiply each %Filled value with the correlation coefficient of the respective axis/boundary/z-threshold combination to obtain
# weighted %filled. Add into a new column in the dataframe
new_df['Weighted_percent_filled'] = new_df['Percent_Filled'] * new_df['Coefficient']
print(new_df)
new_df.to_csv('steric_descriptors-training-with-coefficient.csv')


# Organize table according to the name of the training sample. 
# Obtain average of all weighted %filled values across various combinations of axis/boundary/z-threshold

group_df = new_df.groupby(['Name'])
name_list = []
average_weighted_percent_filled = []
for group_name, group_data in group_df:
    group_name_string = group_name[0]   #Convert group name from tuple to string
    group_name = group_name_string
    weighted_values = group_data["Weighted_percent_filled"]  
    name_list.append(group_name)
    average = (sum(weighted_values)/len(weighted_values))
    average_weighted_percent_filled.append(average)
#     print (group_name, average)

avg_weighted_percent_filled_tuples = list(zip(name_list, average_weighted_percent_filled))
avg_weighted_percent_filled_df = pd.DataFrame(avg_weighted_percent_filled_tuples, columns=['Name','Weighted %Filled'])

print(avg_weighted_percent_filled_df)

x_max_steric = x_max_steric.drop(['Axis','Boundary','Threshold','Unfilled','Filled','Nonalt'],axis=1)
x_max_steric = x_max_steric.rename(columns={"%Filled": "%Filled_x"})
y_max_steric = y_max_steric.drop(['Axis','Boundary','Threshold','Unfilled','Filled'],axis=1)
y_max_steric = y_max_steric.rename(columns={"%Filled": "%Filled_y"})


new_merge_3 = pd.merge(pd.merge(pd.merge(new_merge_2,x_max_steric, how="outer", on=["Name"])
                                ,y_max_steric, how="outer", on=["Name"])
                                ,avg_weighted_percent_filled_df, how="outer", on=["Name"])

new_merge_3.to_csv("all_descriptor_data.csv")  






