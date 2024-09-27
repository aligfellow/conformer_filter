import rmsd
import glob
import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-n', '--basename', dest='basename', help='specify BASENAME: leave out _XXXX', metavar='BASENAME')
parser.add_argument('-o', '--output_suffix', dest='output', help='specify suffix of OUTPUT folder for data', metavar='OUTPUT')
parser.add_argument('-j', '--job', dest='job', default='GAUSS', help='specify JOB type: GAUSS, ORCA, CREST, XTB (not implemented)', metavar='JOB')
parser.add_argument('-e', '--energy', dest='energy', default='SCF', help='specify ENERGY metric: SCF, GIBBS', metavar='ENERGY') 
parser.add_argument('-r', '--rmsd', dest='rmsd', default='0.25', type=float, help='specify RMSD cutoff', metavar='RMSD')
parser.add_argument('-E', '--Ecutoff', dest='Ecutoff', default='0.0', type=float, help='specify ENERGY cutoff', metavar='Ecutoff')
parser.add_argument('-N', '--Nconf', dest='Nconf', type=int, help='specify number of conformers to keep', metavar='Nconf')
parser.add_argument('-f', '--filter', dest='filter', choices=['rmsd', 'energy'], help='specify filtering method: rmsd or energy', metavar='FILTER')

args = parser.parse_args()

if args.basename is None:
    print('No BASENAME specified')
    exit()
if args.basename[-4:] == '.xyz':
    args.basename = args.basename[:-4]
if args.filter is None and args.Nconf:
    print('No FILTER specified')
    exit()

# check if the files exist
try:
    files = sorted(glob.glob(args.basename + '_????.xyz'))
except:
    print('No files of this basename found')
    exit()

files = sorted(glob.glob(args.basename + '_????.xyz')) # xyz list

if args.job == 'GAUSS':
    try:
        outfiles = sorted(glob.glob(args.basename + '_????.log'))
    except:
        print('No .log files found')
elif args.job == 'ORCA':
    try:
        outfiles = sorted(glob.glob(args.basename + '_????.out'))
    except:
        print('No .out files found')
elif args.job == 'CREST':
    try:
        outfiles = sorted(glob.glob(args.basename + '_????.xyz'))
    except:
        print('No .xyz files found (read energies from .xyz header)')
elif args.job == 'XTB':
    try:
        outfiles = sorted(glob.glob(args.basename + '_????.out'))
    except:
        print('No .out files found')
        

# Initialise matrices
rmsd_matrix = np.zeros((len(files),len(files)))
#erase = np.zeros((n_files,n_files))
energies = np.zeros(len(files))

if args.job == 'GAUSS':
    if args.energy == 'SCF':
        for i in range(len(outfiles)):
            with open(outfiles[i], 'r') as f:
                lines = f.readlines()
            for line in lines:
                if 'SCF Done' in line:
                    energies[i] = float(line.split()[4]) # implicitly will get the last value
                    energies[i] = energies[i] * 627.5 # convert to kcal/mol

elif args.job == 'ORCA':
    if args.energy == 'GIBBS':
        for i in range(len(outfiles)):
            with open(outfiles[i], 'r') as f:
                lines = f.readlines()
            for line in lines:
                if 'Final Gibbs free energy' in line:
                    energies[i] = float(line.split()[5])
                    energies[i] = energies[i] * 627.5 # convert to kcal/mol
    elif args.energy == 'SCF':
        for i in range(len(outfiles)):
            with open(outfiles[i], 'r') as f:
                lines = f.readlines()
            for line in lines:
                if 'FINAL SINGLE POINT ENERGY' in line:
                    energies[i] = float(line.split()[4])
                    energies[i] = energies[i] * 627.5 # convert to kcal/mol

elif args.job == 'CREST':
    if args.energy == 'SCF':
        for i in range(len(files)):
            with open(files[i], 'r') as f:
                lines = f.readlines()
            energies[i] = float(lines[1]) * 627.5 # convert to kcal/mol
            
# Calculate RMSD for each pair of conformers.
# Uses a version of https://github.com/charnley/rmsd with return(rmsd) at the end instead of print(rmsd)
for i in range(len(files)):
    for j in range(i+1,len(files)):
        rmsd_matrix[i,j] = float(rmsd.main([str(files[i]),str(files[j]),'-e']))
        rmsd_matrix[j,i] = rmsd_matrix[i,j]
        
to_delete = [] # list of files to delete
to_keep = [] # list of files to keep

# make the energies relative to the lowest energy
energies = energies - min(energies)

# Delete conformers with energy above threshold
if args.Ecutoff != 0.0:
    for i in range(len(files)):
        if ( energies[i] > args.Ecutoff ):
            to_delete.append(files[i])

# Delete conformers with RMSD below threshold
for i in range(len(files)):
    for j in range(i+1,len(files)):
        if ( rmsd_matrix[i,j] < args.rmsd ): 
            if energies[i] > energies[j]:
                to_delete.append(files[i])
            else:
                to_delete.append(files[j])
to_delete = list(dict.fromkeys(to_delete))
to_keep = [file for file in files if file not in to_delete]

if args.Nconf:
    if args.filter == 'rmsd':
        # incrementally remove conformers with lowest RMSD until Nconf is reached
        current_rmsd = args.rmsd
        while len(to_keep) > args.Nconf:
            current_rmsd += 0.05 # increment RMSD threshold
            
            for i in range(len(files)):
                for j in range(i+1,len(files)):
                    if ( rmsd_matrix[i,j] < current_rmsd ): 
                        if energies[i] > energies[j]:
                            to_delete.append(files[i])
                        else:
                            to_delete.append(files[j])
            to_delete = list(dict.fromkeys(to_delete))
            to_keep = [file for file in files if file not in to_delete]
        print(f"Final RMSD cutoff: {current_rmsd:.2f}")
    elif args.filter == 'energy':
        sorted_indices = np.argsort(energies)
        to_delete = [files[i] for i in sorted_indices[args.Nconf:]]
        to_keep = [files[i] for i in sorted_indices[:args.Nconf]]
        
# Prints how many remaining conformers you have with the chosen thresholds
print('N files remaining: ' + str(len(to_keep)))

# setup folder for the data 
import os
if args.output:
    output_folder= args.output
else: 
    output_folder='RMSD_' + args.job + '_' + args.energy + '_' + str(args.rmsd) + '_Nconf' + str(len(to_keep))
    if args.Ecutoff != 0.0:
        output_folder += '_E' + str(args.Ecutoff)
    if args.Nconf:
        output_folder += '_FILTER' + str(args.filter)  

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
else:
    print('Output folder already exists. Overwriting files.')

with open(f'{output_folder}/delete_files.txt', 'w') as f:
    for file in to_delete:
        f.write(file + '\n')
with open(f'{output_folder}/keep_files.txt', 'w') as f:
    for file in to_keep:
        f.write(file + '\n')
# also make a txt file with the metadata 
with open(f'{output_folder}/metadata.txt', 'w') as f:
    f.write(f'BASENAME: {args.basename}\n')
    f.write(f'JOB: {args.job}\n')
    f.write(f'ENERGY: {args.energy}\n')
    f.write(f'RMSD: {args.rmsd}\n')
    f.write(f'Ecutoff: {args.Ecutoff}\n')
    f.write(f'Nconf_initial {len(files)}\n')
    f.write(f'Nconf: {len(to_keep)}\n')
    f.write(f'FILTERING: {args.filter}\n')

# print the rmsd matrix to a text file
np.savetxt(f'{output_folder}/rmsd_matrix.txt', rmsd_matrix, fmt='%.4f')

# print the energies to a text file
np.savetxt(f'{output_folder}/energies.txt', energies, fmt='%.4f')

