"""
ToDo
 - implement energy scan where a PyLammps instance is held in memory for each
   structure, and only the forcefield is replaced
 - add parameter logging
 - add multiprocessing or threading
"""
import lammps
import GAOptimizer as GAO
import mergeREAXFFparams as mparams
import numpy as np
import copy
import os
import re


intType = {'off_diagonal':'pair', 'angles':'angle'}

def lammpsStructure(fname: str):
    """Create a restart for the structure in file <fname>.
    Arguments:
     - fname -- filename for xyz file to read
    Returns: None
    """
    print(f"I'm making a PyLammps object with file {fname}")
    L = lammps.lammps(cmdargs=['-screen','none','-log','none','-nocite'])
    L.command('atom_style charge')
    L.command('units real')
    L.command('boundary f f f')
    L.command(f'read_data {fname}')
    L.command('pair_style reax/c NULL')
    L.command('timestep 0.02')
    L.command('pair_coeff * * ffield.reax.Pitman_Shin Cr H O Si')
    L.command('fix 1 all qeq/reax 1 0.0 10.0 1.0e-4 reax/c')
    L.command(f'write_restart {fname}.restart')
    return None


def calcEnergy(fname: str, ffieldname: str) -> float:
    """Evaluate energy of a structure.
    Arguments:
     - fname -- filename for xyz file to read
     - ffieldname -- ffield filename
    """
    L = lammps.lammps(cmdargs=['-screen','none','-log','none','-nocite'])
    L.command(f'read_restart {fname}')
    L.command('pair_style reax/c NULL')
    L.command(f'pair_coeff * * {ffieldname} Cr H O Si')
    L.command('fix 1 all qeq/reax 1 0.0 10.0 1.0e-4 reax/c')
    L.command('run 0')
    return L.get_thermo('pe')


def read_train_data(train_file):
    """Read in training set data.
    Arguments:
     - train_file (str) -- name of the training set file
    
    Returns:
     - structures (list) -- structure names
     - values (numpy array) -- training set value for structure at same index
    """
    
    structures = []
    values = []
    with open(train_file) as f:
        for line in f:
            splitted = line.split()
            structures.append(splitted[0]) 
            values.append(float(splitted[1]))
    return structures, np.array(values)
       

def param_loop(param_tree, func):
    """Loop over parameter element tree and perform some function"""
    root = param_tree.getroot()
    for section in root:
        s_tag = section.tag
        for interaction in section:
            i_tag = interaction.tag
            i_attrib = f'{interaction.get("order")}'
            for group in interaction:
                g_tag = group.tag
                g_attrib = f'{group.get("order")}'
                yield func(s_tag,i_tag,i_attrib,g_tag,g_attrib)


if __name__ == "__main__":
    import geo2data as g2d
    import numpy as np
    import copy
    import xml.etree.ElementTree as et
    import argparse

    
    print("Starting GA program...", flush=True)
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--train', nargs='+', help='Training set files')
    parser.add_argument('-f', '--ffield', help='initial force field file',
                        required=True, type=str)
    parser.add_argument('-l', '--log', help='log parameter values',
                        action='store_true')
    parser.add_argument('-n', '--nMembers', help='Number of members in population',
                        required=False, type=int, default=50)
    parser.add_argument('-g', '--generations', help='Number of generations',
                        required=True, type=int)
    parser.add_argument('-r', '--restart', type=str, 
                        help='prefix for force field files to use in restart')
    
    args = parser.parse_args()
    train_files = args.train
    ff_init = args.ffield
    log = args.log
    nMembers = args.nMembers
    generations = args.generations
    restart = args.restart
    
    bounds = et.parse('bounds.xml')
    with open(ff_init) as f_init:
        xml_ffield = mparams.pfile2xml.file2xml(f_init)
    
    scans = {}
    training_scans = np.array([])
    for train in train_files:
        structures, training_vals = read_train_data(train)
        scans[train] = (structures)
        training_scans = np.concatenate((training_scans,training_vals),axis=None)
    
    # Create LAMMPS restarts for quicker loading
    print("Making scans...", flush=True)
    for scan in scans:
        for structure in scans[scan]:
            lammpsStructure(structure)

    def fitnessfunc(params):
        """Calculate the fitness of a given population member."""
        mparams.pfile2xml.xml2file(params,'ffield.reax.temp')
        energies = np.array([])
        # for dfile in g2d.geo2data('header.txt', 'geo.in'):
        for scan in scans:
            scan_energies = []
            for dfile in scans[scan]:
                scan_energies.append(calcEnergy(f"{dfile}.restart",'ffield.reax.temp'))
            scan_energies = np.array(scan_energies)
            scan_energies = scan_energies - scan_energies[0]
            energies = np.concatenate((energies, scan_energies), axis=None)
        # return np.sum((training_vals[1:] - energies[1:])**2/training_vals[1:]**2)
        return np.sum((training_scans - energies)**2)

    print("Making initial population...", flush=True)
    if not restart:
        Population = GAO.population.population.new_population(nMembers, bounds, xml_ffield,
                                           fitnessfunc)
    else:
        popxml = [mparams.pfile2xml.file2xml(open(i.name)) \
                  for i in os.scandir() if re.match(fr'{restart}\.\d+',i.name)]
        Population = GAO.population.population.read_population(bounds, xml_ffield, popxml, fitnessfunc)
    
    print("initial population made", flush=True)
    def genParamLog(*args):
        if restart:
            mode = 'a'
        else:
            mode = 'w'
        xpath = f'./{args[0]}/'
        xpath += f'{args[1]}[@order="{args[2]}"]/'
        xpath += f'{args[3]}[@order="{args[4]}"]'
        fname = 'paramlog/'
        for i in args:
            fname += i
        fname += ".txt"
        return((xpath,open(fname,mode,1)))

    if log == True:
        files = {}
        for i in param_loop(bounds,genParamLog):
            files[i[0]] = i[1]
        for member in Population.members:
            for chrom in member.chromosome:
                files[chrom].write(member.params.find(chrom).text+'\t')
        for f in files:
            files[f].write('\n')

    print(f"Population size: {Population.nMembers}")
    print("Starting evolution", flush=True)
    for i in range(generations):
        if i % (generations//10) == 0:
            print(f"generation {i}", flush=True)
            for i in range(Population.nMembers):
                mparams.pfile2xml.xml2file(Population.members[i].params, f'ffield.reax.{i}')
                for scan in scans:
                    energies = []
                    for structure in scans[scan]:
                        energies.append(calcEnergy(f"{structure}.restart",f'ffield.reax.{i}'))
                    np.savetxt(f'energies/{scan}_energies_chrom{i}.txt',energies)
        Population.evolve()
        if log == True:
            for member in Population.members:
                for chrom in member.chromosome:
                    files[chrom].write(member.params.find(chrom).text+'\t')
            for f in files:
                files[f].write('\n')

    print("Evolution finished")
    print("Writing final state...", flush=True)
    for i in range(Population.nMembers):
        mparams.pfile2xml.xml2file(Population.members[i].params, f'ffield.reax.{i}')
        for scan in scans:
            energies = []
            for structure in scans[scan]:
                energies.append(calcEnergy(f"{structure}.restart",f'ffield.reax.{i}'))
            np.savetxt(f'energies/{scan}_energies_chrom{i}.txt',energies)
    print("Program complete!")
