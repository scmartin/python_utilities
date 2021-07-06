#!/usr/bin/env python

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

def lammpsStructure(fname: str, ffieldname: str, periodic: bool):
    """Create a restart for the structure in file <fname>.
    Arguments:
     - fname -- filename for xyz file to read
    Returns: None
    """
    print(f"I'm making a restart with file {fname}")
    L = lammps.lammps(cmdargs=['-screen','none','-log','none','-nocite'])
    L.command('atom_style charge')
    L.command('units real')
    if periodic:
        L.command('boundary p p p')
    else:
        L.command('boundary f f f')
    L.command(f'read_data {fname}')
    L.command('pair_style reax/c NULL')
    L.command('timestep 0.02')
    L.command(f'pair_coeff * * {ffieldname} Cr H O Si')
    L.command('fix 1 all qeq/reax 1 0.0 10.0 1.0e-4 reax/c')
    L.command(f'write_restart {fname}.restart')
    return None


def calcEnergy(fname: str, ffieldname: str) -> float:
    """Evaluate energy of a structure.
    Arguments:
     - fname -- filename for restart file to read
     - ffieldname -- ffield filename
    """
    L = lammps.lammps(cmdargs=['-screen','none','-log','none','-nocite'])
    L.command(f'read_restart {fname}')
    L.command('pair_style reax/c NULL')
    L.command(f'pair_coeff * * {ffieldname} Cr H O Si')
    L.command('fix 1 all qeq/reax 1 0.0 10.0 1.0e-4 reax/c')
    L.command('run 0')
    return L.get_thermo('pe')


def calcForce(fname: str,ffieldname: str):
    """Evaluate atomic forces of a structure.
    Arguments:
     - fname -- filename for restart file to read
     - ffieldname -- ffield filename
    """
    L = lammps.lammps(cmdargs=['-screen','none','-log','none','-nocite'])
    L.command(f'read_restart {fname}')
    L.command('pair_style reax/c NULL')
    L.command(f'pair_coeff * * {ffieldname} Cr Si')
    L.command('fix 1 all qeq/reax 1 0.0 10.0 1.0e-4 reax/c')
    L.numpy.extract_atom('f')

def calcLattice(fname: str, ffieldname: str):
    """Evaluate the relaxed lattice parameters for a structure.
    Arguments:
     - fname -- filename for a restart file to read
     - ffieldname -- filename of a force field file
    """
    L = lammps.lammps(cmdargs=['-log','none','-nocite','-screen','none'])
    L.command(f'read_restart {fname}')
    L.command('replicate 3 3 3')
    L.command('pair_style reax/c NULL')
    L.command(f'pair_coeff * * {ffieldname} Cr H O Si')
    L.command('fix 1 all qeq/reax 1 0.0 10.0 1.0e-4 reax/c')
    L.command('fix 2 all box/relax iso 0.0')
    L.command('minimize 1e-6 1e-6 100000 1000000')
    box = L.extract_box()
    lat = np.array([box[1][0] - box[0][0],
                    box[1][1] - box[0][1],
                    box[1][2] - box[0][2],
                    box[2],box[3],box[4]])
    return(lat/3.0)


def calcPtensor(fname: str, ffieldname: str):
    """Evaluate the pressure tensor for a structure.
    Arguments:
     - fname -- filename for the restart file to read
     - ffieldname -- filename of the force field file
    Returns:
     - p_tensor -- numpy array with the pressure tensor
       components [pxx, pyy, pzz, pxy, pxz, pyz]
    """
    L = lammps.lammps(cmdargs=['-log','none','-nocite','-screen','none'])
    L.command(f'read_restart {fname}')
    L.command('replicate 10 10 10')
    L.command('pair_style reax/c NULL')
    L.command(f'pair_coeff * * {ffieldname} Cr H O Si')
    L.command('fix 1 all qeq/reax 1 0.0 10.0 1.0e-4 reax/c')
    L.command('thermo_style custom pxx pyy pzz pxy pxz pyz')
    L.command('run 0')
    p_tensor = []
    for i in ['pxx','pyy','pzz','pxy','pxz','pyz']:
        p_tensor.append(L.get_thermo(i))
    return np.array(p_tensor)


def read_train_data(train_file):
    """Read in energy training set data.
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
    parser.add_argument('-e', '--energy', nargs='*', help='Energy difference based training set files')
    parser.add_argument('-a', '--lattice', nargs='*', help='Lattice paramter training files')
    parser.add_argument('-F', '--force', nargs='*', help='Training structure for force minimization.' \
                                                        +'Must be in LAMMPS data format')
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
    E_train_files = args.energy
    Lattice_train_files = args.lattice
    F_train_files = args.force
    ff_init = args.ffield
    log = args.log
    nMembers = args.nMembers
    generations = args.generations
    restart = args.restart
    
    bounds = et.parse('bounds.xml')
    with open(ff_init) as f_init:
        xml_ffield = mparams.pfile2xml.file2xml(f_init)
   
    # Read in \delta E training data
    scan_names = []
    scans = {}
    training_scans = np.array([])
    try:
        for train in E_train_files:
            structures, training_vals = read_train_data(train)
            scan_names.append(train)
            scans[train] = (structures)
            training_scans = np.concatenate((training_scans,training_vals[1:]),axis=None)
        # Create LAMMPS restarts for quicker loading
        print("Making energy scans...", flush=True)
        for scan in scan_names:
            for structure in scans[scan]:
                lammpsStructure(structure, ff_init, periodic=False)
    except TypeError as error:
        print("No energy training data given")

    # Read lattice training data
    lattice_names = []
    training_p_tensors = []
    try:
        for train in Lattice_train_files:
            with open(train) as tf:
                for line in tf:
                    lat = line.split()
                    lattice_names.append(lat[0])
                    training_p_tensors.append(lat[1:])
        for structure in lattice_names:
            lammpsStructure(structure, ff_init, periodic=True)
        training_p_tensors = np.array(training_p_tensors).astype(float)
    except TypeError:
        print("No lattice training data given")
    
    # Read force minimization training data
    # try:
    #     for structure in F_train_files:
    #         lammpsStructure(structure,)
    # except TypeError:
    #     print("No force minimization training data given")
    


    def fitnessfunc(params):
        """Calculate the fitness of a given population member."""
        mparams.pfile2xml.xml2file(params,'ffield.reax.temp')
        energies = np.array([])
        # for dfile in g2d.geo2data('header.txt', 'geo.in'):
        if scan_names:
            for scan in scan_names:
                scan_energies = []
                for dfile in scans[scan]:
                    scan_energies.append(calcEnergy(f"{dfile}.restart",'ffield.reax.temp'))
                scan_energies = np.array(scan_energies)
                scan_energies = scan_energies - scan_energies[0]
                energies = np.concatenate((energies, scan_energies[1:]), axis=None)
            energy_fitness = np.sum((training_scans - energies)**2
                                    /training_scans**2)/energies.shape[0]
        else:
            energy_fitness = 0.0
        if lattice_names:
            lattice_stresses = []
            for lat in lattice_names:
                lattice_stresses.append(calcPtensor(f'{lat}.restart','ffield.reax.temp'))
            lattice_stresses = np.array(lattice_stresses)
            lattice_fitness = (np.sum(np.log(1.0+(training_p_tensors 
                                              - lattice_stresses)**2))
                               /lattice_stresses.shape[0])
        else:
            lattice_fitness = 0.0
        # lattice_vols = []
        # for lat in lattice_names:
        #     lattice_vols.append(np.prod(calcLattice(f'{lat}.restart','ffield.reax.temp')[:3]))
        # try:
        #     lattice_fitness = np.sum((np.prod(training_lattices[:,:3],axis=1) -
        #                 np.array(lattice_vols))**2/np.prod(training_lattices[:,:3],axis=1)**2)
        # except IndexError:
        #     lattice_fitness = np.sum((np.prod(training_lattices[:3]) -
        #                 np.array(lattice_vols))**2/np.prod(training_lattices[:3])**2)
        #
        # forces = []
        # for structure in F_train_files:
        #     forces.append(calcForce(f'{structure}','ffield.reax.temp'))
        # forces = np.array(forces)
        # print(forces)
        print(f'Efitness: {energy_fitness}\nLfitness: {lattice_fitness}')
        print(f'Tfitness: {energy_fitness + lattice_fitness} \n')
        return (energy_fitness + lattice_fitness)

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
        try:
            os.mkdir('paramlog')
        except FileExistsError:
            pass
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

    try:
        os.mkdir('energies')
    except FileExistsError:
        pass
    try:
        os.mkdir('lattices')
    except FileExistsError:
        pass
    
    print(f"Population size: {Population.nMembers}")
    print("Starting evolution", flush=True)
    for i in range(generations):
        if i % (generations//10) == 0:
            print(f"generation {i}", flush=True)
            for i in range(Population.nMembers):
                mparams.pfile2xml.xml2file(Population.members[i].params, f'ffield.reax.{i}')
                for scan in scan_names:
                    energies = []
                    for structure in scans[scan]:
                        energies.append(calcEnergy(f"{structure}.restart",f'ffield.reax.{i}'))
                    np.savetxt(f'energies/{scan}_energies_chrom{i}.txt',energies)
                for lat in lattice_names:
                    p_tensor = calcPtensor(f'{lat}.restart',f'ffield.reax.{i}')
                    np.savetxt(f'lattices/{lat}_pressure_chrom{i}.txt',p_tensor)
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
        for scan in scan_names:
            energies = []
            for structure in scans[scan]:
                energies.append(calcEnergy(f"{structure}.restart",f'ffield.reax.{i}'))
            np.savetxt(f'energies/{scan}_energies_chrom{i}.txt',energies)
        for lat in lattice_names:
            p_tensor = calcPtensor(f'{lat}.restart',f'ffield.reax.{i}')
            np.savetxt(f'lattices/{lat}_pressure_chrom{i}.txt',p_tensor)
    print("Program complete!")
