from mergeREAXFFparams import pfile2xml
import xml.etree.ElementTree as et
import copy
import os
import itertools

def merge_atoms(et_list):
    """Make a merged list of atoms from multiple ElementTrees.
    Parameters:
     - et_list -- a list containing force field element trees
    Returns:
     - atoms   -- Element object with atoms as children
    """
    atoms = et.Element('atoms')
    innerWall = True
    asymbols = []
    atomcount = 1
    for ff in et_list:
        for atom in ff.getroot().findall('atoms/atom'):
            if atom.text not in asymbols:
                if float(atom.find('param[@order="32"]').text) == 0.0:
                    print(f"{atom.text} has no inner wall. Parameters will be"
                    + " set accordingly for all atoms")
                    innerWall = False
                newatom = copy.deepcopy(atom)
                newatom.set("order", str(atomcount))
                atoms.append(newatom)
                asymbols.append(atom.text)
                atomcount += 1
            else:
                #ToDo: Add consistency check
                print(f"atom {atom.text} is duplicated")
    atoms.set('n', str(len(atoms)))
    if innerWall == False:
        for atom in atoms:
            for order in [30,31,32]:
                atom.find(f'./param[@order="{order}"]').text = f'{0.0:f}'
    return atoms

def merge_nbody(et_list, section):
    """Merge N-body interactions from multiple ElementTrees.
    Parameters:
     - et_list -- a list containing force field element trees
    Returns:
     - nbody   -- Element object with atoms as children
    """
    nbody = et.Element(section[0])
    previous_connections = []
    os.makedirs('duplicates',exist_ok=True)
    with open(f"duplicates/{section[0]}duplicates.txt", "w") as dup:
        for ff in et_list:
            current_connections = []
            for item in ff.findall(f'{section[0]}/{section[1]}'):
                connection = nbodylist_fromxml(item)
                if (connection not in previous_connections and connection[::-1] not in previous_connections):
                    current_connections.append(connection)
                    nbody.append(item)
                else:
                    #ToDO: write consistency check
                    dup.write(f'{connection} is already here.\n')
            previous_connections += current_connections
    return nbody

def nbodylist_fromxml(nbody_xml):
    """Create ordered list of atoms in nbody interaction"""
    connection = []
    count = 0
    atom = nbody_xml.find(f'./atom[@order="{count}"]')
    while atom is not None:
        connection.append(atom.text)
        count += 1
        atom = nbody_xml.find(f'./atom[@order="{count}"]')
    return connection

def missing_nbody(root, section, possible):
    """Generate xml element for missing n-body interactions"""
    connections = []
    for item in root.findall(f'./{section[0]}/{section[1]}'):
        connection = nbodylist_fromxml(item)
        connections.append(connection)
    missing = []
    for nbody in possible:
        if (nbody not in connections and nbody[::-1] not in connections):
            missing.append(nbody)
    return missing
    

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="This code merges ReaxFF parameter files")
    parser.add_argument('-f', '--filenames', required=True, type=str, nargs='+')
    parser.add_argument('-o', '--output', type=str, default='ffield.merged')
    parser.add_argument('-t', '--title', type=str, default='Merged forcefield file.\n')
    args = parser.parse_args()

    forcefields = []
    for fname in args.filenames:
        with open(fname) as f:
            forcefields.append(pfile2xml.file2xml(f))
    new_ffield = et.ElementTree(et.Element('ffield'))
    missing_ff = et.ElementTree(et.Element('ffield'))
    newroot = new_ffield.getroot()
    missingroot = missing_ff.getroot()
    newroot.set('title', args.title)
    missingroot.set('title', "missing interactions for "+args.title)

    # ToDo: Write comparison function to check general parameters
    newroot.append(forcefields[0].getroot().find('general'))
    
    newroot.append(merge_atoms(forcefields))
    newroot.find('atoms').text = forcefields[0].getroot().find('atoms').text
    
    missingroot.append(forcefields[0].getroot().find('general'))
    
    missingroot.append(newroot.find('atoms'))
    
    sections = [('bonds','bond'),
                ('off_diagonal','pair'),
                ('angles','angle'),
                ('torsions','torsion'),
                ('hbonds','hbond')]
    # ToDo: Generate possible interactions except h-bonds
    possible_interactions = {'bonds':(2,[]), 'off_diagonal':(2,[]), 'angles':(3,[]),
                'torsions':(4,[]), 'hbonds':(3,[])}
    for key, val in possible_interactions.items():
        if key == 'off_diagonal':
            possible_interactions[key] = [list(i) for i in itertools.combinations_with_replacement([atom.text for atom in newroot.find('atoms')],val[0]) if i[0] != i[1]]
        elif key == 'hbonds':
            possible_interactions[key] = [list(i) for i in itertools.combinations_with_replacement([atom.text for atom in newroot.find('atoms')],val[0]) if 'H' in i]
        else:
            possible_interactions[key] = [list(i) for i in itertools.combinations_with_replacement([atom.text for atom in newroot.find('atoms')],val[0])]
    
    # ToDo: Test actual interactions vs possible
    os.makedirs('missing',exist_ok=True)
    for section in sections:
        new_section = merge_nbody(forcefields, section)
        newroot.append(new_section)
        newroot.find(section[0]).text = forcefields[0].getroot().find(section[0]).text
        newroot.find(section[0]).set('n', str(len(newroot.find(section[0]))))
        with open(f'missing/{section[0]}missing.txt', 'w') as f:
            for missing in missing_nbody(newroot, section, possible_interactions[section[0]]):
                f.write(f'{missing}\n')
    # new_ffield.write('mergetest.xml')
    pfile2xml.xml2file(new_ffield, 'ffield.mergetest')
