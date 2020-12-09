import xml.etree.ElementTree as et
from typing import TextIO

# ToDo next: make atoms elements for multiparticle interactions

def genemptyXML():
    "Generate the empty XML representation of reaxFF forcefield file"

    emptyTree = et.ElementTree(et.Element('ffield'))
    root = emptyTree.getroot()
    sections = ['general','atoms','bonds','off_diagonal','angles','torsions','hbonds']
    for i in sections:
        et.SubElement(root,i)
    return emptyTree


def file2xml(f: TextIO) -> et.ElementTree: 
    """Create xml document from ReaxFF parameter file object.

    Arguments:
    f -- file object refering to the parameter file

    Returns:
    xparams -- element tree of the parameter file
    """

    xparams = genemptyXML()
    root = xparams.getroot()
    title = f.readline()
    root.set('title',title)
    
    #set the general parameters
    general = f.readline().split('!')
    elem_general = root.find("./general")
    elem_general.set('n',general[0].strip())
    elem_general.text = general[1].rstrip()

    line = f.readline()
    count = 0
    # iterate through general parameter lines
    while "Nr of atoms" not in line:
        count += 1
        splitLine = line.split('!')
        sub = et.SubElement(elem_general, "param")
        sub.set("order", str(count))
        sub.text = splitLine[0].strip()
        sub.set("description", splitLine[1].rstrip())
        line = f.readline()
    if count != int(elem_general.get('n')):
        print("Number of general parameters doesn't match")

    # generate the atoms section of the xml doc
    splitLine = line.split('!')
    atoms = root.find("./atoms")
    atoms.set("n", splitLine[0].strip())
    text = splitLine[1]
    for i in range(3):
        text += f.readline().lstrip()
    atoms.text = text
    line = f.readline()

    # iterate through atoms
    atomcount = 1
    while "Nr of bonds" not in line:
        count = 0
        splitLine = line.split()
        atom = et.SubElement(atoms, "atom")
        atom.set("order", str(atomcount))
        atom.text = splitLine[0]
        for item in splitLine[1:]:
            count += 1
            p = et.SubElement(atom, "param")
            p.set("order", str(count))
            p.text = item
        for i in range(3):
            splitLine = f.readline().split()
            for item in splitLine:
                count += 1
                p = et.SubElement(atom, "param") 
                p.set("order", str(count))
                p.text = item
        line = f.readline()
        atomcount += 1
    
    # generate the bonds section of the xml doc
    splitLine = line.split('!')
    bonds = root.find("./bonds")
    bonds.set("n", splitLine[0].strip())
    text = splitLine[1]
    text += f.readline().lstrip()
    bonds.text = text
    line = f.readline()

    # iterate through bonds
    bondcount = 1
    while "Nr of off-diagonal" not in line:
        count = 0
        splitLine = line.split()
        elems = [atoms.find(f"./atom[@order='{etype}']").text for etype in splitLine[:2]]
        bond = et.SubElement(bonds, "bond")
        bond.set("order",str(bondcount))
        for ind,elem in enumerate(elems):
            atom = et.SubElement(bond, "atom")
            atom.set('order', str(ind))
            atom.text = elem
        for item in splitLine[2:]:
            count += 1
            p = et.SubElement(bond, "param")
            p.set("order", str(count))
            p.text = item
        splitLine = f.readline().split()
        for item in splitLine:
            count += 1
            p = et.SubElement(bond, "param") 
            p.set("order", str(count))
            p.text = item
        bondcount += 1
        line = f.readline()
    
    # generate the off-diagonal section of the xml doc
    splitLine = line.split('!')
    od = root.find("./off_diagonal")
    od.set("n", splitLine[0].strip())
    text = splitLine[1]
    od.text = text
    line = f.readline()

    # iterate through pairs
    odcount = 1
    while "Nr of angles" not in line:
        count = 0
        splitLine = line.split()
        elems = [atoms.find(f"./atom[@order='{etype}']").text for etype in splitLine[:2]]
        pair = et.SubElement(od, "pair")
        pair.set("order", str(odcount))
        for ind,elem in enumerate(elems):
            atom = et.SubElement(pair, "atom")
            atom.set('order', str(ind))
            atom.text = elem
        for item in splitLine[2:]:
            count += 1
            p = et.SubElement(pair, "param")
            p.set("order", str(count))
            p.text = item
        odcount += 1
        line = f.readline()
    
    # generate the angles section of the xml doc
    splitLine = line.split('!')
    angles = root.find("./angles")
    angles.set("n", splitLine[0].strip())
    text = splitLine[1]
    angles.text = text
    line = f.readline()

    # iterate through angles
    anglecount = 1
    while "Nr of torsions" not in line:
        count = 0
        splitLine = line.split()
        elems = [atoms.find(f"./atom[@order='{etype}']").text for etype in splitLine[:3]]
        angle = et.SubElement(angles, "angle")
        angle.set("order", str(anglecount))
        for ind,elem in enumerate(elems):
            atom = et.SubElement(angle, "atom")
            atom.set('order', str(ind))
            atom.text = elem
        for item in splitLine[3:]:
            count += 1
            p = et.SubElement(angle, "param")
            p.set("order", str(count))
            p.text = item
        anglecount += 1
        line = f.readline()
    
    # generate the torsions section of the xml doc
    splitLine = line.split('!')
    torsions = root.find("./torsions")
    torsions.set("n", splitLine[0].strip())
    text = splitLine[1]
    torsions.text = text
    line = f.readline()

    # iterate through torsions
    torsioncount = 1
    while "Nr of hydrogen" not in line:
        count = 0
        splitLine = line.split()
        elems = [atoms.find(f"./atom[@order='{etype}']").text if etype != '0' else '0' for etype in splitLine[:4]]
        tors = et.SubElement(torsions, "torsion")
        tors.set("order", str(torsioncount))
        for ind,elem in enumerate(elems):
            atom = et.SubElement(tors, "atom")
            atom.set('order', str(ind))
            atom.text = elem
        for item in splitLine[4:]:
            count += 1
            p = et.SubElement(tors, "param")
            p.set("order", str(count))
            p.text = item
        torsioncount += 1
        line = f.readline()

    # generate h-bond section of xml doc
    splitLine = line.split('!')
    hbonds = root.find("./hbonds")
    hbonds.set("n", splitLine[0].strip())
    text = splitLine[1]
    hbonds.text = text
    line = f.readline()

    # iterate through h-bonds
    hbondcount = 1
    while line:
        count = 0
        splitLine = line.split()
        elems = [atoms.find(f"./atom[@order='{etype}']").text for etype in splitLine[:3]]
        hbond = et.SubElement(hbonds, "hbond")
        hbond.set("order", str(hbondcount))
        for ind,elem in enumerate(elems):
            atom = et.SubElement(hbond, "atom")
            atom.set('order', str(ind))
            atom.text = elem
        for item in splitLine[3:]:
            count += 1
            p = et.SubElement(hbond, "param")
            p.set("order", str(count))
            p.text = item
        hbondcount += 1
        line = f.readline()


    return xparams


def xml2file(tree, filename):
    """generate a ReaxFF forcefield file from an xml representation
    """
    f = open(filename, "w")
    root = tree.getroot()
    f.write(root.get('title'))
    general = root.find('general')
    f.write(f"{int(general.get('n')):3d}{' '*7}!{general.text}\n")
    for param in general:
        f.write(f"{float(param.text):10.4f} !{param.get('description')}\n")
    sections = ['atoms', 'bonds', 'off_diagonal', 'angles', 'torsions', 'hbonds']
    formats = {'atoms': [(1,"{:>3s}"),(8,"{:9.4f}"),4],
               'bonds': [(2,"{:>3s}"),(8,"{:9.4f}"),2],
               'off_diagonal': [(2,"{:>3s}"),(6,"{:9.4f}"),1],
               'angles': [(3,"{:>3s}"),(7,"{:9.4f}"),1],
               'torsions': [(4,"{:>3s}"),(7,"{:9.4f}"),1],
               'hbonds': [(3,"{:>3s}"),(4,"{:9.4f}"),1]}
    atoms = {}
    for section in sections:
        f.write("{:5s}".format(root.find(f"{section}").get('n')))
        f.write("!"+root.find(f"{section}").text.rstrip()+"\n")
        formatter = formats[section]
        for child in root.find(section):
            if section == 'atoms':
                atoms[child.text] = child.get("order")
            else:
                atom_order = {}
                for atom in child.findall("atom"):
                    atom_order[atom.get('order')] = atom.text
            line = ''
            linecount = 0
            paramcount = 1
            while linecount < formatter[2]:
                for i in range(formatter[0][0]):
                    if linecount == 0 and section == 'atoms':
                        line += formatter[0][1].format(child.text)
                    elif linecount == 0 and section != 'atoms':
                        try:
                            line += formatter[0][1].format(atoms[child.find(f"atom[@order='{i}']").text])
                        except KeyError:
                            line += formatter[0][1].format("0")
                    else:
                        line += formatter[0][1].format(" ")
                for j in range(formatter[1][0]):
                    line += formatter[1][1].format(float(child.find(f"param[@order='{paramcount}']").text))
                    paramcount += 1
                linecount += 1
                f.write(line+'\n')
                line = ""
    f.close()


if __name__ == "__main__":
    pass
