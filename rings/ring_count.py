#!/usr/bin/env python

def find_next(nodes, count, cycles, sisi):
    if count <= cutoff:
        for k in sisi:
            if set(k) == set((nodes[0], nodes[-1])) \
            and set(nodes) not in [set(i) for i in cycles[count]]:
                cycles[count].append(nodes)
            elif set(k) & set(nodes) == set([nodes[-1]]):
                (new_node,) = set(k) - set([nodes[-1]])
                find_next(nodes + [new_node], count+1, cycles, sisi)

def check_primitive(cycle, cycles):
    checked = [i for i in range(3,len(cycle))]
    removed = False
    for cycle_size in checked:
        if removed:
            break
        for primitive in cycles[cycle_size]:
            overlap = len(set(cycle) & set(primitive))
            if (overlap > 2):
                other_paths = [len(cycle) - overlap + 2,
                               cycle_size - overlap + 2]
                if overlap in other_paths:
                    other_paths.remove(overlap)
                    if other_paths[0] < overlap:
                        cycles[len(cycle)].remove(cycle)
                        removed = True
                        break
                elif any(i < overlap for i in other_paths):
                    cycles[len(cycle)].remove(cycle)
                    removed = True
                    break

def ring(datafile, particles, cutoff):
    for i in range(2,cutoff+1):
        try:
            with open(f'primitive_{i}_member_rings.dat', 'w') as rf:
                pass
            with open(f'{i}_member_rings.dat', 'w') as rf:
                pass
        except FileNotFoundError:
            print("Ring lists don't exist. I will create new files")
    sio = []
    #sio_mat = np.zeros((particles,particles))
    linecount = 0
    with open(datafile) as f:
        for line in f:
            linecount += 1
            part = line.split()
            if part[0].isdigit():
                if float(part[1]) == 1:
                    sio.append((int(part[0]), 
                                [int(i) for i in part[3:3+int(part[2])]]))
#print("Si-O connection", sio)
#___________________________________________________________________________________
#All SiSi connections, iD in ascending order (makes it easier for later sorting)
            if linecount % (particles + 7) == 0:
                sisi = []
                cycles = {i: [] for i in range(2,cutoff+1)}
                for i in range(len(sio)):
                    for j in range(i+1, len(sio)):
                        if set(sio[i][1]) & set(sio[j][1]):
                            m1 = set(sio[i][1]) & set(sio[j][1])
                            if len(m1) == 2:    # if 2 Si atoms are connected by 2 O atoms
                                cycles[2].append([sio[i][0], sio[j][0]])
                            sisi.append((sio[i][0],sio[j][0]))
                sio = []
                
                for i in sisi:
                    for j in sisi:
                        if set(i) & set(j) and set(j) != set(i):
                            (m1,) = set(i) & set(j)
                            (rj,) = set(j) - set([m1])
                            (ri,) = set(i) - set([m1])
                            nodes = [ri, m1, rj]
                            counter = 3
                            find_next(nodes, counter, cycles, sisi)
                for i in cycles:
                    print(f"{i} member: {len(cycles[i])}")
                    with open(f'{i}_member_rings.dat', 'a') as rf:
                        for cycle in cycles[i]:
                            rf.write(str(cycle)+'\n')
                        rf.write('\n')
                print('\n')
                for i in range(4,cutoff+1):
                    for cycle in cycles[i][::-1]:
                        check_primitive(cycle, cycles)

                for i in cycles:
                    print(f"{i} member: {len(cycles[i])}")
                    with open(f'primitive_{i}_member_rings.dat', 'a') as rf:
                        for cycle in cycles[i]:
                            rf.write(str(cycle)+'\n')
                        rf.write('\n')
                print(f'\n\n')

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, type=str,
                        help="REAX connectivity file")
    parser.add_argument('-p', '--particles', required=True, type=int,
                        help="number of particles per frame")
    parser.add_argument('-c', '--cutoff', required=True, type=int,
                        help="largest ring size to consider")
    args = parser.parse_args()
    datafile = args.input
    particles = args.particles
    cutoff = args.cutoff

    ring(datafile, particles, cutoff)
