# ToDo:
# - redesign evolution algorithm to create a new population
#   without replacing any of the old population. Then, take
#   top 50% of parent and child populations
# - add threading for creating child population
# - add convergence criteria
#

import numpy as np
from abc import ABC, abstractmethod
import xml.etree.ElementTree as et
import copy

class population:
    """defines a population with n_members individuals. The fitness function is set at the population level so that all individuals
    share the same fitness function. The fitness function is to be minimized, so the most fit individual will have the lowest fitness
    value. """
    
    def __init__(self, bounds, ffxml, fitnessFunction):
        self._generation = 0
        # individuals are generated with the _member initializer
        self._initialff = ffxml
        self._fitnessFunc = fitnessFunction
        self._bounds = bounds

    @classmethod
    def new_population(cls, n_members: int, bounds, ffxml, fitnessFunction):
        instant = cls(bounds, ffxml, fitnessFunction)
        instant._nMembers = n_members
        instant._generation = 0
        # individuals are generated with the _member initializer
        instant._members = [_member.make_member(bounds, ffxml, False) for i in range(n_members)]
        for i in instant._members:
            i.fitness = instant._fitnessFunc(i.params)
        
        # members are sorted so that the least fit member is easily accessed
        instant._members.sort(key=lambda member: member.fitness)
        return instant

    @classmethod
    def read_population(cls, bounds, base_xml, ffxml, fitnessFunction):
        instant = cls(bounds, base_xml, fitnessFunction)
        # individuals are generated with the _member initializer
        instant._members = [_member.make_member(bounds, i, True) for i in ffxml]
        instant._nMembers = len(instant._members)
        for i in instant._members:
            i.fitness = instant._fitnessFunc(i.params)
        
        # members are sorted so that the least fit member is easily accessed
        instant._members.sort(key=lambda member: member.fitness)
        return instant

    @property
    def bounds(self):
        return self._bounds

    @property
    def initialff(self):
        return self._initialff
    @property
    def generation(self):
        return self._generation
    
    @property
    def nMembers(self):
        return self._nMembers
    
    @property
    def members(self):
        return self._members
    
    def mutate_member(self, idx):
        """Choose a member and mutate it."""
        
        self._members[idx].mutate(self._bounds,1.0)
        self._members[idx]._fitness = self._fitnessFunc(self._members[idx].params)
        self._members.sort(key=lambda member: member.fitness)

        
    def mate(self):
        """Choose two individuals from the population and create a child. If the child is fitter than the
        least fit current individual in the population, replace the individual with the child"""
        
        parents = np.random.randint(self._nMembers, size=2)
        child = _member.make_member(self.bounds, self.initialff, True)
        root = self.bounds.getroot()
        for section in root:
            s_tag = section.tag
            for interaction in section:
                i_tag = interaction.tag
                i_attrib = f'[@order="{interaction.get("order")}"]'
                for group in interaction:
                    parent = self.members[parents[np.random.randint(2)]].params
                    where = f'./{s_tag}/{i_tag}{i_attrib}/{group.tag}[@order="{group.get("order")}"]'
                    param = parent.find(where).text
                    child.params.find(where).text = param

        #Mutate the child to widen parameters outside the parents
        if np.random.random_sample() < 0.6:
            child.mutate(self._bounds, 1.0)
        child._fitness = self._fitnessFunc(child.params)
        return child
        # if (child.__lt__(self._members[-1])):
        #     self._members[-1] = child
        
        # self._members.sort(key=lambda member: member.fitness)
    
    def evolve(self):
        """Move the population forward a generation by creating population.nMembers new potential children.
        Whether children survive is determined by the population.mate() function."""
        children = []
        for whichmember in range(self.nMembers):
            children.append(self.mate())
        children.sort(key=lambda member: member.fitness)
        self._members = self._members[:self.nMembers//2] + children[:(self.nMembers+1)//2]
        self._members.sort(key=lambda member: member.fitness)
        self._generation += 1


class _member(ABC):
    """Define an individual. Each individual has a parameter vector and a fitness."""
    
    def __init__(self):
        self._fitness = 0.

    def make_member(bounds,ffxml,restart):
        creator = choose_type(bounds)
        return(creator(bounds,ffxml,restart))
    
    def __gt__(self, othermember):
        """Define greater-than comparison function."""
        return self.fitness > othermember.fitness
    
    def __lt__(self, othermember):
        """Define less-than comparison function."""
        return self.fitness < othermember.fitness
        
    @property
    def fitness(self):
        return self._fitness

    @fitness.setter
    def fitness(self, value):
        self._fitness = value
    
    @abstractmethod
    def mutate(self):
        pass
   

# class param_list_member(_member):
#     """Parameters are given as a simple list or numpy array"""
# 
#     def __init__(self, bounds):
#         super().__init__()
#         self._params = []
#         for i in range(len(bounds)):
#             self._params.append(bounds[0] + np.random.random()
#                                 *(bounds[1] - bounds[0]))
# 
#     def mutate(self, bounds, scale):
#         """Mutate a member by choosing a single parameter and modifying it. If the new
#         parameter value falls outside the bounds, randomly assign a value within the bounds."""
#         idx = np.random.randint(len(self.params))
#         bounds = bounds[idx]
#         self.params[idx] = self.params[idx] + (np.random.random() - 0.5)*(bounds[1] - bounds[0])*scale
#         if (self.params[idx] < min(bounds) or self.params[idx] > max(bounds)):
#             self.params[idx] = np.random.random()*(bounds[1] - bounds[0]) + bounds[0]
            

class reaxFF_member(_member):
    """Parameters are defined by a nested dictionary and xml representation of
    the forcefield file"""

    def __init__(self, bounds, ffxml, restart):
        super().__init__()
        self._chromosome = []
        self.params = copy.deepcopy(ffxml)
        root = bounds.getroot()
        for section in root:
            s_tag = section.tag
            for interaction in section:
                i_tag = interaction.tag
                i_attrib = f'[@order="{interaction.get("order")}"]'
                for group in interaction:
                    where = f'./{s_tag}/{i_tag}{i_attrib}/{group.tag}[@order="{group.get("order")}"]'
                    self._chromosome.append(where)
                    if not restart:
                        lower = float(group.get("lower"))
                        upper = float(group.get("upper"))
                        param = lower + np.random.random()*(upper - lower)
                        self.params.find(where).text = str(param)
        
    @property
    def chromosome(self):
        return self._chromosome

    def mutate(self, bounds, scale):
        """Mutate a member by choosing a single parameter and modifying it. If the new
        parameter value falls outside the bounds, randomly assign a value within the bounds."""
        gene = self.chromosome[np.random.randint(len(self.chromosome))]
        root = bounds.getroot()
        lower = float(root.find(gene).get('lower'))
        upper = float(root.find(gene).get('upper'))
        mutation = (np.random.random()-0.5)*(upper - lower)*scale
        param = float(self.params.find(gene).text)
        newparam = param + mutation
        if newparam < lower or newparam > upper:
            newparam = lower + np.random.random()*(upper - lower)
        self.params.find(gene).text = str(newparam)


def choose_type(bounds):
    if type(bounds) == list:
        return param_list_member
    elif type(bounds) == et.ElementTree:
        return reaxFF_member
    else:
        raise ValueError(f'Bounds are of type {type(bounds)}, which is not of a valid type')

