import numpy as np


class population:
    """defines a population with n_members individuals. The fitness function is set at the population level so that all individuals
    share the same fitness function. The fitness function is to be minimized, so the most fit individual will have the lowest fitness
    value. """
    
    def __init__(self, n_members: int, bounds: float, fitnessFunction: function):
        self._nMembers = n_members
        self._generation = 0
        
        # individuals are generated with the pop_member initializer
        self._members = [pop_member([np.random.random()*(bounds[j][1] - bounds[j][0]) + bounds[j][0] for j in range(len(bounds))]) for i in range(n_members)]
        self._fitnessFunc = fitnessFunction
        self._bounds = bounds
        for i in self._members:
            i._fitness = self._fitnessFunc(i.params)
        
        # members are sorted so that the least fit member is easily accessed
        self._members.sort(key=lambda member: member.fitness)
        
    @property
    def bounds(self):
        return self._bounds
    
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
        
        self._members[idx].mutate(self._bounds, 1.0/(self._generation+1))
        self._members[idx]._fitness = self._fitnessFunc(self._members[idx].params)
        self._members.sort(key=lambda member: member.fitness)

        
    def mate(self):
        """Choose two individuals from the population and create a child. If the child is fitter than the
        least fit current individual in the population, replace the individual with the child"""
        
        parents = np.random.randint(self._nMembers, size=2)
        child = pop_member(np.choose(np.random.randint(2, size=len(self._bounds)), [self._members[parents[0]].params, 
                                                                                            self._members[parents[1]].params]))
        #Mutate the child to widen parameters outside the parents
        child.mutate(self._bounds, 1.0/(self._generation+1))
        child._fitness = self._fitnessFunc(child.params)

        if (child.__lt__(self._members[-1])):
            self._members[-1] = child
        
        self._members.sort(key=lambda member: member.fitness)
    
    def evolve(self):
        """Move the population forward a generation by creating population.nMembers new potential children.
        Whether children survive is determined by the population.mate() function."""
        for whichmember in range(self.nMembers):
            self.mate()
        self._generation += 1


class pop_member:
    """Define an individual. Each individual has a parameter vector and a fitness."""
    
    def __init__(self, params):
        self.params = np.array(params)
        self._fitness = 0.

    def __gt__(self, othermember):
        """Define greater-than comparison function."""
        return self.fitness > othermember.fitness
    
    def __lt__(self, othermember):
        """Define less-than comparison function."""
        return self.fitness < othermember.fitness
        
    def mutate(self, bounds, scale):
        """Mutate a member by choosing a single parameter and modifying it. If the new
        parameter value falls outside the bounds, randomly assign a value within the bounds."""
        idx = np.random.randint(len(self.params))
        bounds = bounds[idx]
        self.params[idx] = self.params[idx] + (np.random.random() - 0.5)*(bounds[1] - bounds[0])*scale
        if (self.params[idx] < min(bounds) or self.params[idx] > max(bounds)):
            self.params[idx] = np.random.random()*(bounds[1] - bounds[0]) + bounds[0]
            
    @property
    def fitness(self):
        return self._fitness
    