################################################################################
# Name:   SimulationObjects
# Author: Douglas E. White
# Date:   10/17/2013
################################################################################
from simulationMath import *
import random as rand
import math as math

class SimulationObject(object):
    """ Base class from which all simulation obejcts must inherit
    """

    def __init__(self, location, radius, ID, owner_ID, sim_type):
        """ Base class which defines properties all sim objects MUST have
            location - the location of the sphere
            radius - the radius of the sphere
            ID - the ID of the object #WARNING# This ID is used ot hash objects
                  it MUST BE UNIQUE
            owner_ID - usually the same as the ID, also must be unique between
                       agents, unless all agents are part of a larger group
                       i.e. this is the mechanism for multi-agent agents
            sim_type - the type of object the simulation object is
        """
        self.location = location
        self.radius = radius
        self.sim_type = sim_type
        self.ID = ID
        self.owner_ID = owner_ID
        #keep track of the opt and col vecs
        self._disp_vec = [0,0,0]
        self._fixed_contraint_vec = [0,0,0]
        self._v = [0,0,0]
        #keep track of production consumptions values
        self.gradient_source_sink_coeff = dict()
        #keep track of the relative indices in the gradient array
        self.gradient_position = dict()
        #keep track of the value of the gradient associated with these agents
        self.gradient_value = dict()
         
    def update(self, sim, dt):
        """ Updates the simulation object
        """
        pass


    def get_max_interaction_length(self):
        """ Get the max interaction length of the object
        """
        return self.radius*2.0 #in um

    def get_interaction_length(self):
        return self.radius #in um

    def get_spring_constant(self, other):
        """ Gets the spring constant of the object
            Returns: 1.0 by default
            NOTE: Meant to be overwritten by a base class if more
                  functionality is required
        """
        return 0.25

    def add_displacement_vec(self, vec):
        """ Adds a vector to the optimization vector
        """
        self._disp_vec = AddVec(self._disp_vec, vec)

    def add_fixed_constraint_vec(self, vec):
        """ Adds a vector to the optimization vector
        """
        self._fixed_contraint_vec = AddVec(self._fixed_contraint_vec, vec)

    def set_gradient_source_sink_coeff(self, name, source, sink):
        """ Adds a production/consumption terms to the dicationary based on
            the gradient name
        """
        #overwrite exsisting data
        self.gradient_source_sink_coeff[name] = (source, sink)

    def get_gradient_source_sink_coeff(self, name):
        """ Gets the graident terms (source, sink) for the given gradient name
            returns - a tuple of the (source, sink) values. If these are not
                      in the dictionary, returns (0,0)
        """
        #returns the gradient values for the source and sink
        if(name in self.gradient_source_sink_coeff.keys()):
            #name is in the dictionary
            return self.gradient_source_sink_coeff[name]
        else:
            #the name is not in the dictionary
            return (0,0)

    def set_gradient_location(self, name, location):
        """ Adds the location of the agent on the grid for the gradient
            whose name is specified by name
        """
        self.gradient_position[name] = location

    def get_gradient_location(self, name):
        """ Return the location of the agent on the grid for the gradient
            specified by the name, name
        """
        if(name in self.gradient_position.keys()):
            #name is in the dictionary
            return self.gradient_position[name]
        else:
            #the name is not in the dictionary
            return None

    def set_gradient_value(self, name, value):
        """ Adds the value of the gradient at this agent
        """
        self.gradient_value[name] = value

    def get_gradient_value(self, name):
        """ Return the value fo the gradient at this agent
        """
        if(name in self.gradient_value.keys()):
            #name is in the dictionary
            return self.gradient_value[name]
        else:
            #the name is not in the dictionary
            return -1
        
    def update_constraints(self, dt):
        """ Updates all of the contraints on the object
        """
        #first update the posiiton by the col and opt vectors
        #make sure neither of these vectors is greater than error
        mag = Mag(self._disp_vec)
        if(mag > 5):
            n = NormVec(self._disp_vec)
            self._disp_vec = ScaleVec(n, 5.0)
        self.location = AddVec(self.location, self._disp_vec)
        #then clear it
        self._disp_vec = [0,0,0]
        
        #then update the the pos using the fixed vectors
        mag = Mag(self._fixed_contraint_vec)
        if(mag > 5):
            n = NormVec(self._disp_vec)
            self._fixed_contraint_vec = ScaleVec(n, 5.0)
        self.location = AddVec(self.location, self._fixed_contraint_vec)
        #htne clear it
        self._fixed_contraint_vec = [0,0,0]
        
    def __repr__(self):
        """ Returns a string representation of the object
        """
        return self.sim_type+": "+repr(self.ID)+" "+repr(self.location)

    def __eq__(self, other):
        """ Handles the equal operator for the object
        """
        if(isinstance(other, SimulationObject)):
            return self.ID == other.ID
        #otherwise
        return False

    def __hash__(self):
        """ Handles the hashing operator for the object
        """
        return hash(self.ID)

class DividingCell(SimulationObject):
    """ A stem cell class
    """
    def __init__(self, location, radius, ID,state="U",params=None,
                 division_set = 0.0,
                 division_time = 14.0,
                 owner_ID = None):
        """ Constructor for a stem cell
            location - the location fo the stem cell
            radius - the size of the stem cell
            ID - the unique ID for the agent
            state - the state of the stem cell
            division_set - the initial division set for the cell
            division_time - the time it takes the cell to divide
            owner_ID - the ID associated with the owner of this agent
        """
        #define some variables
        if(owner_ID == None):
            owner_ID = ID
        self.division_timer = division_set
        self.division_time = division_time
        self.state = state #Demarcus
        #call the parent constructor
        super(DividingCell, self).__init__(location,
                                       radius,
                                       ID,
                                       owner_ID,
                                       "dividingcell")
        self.params = params  #Demarcus, save params

    def update(self, sim, dt):
        """ Updates the stem cell to decide wether they differentiate
            or divide
        """
        #growth kinetics
        #print("Test: Updating Dividingcell class")
        self.division_timer += dt
        if(self.division_timer >= self.division_time):
            #now you can divide
            #get the location
            #pick a random point on a sphere
            #TODO-Demarcus: Make sure cell is dividing.
            print("Test: Stem cell is dividing")
            location = RandomPointOnSphere()*self.radius/2.0 + self.location
            #get the radius
            radius = self.radius
            #get the ID
            ID = sim.get_ID()
            #make the object
            sc = DividingCell(location, radius, ID)
            #add it to the imsulation
            sim.add_object_to_addition_queue(sc)
            #reset the division time
            self.division_timer = 0


    def get_interaction_length(self):
        """ Gets the interaciton elngth for the cell. Overiides parent
            Returns - the length of any interactions with this cell (float)
        """
        return self.radius + 2.0 #in um
    
class StemCell(SimulationObject):
    """ A stem cell class
    """
    def __init__(self, location, radius, ID, state,
                 params = None,
                 division_set = 0.0,
                 division_time = 19.0,
                 owner_ID = None,
                 differentiation_rule = 'competing_feedback_2015'):
        """ Constructor for a stem cell
            location - the location fo the stem cell
            radius - the size of the stem cell
            ID - the unique ID for the agent
            state - the state of the stem cell
            division_set - the initial division set for the cell. cells are not synchronized.
            division_time - the time it takes the cell to divide
            owner_ID - the ID associated with the owner of this agent
        """
        #define some variables
        if(owner_ID == None):
            owner_ID = ID
        #set thet state
        self.state = state
        self.division_timer = division_set
        self.division_time = division_time
        #call the parent constructor
        super(StemCell, self).__init__(location,
                                       radius,
                                       ID,
                                       owner_ID,
                                       "stemcell")
        #save the parameters
        self._params = params
        
    def update(self, sim, dt):
        """ Updates the stem cell to decide wether they differentiate
            or divide
            
        returns: 1 if in D state, 0 if in T or U state. For convergence test.
        """
        #growth kinetics
        self.division_timer += dt
        #you can grow unless you are in the A state meaning apoptosis
        if(self.division_timer >= self.division_time):
            #check state
            #now you can divide
            if(self.state == "T"):
                #change the current sytate to D
                self.state = "D"
                self.division_time = 51 #in hours
            
            #get the location, random point on a sphere
            location = RandomPointOnSphere()*self.radius/2.0 + self.location
            radius = self.radius  #get the radius
            ID = sim.get_ID()     #get the ID
            
            #make the object
            sc = StemCell(location, radius, ID, self.state,
                          params = self._params,
                          division_time = self.division_time)
            #add it to the simulation, Demarcus
            sim.add_object_to_addition_queue(sc)
            #reset the division time
            self.division_timer = 0
        
        if(self.state == "U"):
            #then the stem cell is still a stem cell
            #HANDLE DIFFERENTIATION

            #RANDOM RULE
            x = rand.random()
            prob = self._params[0]
            #longer before the differentiation starts
            if(x < prob):
                #differentiation occurs
                self.state = "T"

            #get the neighboring states
            nbs = sim.network.neighbors(self)
            tot = len(nbs)
            if(tot > 0):
                u = 0
                d = 0
                for i in range(0, tot):
                    if(nbs[i].state == "U" or nbs[i].state == "T"):
                        u += 1
                    if(nbs[i].state == "D"):
                        d += 1
                        
                #set the hill coefficients
                knf = self._params[1]
                kpf = self._params[2]
                n_p = self._params[3]
                n_n = self._params[4]
                
                #now calcualte the values
                norm_u = float(u) / float(tot)
                norm_d = float(d) / float(tot)
                
                ### Demarcus 07292015. Add 2013 Competing Feedback Rule
                #x1 = rand.random()
                #competing_feedback = 1. / (1. + math.expm1(float(u)-float(d)))
                #if (x1 < competing_feedback):
                #    self.state = "T"
                
                #use as a ratio of u to d
                #compute the activator porbability
                negative_feedback = 1. / (1. + (norm_u/knf)**n_n)
##                negative_feedback = 0
                #compute the inhibitor probability
                postive_feedback = (1*norm_d**n_p)/(kpf**n_p + norm_d**n_p)
##                postive_feedback = 0
                
                #Now compute them seperately
                x1 = rand.random()
                x2 = rand.random()
                
                #IMPORTANT: normally this is an OR not AND
                if(x1 < negative_feedback or x2 < postive_feedback):
                    #put yourself into a differentiating state
                    self.state = "T"
                    
        #Convergence updating statistic. Demarcus
        # need to know what percent of cells are in D state.
        # not limited by time_step of cell count
        if(self.state == "D"):
            return 1
        else:
            return 0


    def get_interaction_length(self):
        """ Gets the interaciton elngth for the cell. Overiides parent
            Returns - the length of any interactions with this cell (float)
        """
        return self.radius + 2.0 #in um

