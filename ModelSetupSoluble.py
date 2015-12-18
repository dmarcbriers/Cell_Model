from __future__ import division
import sys
import os
import random as r
import time
import argparse

from scipy.spatial.distance import euclidean
import networkx as nx

from StemCellSoluble import *
from Simulation import *
from Gradient import Gradient
#import simAnaylsisScriptOct4Nanog as sa
import platform
import simulationUtils as su

#runs from the command line and takes arguments in the follwing format
#-0
#-1 structure sim save file path
#-2 simulation base save file path
#-3 simulation number
#-4 start time
#-5 end time
#-6 time step

if(__name__ == '__main__'):
    
    ###############################
    ### Command Line Arguments ####
    ###############################

    #Use argparse Module to handle command line arguments and gi
    parser = argparse.ArgumentParser(description=\
    "This is an agent-based model that can explore spatial differentiaion patterns\n" + \
    "in aggregates of stem cells known as Embroid Bodies. See White et al Integrative Biol (2015).")
    
    #Set command line arguments
    parser.add_argument("sim_id",
                        help="There are several starting simulation structures numbered 0-99 to " +\
                             "capture biological variability. Choose any structure numbered 0-99.",
                        type=str)
    parser.add_argument('a',help='basal differentiaion probability. Range: .01,.005,.002',type=float)
    parser.add_argument('k1',help='tuning paramter for negative feedforward rule. Range: [0,1]',type=float)
    parser.add_argument('n1',help='Hill coefficient for negative feedforward rule. Range: 10,25,50. Max ~300',type=float)
    parser.add_argument('k2',help='Tuning parameter for positive feedback rule. Range: [0,1]',type=float)
    parser.add_argument('n2',help='Hill coefficient for positive feedback rule. Range: 10,25,50. Max ~300',type=float)
    parser.add_argument('time_end',help='Maximum number of time steps the model will run.(Default=200)',
                        nargs='?',type=int,default=200)

    args = parser.parse_args()


    #struct_path = str(args[1])
    #sim_base_path = str(args[2])
    #sim_id = int(args[3])
    #ts = float(args[4])
    #te = float(args[5])
    #dt = float(args[6])
    #p1 = float(args[7])
    #p2 = float(args[8])
    #p3 = float(args[9])

    #Keep track of runtime
    start_time = time.time()

    ts = 0
    te = 100
    dt = 1
    sim_id = args.sim_id #'9'
    #sim_base_path = '/home/dbriers/hyness/EBICS/Cell_Model' #"C:\\users\\doug\\desktop\\"
    #struct_path = os.path.join(sim_base_path,"1000",str(sim_id))
    #struct_path = sim_base_path + "structs\\EB\\1000\\25\\"
    #p1 = args.a  #.001
    #p2 = args.k1 #.06
    #p3 = args.n1 #10
    #p4 = args.k2 #.1
    #p5 = args.n2 #100 #p5 > 340 needs type bigFloat,Decimal. #1000.
####    p1 = .001
##    p1 = 0.000
####    p2 = .007
####    p2 = 0.007
##    p2 = 0.0075
##    p3 = 10
##    p4 = 5E-4
##    p5 = 10.0
    
    #OUTPUT Folder as relative directory.
    list_of_params = ["results",str(args.a),str(args.k1),str(args.n1),str(args.k2),str(args.n2)]
    sim_base_path = '_'.join(list_of_params)
    struct_path =  os.path.join("1000", sim_id) ##TODO: Take from argparse.
    
    # Output folder and LOG file.
    su.make_sure_path_exists(os.path.join(sim_base_path,sim_id))
    log_file = open( os.path.join(sim_base_path,sim_id,'simulation.log') ,"w")
    sys.stdout = log_file
 

    #################################
    ### Simulation Configurations ###
    #################################

    #Now make a new simulation
    sim = Simulation(sim_id,
                     sim_base_path,
                     ts,
                     te,
                     dt)
    #add some gradients
##    #from Van Winkle et al 2012 for 02
##    D = 1.7e-9 / (1e-12) #um^2/sec
##    consump_rate = (4.0e-17)#*(10**6) mols/cell/sec #umols/cell/sec
##    #This value later gets normalized using the 
##    #set the outside concentration
##     outside_c = 1.04e-4*10**6 #(umol/L) or uM
    #from War-tenberg et al., 2001

    D = 10.0 # um^2/sec
    consump_rate = 2.0e-20 #mol/cell/sec
    production_rate = 2.0e-20
##    outside_c = .1 #umol/L or uM
    outside_c = 0.0 #umol/L
    tnf = Gradient("TNF", D, 700.0, 700.0, 700.0, 15, 15, 15,
                 outside_c = outside_c)
    #add this to the simulation
    sim.add_gradient(tnf)
    lif = Gradient("LIF", D, 700.0, 700.0, 700.0, 15, 15, 15,
                 outside_c = outside_c)
    #add this to the simulation
    sim.add_gradient(lif)
##    while(True):
##        print(struct_path)
##        print(sim_base_path)
##        print(sim_id)
    #find the only file in the dir
    f = os.listdir(struct_path)
##    while(True):
##        print(struct_path + f[0])
    #then load this file
    #open the input structure file
    net = nx.read_gpickle(os.path.join(struct_path,f[0]))
    cent = (0,0,0)
    placed = False
    nodes = net.nodes()
    np.random.shuffle(nodes)
    id_map = dict()
    for i in range(0, len(nodes)):
        node = nodes[i]
        #randmoize the division set
        #so that cells divide at different times
        div_set = r.random()*19
        #now make a new stem cell
        sim_obj = StemCell(node.location, node.radius, node.ID,
                           "U", division_set = div_set,
                           params = [args.a, args.k1, args.n1, args.k2, args.n2])

 
        #add some gradient values
        sim_obj.set_gradient_source_sink_coeff(lif.name,
                                               production_rate,
                                               0)
        sim_obj.set_gradient_source_sink_coeff(tnf.name,
                                               production_rate,
                                               production_rate/2.)
        #keep track fo the ID dict for the connection mapping
        id_map[node.ID] = sim_obj
        #add it to the sim
        sim.add(sim_obj)

    #also import the connections to use as well
    cons = net.edges()
    for i in range(0, len(cons)):
        con = cons[i]
        n1, n2 = con
        s1 = id_map[n1.ID]
        s2 = id_map[n2.ID]
        sim.network.add_edge(s1, s2)
    #Run the simulation
    sim.run()


    ################################
    
    #Print runtime of the simulation.
    print 'Model took', (time.time()-start_time)/60, 'minutes.'
    log_file.close()

    
    #########################
    ### analyze the data ####
    #########################

    #print(sim_base_path)
    #sa.get_simulation_metric_data(sim_base_path + os.sep + repr(sim_id) + os.sep)

