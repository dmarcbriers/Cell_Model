###### This guide will explain how to run the model.
##### Background info:

See white et al 2013 for motivation of spatial patterns
See white et al 2015 integrative biology for descriptio of current model

The model consists of two parts.
1. Probablistic model stored as a networkX graph network.
2. A visualizaiton using vpython/matplotlib, for visual pattern recognition.


#### REQUIRED SOFTWARE
python 2.7+
## Required Python Modules for Simulation
scipy
numpy
networkx
## Required Python Modules for visualization
visual(vpython) - for making 3D images.
PIL - for capturing 3D images made by vpython.
matplotlib - for making 2D images.


#####
### Running the probablistic model of local cell-cell interactons:
#####
    #1A. Single simulations
### ModelSteup.py is the file you need to run. Here are the parameters:

'python ModelSetup.py SIM_ID a k1 n1 k2 n2'

  SIM_ID - write a bash script to pass the number 1-99 to the model. You can run the model
	multiple time with the same set of parameters.
  a - random differentiation [.01, .001]
  k1 - normalization param for negative feedback. ranges from 0-1. [0.1, .5, .9]
  n1 - hill coefficient for negative feedback. Try 1-100? [10, 50]
  k2 - normalization param for positive feedback. ranges from 0-1 [.1, .5, .9]
  n2 - hill coefficient for posititve feedback. Try 1-100? [10, 50]


   #1B. Multiple simulaitons.
   To run multiple simulations, modify the bash script 'run_parallel_simulations.sh'. Make
   sure you have the gnu utility 'parallel' installed. This runs REPLICATES of a simulation.


   If you are using a cluster you could also launch a job array to queue simulations.

### Example parameters.
python ModelSetup.py $SIM_ID .01 .3 25 .5 25


#### Extended Example of parameter ranges:
a - [.001, .01] recommended: .001,.005, or .01
k1 - [0,1] recommended: 0.1,0.3,0.6,0.9
n1 - [1,100] recommended: 25(white 2013) 10,50(White 2015)
k2 - [0,1] recommended: range of [.7,1]. but White explored [0.1,0.3,0.6,0.9]
n2 - [1,100] recomended: 25(White 2013) 10,50(White 2015)

