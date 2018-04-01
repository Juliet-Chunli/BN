# quantifyBN
A platform to analyze nuclear power plant safety by solving Bayesian Network dynamically.

Major computation work is done in computeGraph.py which generates a graph based on input structure and supports its quantification as long as all conditional dependency relationships are specified as required. 

node.py defines node class in computeGraph which has attributes corresponding to its own conditional and non-conditional probability distribution. The graph structure is defined by adding child and parent for each node instance.

event_tree.py a class to quantify event trees which is specified to nuclear power plant application.

timeDepTemplate.py is a template for time dependent scenario quantification which could be run and get output the change of variables of interest according to time.

All others are support function or input data for conditional probability distribution. 
