Adding new cell types
=====================

Here are the steps to add a new cell type:
    
1. Create a source python file in cnmodel/cells, likely based on the bushy.py source, renaming variables as necessary. The main routines in the class however, should maintain their present names and calling parameters.
    
1. Add the values for the new cell type to the data tables (in cnmodel/data). Note that many tables will need to be updated with new columns for the new cell type. Use the existing tables when possible; add new tables when necessary. This step can be quite difficult if you do not know what kinds of ion channels or synaptic currents need to be incorporated into the model.

1. Run the model and make sure the new cell type is performing as desired. You can extend the examples in the examples directory to handle the new files. For example, plot the IV curve using test_cells. The cells should be verified that they meet key target parameters that have been identified in advance. 

1.  Once the model is performing as expected, update the unit tests to include the new cell type. The unit tests will first fail; run with the --audit flag to show the simulations and results parameters. You can then "commit" the results of the simulation to the database, so that they will be used as the benchmark in future tests.



