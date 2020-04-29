Building Simulations with CNModel
=================================

CNModel is meant to be used as an imported package under Python. This means that it does not provide stand-alone funcitonality on its own. The files in the examples directory demonstrate how this the package is meant to be used for different levels and kinds of simulation, as well as providing tests and demonstrations of components.

Typically, we create a separate directory (a "simulation" directory) that holds the code that uses cnmodel for simulations. This directory will be at the same directory level as cnmodel or possibly elsewhere. Do not place the simulation directory inside cnmodel.

This is a minimal example for a bushy cell current-voltage relationship

Importing
---------

The usual setup for a model will start with::
    
    import sys
    import pyqtgraph as pg  # graphics, requires pyqt
    import cnmodel
    from cnmodel.protocols import IVCurve, VCCurve

Create a point cell
-------------

Next, we will set some parameters and instantiate a cell::
    
    temperature = 22. # degrees C
    cell = cnmodel.cells.Bushy.create(model='RM03', species='guineapig, modelType='II',
        ttx=False, nach=None)
    cell.set_temperature(float(temperature))  # confirm temperature is valid for this model
        
Run and plot a simulation
-------------------------

Then we can set up a simulation and run it::

    V0 = cell.find_i0(showinfo=True)
    iv = IVCurve()
    iv.run({'pulse': [(-1, 1.2, 0.05)]},
           cell, durs=[10, 100, 20], 
           sites=None, reppulse=None, temp=float(temperature))
    result = iv.input_resistance_tau()
    print('    From IV: Rin = {:7.1f}  Tau = {:7.1f}  Vm = {:7.1f}'.format(result['slope'], result['tau'], result['intercept']))
    iv.show(cell=cell)

    if sys.flags.interactive == 0:  # keeps application alive for display of data
        pg.QtGui.QApplication.exec_() 


Simulating cells with morphology defined in swc or hoc files
============================================================

For more biophysically realistic stimulations, a cable representation of a cell may be used and decorated with channels. The density of channels in different compartments (or compartment types, for example, axon, axon hillock, dendrite, etc) is defined in the data/channels.py. All that is necessary is to add the name of the .hoc file to the call that creates the cell (as above), using the morphology keyword. The names in the hoc file must match exactly the compartment type names in the data table, and the "soma" compartment should have a lower-case name. Note that it is easy with the tables to more finely decorate compartments if the reconstructions include the appropriate morphological information. This should be done as described below in the instructions for defining simulation parameters.

A tool for converting swc files from reconstructions to hoc files is included in neuronvis (swc_to_hoc). This tool also comments the hoc files with the number of each swc segment, so that references to the original reconstruction are available.

Defining Simulation Parameters with CNModel
===========================================

The construction of single cell and network models with CNModel is designed to be controlled by text tables. The default set of tables are in cnmodel/data. There are tables for synapses, ionchannels, populations, connectivity. These tables should not be directly modified. If it is desired to change the parameters specified in these tables, it is best to copy those you wish to modfiy into your local model directory in a separate path, and modify them there.

1. The data tables are very strict about column alignment. The first character of the column title and the each of the values in that column must line up directly. It is best/easiest to edit these tables in a programming editor with fixed width fonts and the ability to perform column-based insertions. Changes to the data tables should be annotated appropriately.

2. Channels and receptors are all specified as NEURON .mod files. Some of the parameters in these files are available for external manipulation (in all cases, the maximal conductance "gbar" is available; in some models other parameters are exposed). Many other parameters, such as rate constants for the receptors and ion channels, are specified in the .mod files and are not exposed externally.

3. The connectivity data table can be modified to represent a particular pattern of connectivity, and the populations data table can be modified to change the relative numbers of cells used to build a network. Note that only those cells actually needed according to the connectivity patterns will be instantiated.

4. The tables are set up so that they can be annotated by folllowing a parameter with a [n]. This allows you to keep track of why a parameter was changed, or the provenance of a different value. In exploratory simulations, it might be a good idea to generate different tables with appropriate annotations for different exploratory values.

Adding new mechanisms to a cell requires modification of the code that instantiates a cell to recognize the mechanisms at several points. It is especially important to handle new channel types in cnmodel/cells.py, where knowledge of channel names is needed to compute initial states and resting potentials. In addition, channels must be referenced in the code that creates a given cell type; this occurs where the channels are actually inserted. Channels will also need to be referenced in the relevant data tables. Specific naming conventions should be followed to simplify integration. Although the process can be straightforward, contact the authors for help if you have questions.


Data tables that represent new versions of cells or which incorporate additional mechanisms can then be used in sumulations by incorporating the following code pattern::

        from cnmodel import data
        import data_XM13nacncoop as CHAN  # where data_XM13nacncoop.py is a modified table in the current imulation directory
        # The following takes the table named "XM13nacncoop_channels" from the CHAN table,
        # and overwrites the original table "XM13nacncoop_channels" that is in ionchannels.py. The original file in cnmodel is
        # not modified, only the data in memory. 
        changes = data.add_table_data('XM13nacncoop_channels', row_key='field', col_key='model_type',
                       species='mouse', data=CHAN.ChannelData)
        # The following takes the table indicating how the channel compartments should be decorated from the ChannelComparments
        # table, overwriting the original named table in ionchannels.py
        changes_c = data.add_table_data('XM13nacncoop_channels_compartments', row_key='parameter', col_key='compartment',
                species='mouse', model_type='II', data=CHAN.ChannelCompartments)
        # now print out what was changed!
        data.report_changes(changes)
        data.report_changes(changes_c)

We recommend performing these changes as tests with a minimal model first to simplify debugging.


