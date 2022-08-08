Population analysis
===================

Introduction
------------
The population is analysed from the devoted SoHAPPy output file (named `data.txt` by default).

The available tools consist in a few jupyter notebooks and python script files.
The user is intended to indicate the SoHAPPy folder containing the population file in a file named `parameter.yam`, with the variable "outfoler", like shown here:
..  code-block::

	# outfolder : "../../../output/pop_vis24_moonlight/" 
	outfolder : "../../../output/pop_vis24_fullmoonveto/" 

The execution of Jupyter notebook will then present various results as follows:

- *statistics.iynb* :

	* Open a file (`get_data`), build the GRB population for North, South, North only, South only or both sites, performs internal checks.
	
	open a file and derive the 3 and 5 sigma detections rates. 
	* Extract from the data a few characteristics and check that they are compatible with the parameters in the confguration file (delay after the GRB trigger, minimal allowed alitude for detection)
	* Give the list of GRB detected more than 24 h after the Trigger.

- Another one
	* To be described

   
.. automodapi:: analysis.population.init
   :no-inheritance-diagram:

.. automodapi:: analysis.population.pop_plot
   :no-inheritance-diagram:

.. automodapi:: analysis.population.setup
   :no-inheritance-diagram: