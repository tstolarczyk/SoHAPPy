Building the documentation 
Install Sphinx and the extensions. 
Extensions should be installed separately. As an example the automodapi, that create an output superior to the autodoc standard package, has to be installed (see here): 
conda install -c conda-forge sphinx
conda install -c conda-forge sphinx-automodapi
conda install –c conda-forge nbsphinx
conda install -c conda-forge recommonmark
conda install -c conda-forge numpydoc
conda install -c conda-forge sphinx_rtd_theme

make html


Todo list
#########

Various
=======

- Stack the prompt with the afterglow
- Can we suppress simulation for large S and B. yes, but how well ?
- Inject Fermi data and IRF, compute response
- Interpolate time counts within the time slices;
- Use papermill to have the jupyter notebbok automatically run;
- Use jupyter nbconvert to convert executed jupyter notebook into html files;


Masking the first dataset in a stacked Datasets
===============================================

*From a Gammapy Slack conversation with @adonath, Jan. 2021*

By default the :obj:`mask_safe` is not applied to a dataset so that it can be changed for different analysis with the same data. 

When stacking a set of Dataset, all datasets except the first one are masked.
Up to gammapy 0.17 at least, in order to have the mask applied it had to be done "by hand" with the following function:

.. code-block::

	def get_masked_dataset(ds0):
		"""
		This is a trick to get a masked datset from a dataset.

		Parameters
		----------
		ds0 : Dataset object
			The original complete unmasked dataset.

		Returns
		-------
		masked_dataset : Dataset object
			The original dataset with masking applied.

		"""
		# print("dataset_tools/get_maked_dataset(ds) : masking disabled")
		# return ds0

		e_true = ds0.aeff.energy.edges
		e_reco = ds0.counts.geom.axes[0].edges
		region = ds0.counts.geom.region
		masked_dataset = SpectrumDataset.create(
			e_true=e_true, e_reco=e_reco, region=region, name=ds0.name
		)
		masked_dataset.models = ds0.models
		masked_dataset.stack(ds0)

		return masked_dataset

This functionnality is planned to be released soon as :func:`to_masked()` in `dataset/map.py` and should replace the existing function.

