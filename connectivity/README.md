
To run most scripts in this directory, set your local path to the package as local_path and run startup.py and analysis.py first

required libraries: pymaid, python-igraph, numpy, pandas, matplotlib, scipy, seaborn, etc.

- Reconstructed neurons in the form of skeletons are in *data/compact_sk/* ([field description](https://catmaid.readthedocs.io/en/stable/_static/api/index.html#operation/skeletons_compact-detail_list))

- Connectivity matrix can be set up by running *analysis.py*. Alternatively, it can be extracted from *data/pre_post_info/* ([field description](https://github.com/bocklab/pn_kc/blob/master/mushroom_2to3/connect_path.py#L265)).

-  Skeleton ids for different groups of neurons are in *data/skids*
