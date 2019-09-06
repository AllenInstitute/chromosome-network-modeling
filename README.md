chromosome-network-modeling
===========================

[Dev Branch]

This repository contains the scripts needed to reproduce the analysis published 
in *Parental and stochastic skewing of X-chromosome inactivation predicts 
phenotypic penetrance of fragile X syndrome* (citation link to come)

To reproduce this work:
1. Clone this repository
```
git clone https://github.com/AllenInstitute/chromosome-network-modeling.git
```
2. Install requirements (example uses `virtualenv`, one could also use `conda` 
	 or other enviroment manager):
	- [Install 
		virtualenv](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/#installing-virtualenv) 
		if necessary:
	 ```
	 python3 -m pip install --user virtualenv
	 ```
	- Create a new environment
	 ```
	 cd chromosome-network-modeling
	 python3 -m venv env
	 ```
	- Activate and install dev requirements
	 ```
	 source env/bin/activate
	 pip install -r requirements-dev.txt
	 ```
	- Install mcmodels
	 ```
	 pip install git+https://github.com/AllenInstitute/mouse_connectivity_models
	 ```
2. To reproduce the analysis, run the following scripts:
	- Run scripts to produce connectivity and distance matrices
	 ```
	 python3 build_connectivty_model.py
	 python3 build_full_regional_model.py
	 python3 get_region_distances.py
	 ```
	- Run script to reproduce analysis
	 ```
	 python3 run_test_distance_sample.py
	 ```
