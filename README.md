chromosome-network-modeling
===========================

This repository contains the scripts needed to reproduce the analysis published 
in *Parental and stochastic skewing of X-chromosome inactivation predicts 
phenotypic penetrance of fragile X syndrome*

To reproduce this work:
1. Clone this repository
```
git clone https://github.com/AllenInstitute/chromosome-network-modeling.git
```
2. Install requirements (example uses `virtualenv`, one could also use `conda` 
	 or other enviroment manager):
	 a. [Install 
	 virtualenv](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/#installing-virtualenv) 
	 if necessary:
	 ```
	 python3 -m pip install --user virtualenv
	 ```
	 b. Create a new environment
	 ```
	 cd chromosome-network-modeling
	 python3 -m venv env
	 ```
	 c. Activate and install requirements
	 ```
	 source env/bin/activate
	 pip install -r requirements.txt
	 ```
2. To reproduce the analysis, run the following script:
```
python3 run_test_distance_sample.py
```
