# CryoEVAL: Evaluation Metrics for Atomic Model Building Methods

This project is to evaluate the performance of different atomic model building approaches.


## Pre-requisites

| Dependency                        | Notes              |
|-----------------------------------|--------------------|
| Phenix                            | v1.21 |
| Python                            | 3.10               |
| Numpy                             | 1.26.0     |
| Scipy                             | 1.11.3     |
|torch                              |  2.1.1         |
|Bio                            | 1.8.1     |
|einops                         |  0.7.0 |
|tqdm                           |4.66.2                 |




## Installation

```bash
cd <your_project_path>

git clone https://github.com/chenwei-zhang/cryoEVAL.git
```



## Usage

```bash
cd cryoEVAL/src

python evaluate.py -p <predict_model.pdb> -t <target_model.pdb> -o <output_path>  (--modelangelo True --phenix True)
```

### Commands
- -p : Predicted model in PDB/CIF format
- -t : Reference model in PDB/CIF format
- -o : Output path

#### Optional
- --modelangelo True : Add ModelAngelo Paper's evaulation metrics
- --phenix True : Add phenix.chain_comparison evaluation metrics


### Example

```bash
cd cryoEVAL/src

python evaluate.py -p ../example/3j9s_agl.cif -t ../example/3j9s_ref.pdb -o ../example/eval_result.log
```
The evaluation results can be found in the `../example/eval_result.log` file.


<br>

## Credits

We would like to give credit to the following software projects that were partially reused in this project:

- [US-align](https://zhanggroup.org/US-align/)  
  License: Distributed under the Boost Software License, Version 1.0.

- [ModelAngelo](https://github.com/3dem/model-angelo)  
  License: Distributed under the MIT License.

- [Phenix](https://phenix-online.org/documentation/reference/chain_comparison.html)  
  License: Phenix is available without cost to users for non-profit work (typically academic and non-profit institutions).


## Support

Contact Chenwei Zhang at cwzhang@cs.ubc.ca


## Authors

Chenwei Zhang, CS PhD @UBC; Ex Machine Learning Intern @Amgen