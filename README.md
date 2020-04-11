# Instructions

This code can be used to construct the datasets considered in the paper 

*Solving the Parametric Diffusion Equation by Deep Neural Networks - A Numerical Study*

It requires an installation of [FEniCS](https://fenicsproject.org/), which can, for example, easily be attained via `conda`by running
```conda install -c conda-forge fenics```.

## Recommended Workflow

1. Configure FEM parameters in `utility.py`.
2. Configure dataset parameters in `job_creator_data.sh`.
3. Run `job_creator_data.sh`.
4. Launch newly created jobs.

This results in the creation of a datasets folder, where each individual dataset is saved in a different subfolder, including the corresponding Gram matrix. Additionally, a serialized python dictionary documenting all parameters used in the data creation process is saved. If the jobs terminate successfully, the shell scripts are automatically removed.  

## Notes on Parallelization

This code allows for two ways to speed up the sample generation, which can be either used in conjunction or separately. 

Firstly, the dataset generation is by default broken up into parts configured via the `part_max` parameter. If, for example, `num_sample = 20000` and  `part_max = 5` ist chosen, then five jobs are created each yielding a dataset of size 4000. These are saved in the same folder and numbered accordingly. This way of parallelization is especially useful when running the code on any kind of computing cluster. 

Secondly, the `gen_data.py` script uses by default the python `multiprocessing` module. The default number of processes is set to two but can be adapted to the computing setup. If running this script on a single workstation, it might be useful to use less parts and instead  tailor the number of processes to the CPU.

## Notes on the Use of FEniCS Expressions

The `utility.py` script includes the generation of random samples for all our test-cases. The parameter functions $a_{y}$ are generated as FEniCS expressions. In general, this can be done via subclassing in python, C++ classes or simple C++ expressions. The first of these approaches is too slow for our use case. The other two should be equivalent in runtime.  We used the third approach which allows for compact writing, but sadly at the expense of readability. In particular, indicator functions are created by abusing, amongst other things, the ceil and floor functions on the unit square. Therefore, currently the Squares and Cookies datasets only work on the unit square. However, one could somewhat easily adapt those to other domains via scaling.   
