# TopSol for CPVs  

Topology Optimization of metallization patterns in solar cells

## Introduction
**TopSol** is a MATLAB toolbox for full-scale modeling (using finite element method) and design optimization of metallization patterns in CPVs.

The code provides following features:
* Finite element modeling of rectangular solar cells
* Included resistance components are: front and back finger and busbar resistances, edge resistance, bulk resistance, shunt resistance and contact resistance.
* Allows adding variable temperature and illumination profiles
* Numerical model verified against previous benchmark results
* Large-scale optimization module included
* Several novel designs presented in the paper

Further, the code can be easily adapted for:
* Simple rectangular solar cells, see [here](https://link.springer.com/article/10.1007/s00158-014-1185-9)
* Freeform shapes, see [here](https://repository.tudelft.nl/islandora/object/uuid%3A6fba5499-0012-4708-a25c-c532da59e85c) and [here](https://repository.tudelft.nl/islandora/object/uuid%3Ad6d66ace-ec26-4971-9e3a-024c21a56015)
* Simultaneous optimization of front and rear metallizations, see [here](https://repository.tudelft.nl/islandora/object/uuid%3Ad6d66ace-ec26-4971-9e3a-024c21a56015)

## Authors
* Deepak K. Gupta ([guptadeepak2806[AT]@gmail.com](mailto:guptadeepak2806@gmail.com))*

## License
If you plan to distribute the software (commercially or not), please contact [Deepak Gupta](https://dkgupta90.github.io/topsol/) for more information.

## Dependencies
This framework has been tested on Matlab 2014b and 2016b.


## Usage
User documentation and test examples related to the paper will be added soon.

## Citation

If you use this code please use the following citation
```
@article{Gupta2018868,
title = "CPV solar cell modeling and metallization optimization",
journal = "Solar Energy",
volume = "159",
pages = "868 - 881",
year = "2018",
issn = "0038-092X",
doi = "https://doi.org/10.1016/j.solener.2017.11.015",
url = "http://www.sciencedirect.com/science/article/pii/S0038092X17309957",
author = "Deepak K. Gupta and Marco Barink and Matthijs Langelaar"
}
```
This paper can be freely downloaded after 29-11-2019 the embargo period [here](https://repository.tudelft.nl/islandora/object/uuid:a7b1277f-72f1-4069-a300-915a82d14281?collection=research). In the meantime, contact me, if you need a copy.

## Reporting Bugs
In case you experience any problems, please contact [Deepak Gupta](mailto:guptadeepak2806@gmail.com)
