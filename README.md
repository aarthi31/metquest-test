# MetQuest

MetQuest is a dynamic programming based algorithm for identifying all possible pathways from metabolic networks between the source and the target metabolites.

## Getting started

1. Clone this repository to your computer using ```git``` or [download the repository](https://github.com/aarthi31/MetQuest/) and decompress it. All the executable codes can be found in the folder Codes.  
2. Install [Python 2.7.14](https://www.python.org/downloads/)

Ensure that the following packages (with the specified version) are installed:
[NetworkX 1.9.1](http://networkx.github.io/), [cobrapy 0.5.4](https://github.com/opencobra/cobrapy) and
[libSBML 5.13.0](http://sbml.org/Software/libSBML/docs/python-api/libsbml-downloading.html)


## Input

Folder whose structure is as shown:
```
   
    ├── Example                 # Folder  
    │   ├── SBML model(s) of metabolic networks          # XML files of the metabolic networks (COBRA-compatible)
    │   ├── seed_mets.txt         # Text file containing the seed metabolites separated by a newline
    │   ├── source_mets.txt       # Text file containing the source metabolites separated by a newline
    │   └── target_mets.txt       # Text file containing the target metabolites separated by a newline
    └── ...
 ```

## Running MetQuest


## Authors

* [Aarthi Ravikrishnan](https://github.com/aarthi31)
* Meghana Nasre
* [Karthik Raman](https://github.com/karthikraman)


## License

1. By using the software enclosed in this package (MetQuest), you agree to become bound by the terms of this license. 
2. This software is for your internal use only. Please DO NOT redistribute it without the permission from the authors.
3. This software is for **academic use only**. No other usage is allowed without a written permission from the authors. It cannot be used for any commercial interest.
4. The authors appreciate it if you can send us your feedback including any bug report.
5. The authors do not hold any responsibility for the correctness of this software, though we cross-checked all experimental results.

## Acknowledgments

This work was supported by the Indian Institute of Technology Madras grant ERP/1314/004/RESF/KARH to KR and the INSPIRE fellowship, Department of Science and Technology, Government of India to AR.


