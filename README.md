# Reaction-based library enumeration
Reaction-based library enumeration of small molecules from building blocks

This Code is written by:
## Bikash Debnath

## Purpose:
This code will allow to generate library of small molecules from building blocks based on reactions provided by the user.

## Requirements 

1. Python 3.x
2. RDKit

create a separate virtual environment (my-rdkit-env) and install rdkit



## Tutorial
Step 1: Activate conda environment

$conda activate my-rdkit-env

Step 2: Run Below code to generate library.

$python lib_enumeration.py <BB_1.smi> <BB_2.smi> <Products.smi>

For example:
I will generate thioether from two building blocks (BB-thiol.smi -- a library of thiols; BB-acrylamides.smi -- a library of acrylamides) from a reaction (reaction smarts are provided in the "lib_enumeration.py" script. For new excercise, uses needs to provide custom building blocks and reaction smarts.

$python lib_enumeration.py BB_thiol.smi BB_acrylamides.smi.smi test_prod.smi


## Contributing
Contributions are welcome! If you have ideas for improvements or new features, please open an issue or submit a pull request.

## License
This project is licensed under the  GNU General Public License v3.0 - see the LICENSE file for details.

## Contact
If you have any question please contact: bikashd2000@gmail.com
