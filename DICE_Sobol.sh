#!/bin/bash

#Ask user to input Note
echo "Please input N:"
read number
# Determine which factors are being tested, and their range from parameterization file
python config_reader.py

# Construct input file for SALib sampling algorithm
num_var=$(python construct_input.py)

# Generate sample using Saltelli sampling, implemented in SALib
cd SALib-master/ # hack
python -m SALib.sample.saltelli \
     -n $number \
     -p ./SALib_Input.txt \
     -o CDICE_SAvar.txt \
     --delimiter=' '  \
     --max-order=2
# Note: This assumes that parameters are from a uniform distribution
# Options:
# -p, --paramfile: Your parameter range file (3 columns: parameter name, lower bound, upper bound)
#
# -n, --samples: Sample size. 
#                    Number of model runs is N(2D + 2) if calculating second-order indices (default) 
#        or N(D + 2) otherwise.
#
# -o, --output: File to output your samples into.
# 
# --delimiter (optional): Output file delimiter.
#
# --precision (optional): Digits of precision in the output file. Default is 8.
#
# --max-order (optional): Maximum order of indices to calculate. Choose 1 or 2, default is 2. 
#                                          Choosing 1 will reduce total model runs from N(2D + 2) to N(D + 2)
#                                          Must use the same value (either 1 or 2) for both sampling and analysis.

# Construct CDICE input sequences
cd ..
python CDICE_input_maker.py

# Run the CDICE model
mpirun a.exe "CDICE_input.txt" 744

# Parse CDICE_output.txt
echo "Please input nProcs:"
read nProcs
python -c "from CDICE_output_parser import parsers; parsers($nProcs)"

# Now copy the content of the variable of interest to the SALib master directory
cp SCC.txt SALib-master/SALib_Output.txt

# Compute the sobol indicies
cd SALib-master/ # hack
python -m SALib.analyze.sobol \
     -p ./SALib_Input.txt \
     -Y ./SALib_Output.txt \
     -c 0 \
     --max-order=2 \
     -r 1000

# Options:
# -p, --paramfile: Your parameter range file (3 columns: parameter name, lower bound, upper bound)
#
# -Y, --model-output-file: File of model output values to analyze
#
# -c, --column (optional): Column of model output file to analyze. 
#                If the file only has one column, this argument will be ignored.
#
# --delimiter (optional): Model output file delimiter.
#
# --max-order (optional): Maximum order of indices to calculate.
#               This must match the value chosen during sampling.
# -r, --resamples (optional): Number of bootstrap resamples used to calculate confidence intervals on indices. Default 1000.
#