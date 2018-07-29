#PBS -l nodes=4:ppn=15
#PBS -l walltime=40:00:00
#PBS -j oe
#PBS -o output.txt

cd $PBS_O_WORKDIR

module load python-2.7.5
#Ask user to input Note
#echo "Please input N:"
#read number
number=100
d=29
#echo "Please input nProcs:"
#read nProcs
nProcs=60
echo $nProcs
# Determine which factors are being tested, and their range from parameterization file
echo "Reading Configuration file..."
cd ..
python -c "from config_reader_function import config_reader; config_reader('configuration.txt','configuration.npz')"
python -c "from config_reader_function import config_reader; config_reader('configuration_pol.txt','configuration_pol.npz')"

# Construct input file for SALib sampling algorithm
echo "Constructing Input for SALib file..."
python -c "from construct_SALib_input import construct_input; construct_input(1)"

# Generate sample using Saltelli sampling, implemented in SALib
echo "Running SALib Sample..."
cd SALib-master/ # hack
python -m SALib.sample.saltelli \
     -n $number  \
     -p ./SALib_Input.txt \
     -o CDICE_SAvar.txt \
     --delimiter=' ' \
     --precision=8 \
     --max-order=2

# Options: N(2D + 2)
# Construct CDICE input sequences
cd ..
echo "Constructing CDICE Sequences..."
python  -c "from CDICE_input_maker_function import cdice_input_maker; cdice_input_maker(1)"

mpirun a.exe "CDICE_input.txt" $(expr $number \* $d \* 2 + $number \* 2)

# Parse CDICE_output.txt
echo $nProcs
mpirun python -c "from CDICE_output_parser2 import parsers; parsers(60)"