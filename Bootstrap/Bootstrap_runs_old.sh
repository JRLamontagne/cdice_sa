#PBS -l nodes=4:ppn=15
#PBS -l walltime=60:00:00
#PBS -j oe
#PBS -o output.txt

cd $PBS_O_WORKDIR

module load python-2.7.5

#read nProcs
nProcs=60
echo $nProcs
d=29



# Determine which factors are being tested, and their range from parameterization file
echo "Reading Configuration file..."
cd ..
python -c "from config_reader_function import config_reader; config_reader('configuration.txt','configuration.npz')"
python -c "from config_reader_function import config_reader; config_reader('configuration_pol.txt','configuration_pol.npz')"

# Construct input file for SALib sampling algorithm
echo "Constructing Input for SALib file..."
python -c "from construct_SALib_input import construct_input; construct_input(1)"

#######################################################################################################
#read number
number=100000

# Generate sample using Saltelli sampling, implemented in SALib
echo "Running SALib Sample..."
python -c "from sobol_bootstrap import sobol_sample; sobol_sample('SALib_Input.txt','CDICE_SAvar.txt',100000)"

# Construct CDICE input sequences
echo "Constructing CDICE Sequences..."
python  -c "from CDICE_input_maker_function import cdice_input_maker; cdice_input_maker(1)"

mpirun a.exe "CDICE_input.txt" $(expr $number \* $d \* 2 + $number \* 2)

# Parse CDICE_output.txt
echo $nProcs
./parsing_CDICE.sh

# Compute Sobol Indicies
python -c "from sobol_bootstrap import sobol_index; sobol_index('SALib_Input.txt',100000,100)"
#######################################################################################################