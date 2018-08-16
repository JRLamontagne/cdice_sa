#PBS -l nodes=8:ppn=4
#PBS -l walltime=40:00:00
#PBS -j oe
#PBS -o output.txt

cd $PBS_O_WORKDIR

module load python-2.7.5

#read nProcs
nProcs=32
echo $nProcs
d=25



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
python  -c "from CDICE_input_maker_function import cdice_input_maker; cdice_input_maker(1,'exp2',2050)"

mpirun a.exe "CDICE_input.txt" $(expr $number \* $d \* 2 + $number \* 2)

# Parse CDICE_output.txt
echo $nProcs
./parsing_CDICE.sh

# Compute Sobol Indicies
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Dam.txt',100,100000,25,'Dam')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Abat.txt',100,100000,25,'Abat')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Tatm.txt',100,100000,25,'Tatm')"

python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Dam.txt',100,27,27,100000,25,'NPV_Dam')"
python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Abat.txt',100,27,27,100000,25,'NPV_Abat')"

python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Dam.txt',100,17,17,100000,25,'NPV_Dam2')"
python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Abat.txt',100,17,17,100000,25,'NPV_Abat2')"