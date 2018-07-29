#PBS -l nodes=4:ppn=16
#PBS -l walltime=40:00:00
#PBS -j oe
#PBS -o output.txt

cd $PBS_O_WORKDIR

module load python-2.7.5

#read nProcs
nProcs=64
echo $nProcs
d=27



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
number=1000

# Generate sample using Saltelli sampling, implemented in SALib
echo "Running SALib Sample..."
python -c "from sobol_bootstrap import sobol_sample; sobol_sample('SALib_Input.txt','CDICE_SAvar.txt',1000)"

# Construct CDICE input sequences
echo "Constructing CDICE Sequences..."
python  -c "from CDICE_input_maker_function import cdice_input_maker; cdice_input_maker(1,'exp',2050)"

mpirun a.exe "CDICE_input.txt" $(expr $number \* $d \* 2 + $number \* 2)

# Parse CDICE_output.txt
echo $nProcs
./parsing_CDICE.sh

# Compute Sobol Indicies
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Dam.txt',2,1000,27,'Dam')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Abat.txt',2,1000,27,'Abat')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Tatm.txt',2,1000,27,'Tatm')"

python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Dam.txt',100,17,17,1000,27,'NPV_Dam')"
python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Abat.txt',100,17,17,1000,27,'NPV_Abat')"

mv CDICE_output_*.txt Exp/2050/
mv CDICE_input.txt Exp/2050/
mv NPV_*_17.npz Exp/2050/
mv ygross.txt Exp/2050/
mv caa*.txt Exp/2050/
mv forc.txt Exp/2050/
mv Mat.txt Exp/2050/
mv Tatm.txt Exp/2050/
mv SCC.txt Exp/2050/
mv GWP.txt Exp/2050/
mv Dam.txt Exp/2050/
mv Abat.txt Exp/2050/
mv NPV_*.txt Exp/2050/
mv real_intr.txt Exp/2050/
mv util_disc.txt Exp/2050/
mv pop.txt Exp/2050/
mv util.txt Exp/2050/
mv SALib-master/CDICE_SAvar.txt Exp/2050/

#######################################################################################################
#read number
number=1000

# Generate sample using Saltelli sampling, implemented in SALib
echo "Running SALib Sample..."
python -c "from sobol_bootstrap import sobol_sample; sobol_sample('SALib_Input.txt','CDICE_SAvar.txt',1000)"

# Construct CDICE input sequences
echo "Constructing CDICE Sequences..."
python  -c "from CDICE_input_maker_function import cdice_input_maker; cdice_input_maker(1,'exp',2030)"

mpirun a.exe "CDICE_input.txt" $(expr $number \* $d \* 2 + $number \* 2)

# Parse CDICE_output.txt
echo $nProcs
./parsing_CDICE.sh

# Compute Sobol Indicies
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Dam.txt',2,1000,27,'Dam')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Abat.txt',2,1000,27,'Abat')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Tatm.txt',2,1000,27,'Tatm')"

python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Dam.txt',100,17,17,1000,27,'NPV_Dam')"
python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Abat.txt',100,17,17,1000,27,'NPV_Abat')"

mv CDICE_output_*.txt Exp/2030/
mv CDICE_input.txt Exp/2030/
mv NPV_*_17.npz Exp/2030/
mv ygross.txt Exp/2030/
mv caa*.txt Exp/2030/
mv forc.txt Exp/2030/
mv Mat.txt Exp/2030/
mv Tatm.txt Exp/2030/
mv SCC.txt Exp/2030/
mv GWP.txt Exp/2030/
mv Dam.txt Exp/2030/
mv Abat.txt Exp/2030/
mv NPV_*.txt Exp/2030/
mv real_intr.txt Exp/2030/
mv util_disc.txt Exp/2030/
mv pop.txt Exp/2030/
mv util.txt Exp/2030/
mv SALib-master/CDICE_SAvar.txt Exp/2030/

#######################################################################################################
#read number
number=1000

# Generate sample using Saltelli sampling, implemented in SALib
echo "Running SALib Sample..."
python -c "from sobol_bootstrap import sobol_sample; sobol_sample('SALib_Input.txt','CDICE_SAvar.txt',1000)"

# Construct CDICE input sequences
echo "Constructing CDICE Sequences..."
python  -c "from CDICE_input_maker_function import cdice_input_maker; cdice_input_maker(1,'exp',2040)"

mpirun a.exe "CDICE_input.txt" $(expr $number \* $d \* 2 + $number \* 2)

# Parse CDICE_output.txt
echo $nProcs
./parsing_CDICE.sh

# Compute Sobol Indicies
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Dam.txt',2,1000,27,'Dam')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Abat.txt',2,1000,27,'Abat')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Tatm.txt',2,1000,27,'Tatm')"

python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Dam.txt',100,17,17,1000,27,'NPV_Dam')"
python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Abat.txt',100,17,17,1000,27,'NPV_Abat')"

mv CDICE_output_*.txt Exp/2040/
mv CDICE_input.txt Exp/2040/
mv NPV_*_17.npz Exp/2040/
mv ygross.txt Exp/2040/
mv caa*.txt Exp/2040/
mv forc.txt Exp/2040/
mv Mat.txt Exp/2040/
mv Tatm.txt Exp/2040/
mv SCC.txt Exp/2040/
mv GWP.txt Exp/2040/
mv Dam.txt Exp/2040/
mv Abat.txt Exp/2040/
mv NPV_*.txt Exp/2040/
mv real_intr.txt Exp/2040/
mv util_disc.txt Exp/2040/
mv pop.txt Exp/2040/
mv util.txt Exp/2040/
mv SALib-master/CDICE_SAvar.txt Exp/2040/

#######################################################################################################
#read number
number=1000

# Generate sample using Saltelli sampling, implemented in SALib
echo "Running SALib Sample..."
python -c "from sobol_bootstrap import sobol_sample; sobol_sample('SALib_Input.txt','CDICE_SAvar.txt',1000)"

# Construct CDICE input sequences
echo "Constructing CDICE Sequences..."
python  -c "from CDICE_input_maker_function import cdice_input_maker; cdice_input_maker(1,'linear',2030)"

mpirun a.exe "CDICE_input.txt" $(expr $number \* $d \* 2 + $number \* 2)

# Parse CDICE_output.txt
echo $nProcs
./parsing_CDICE.sh

# Compute Sobol Indicies
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Dam.txt',2,1000,27,'Dam')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Abat.txt',2,1000,27,'Abat')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Tatm.txt',2,1000,27,'Tatm')"

python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Dam.txt',100,17,17,1000,27,'NPV_Dam')"
python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Abat.txt',100,17,17,1000,27,'NPV_Abat')"

mv CDICE_output_*.txt Linear/2030/
mv CDICE_input.txt Linear/2030/
mv NPV_*_17.npz Linear/2030/
mv ygross.txt Linear/2030/
mv caa*.txt Linear/2030/
mv forc.txt Linear/2030/
mv Mat.txt Linear/2030/
mv Tatm.txt Linear/2030/
mv SCC.txt Linear/2030/
mv GWP.txt Linear/2030/
mv Dam.txt Linear/2030/
mv Abat.txt Linear/2030/
mv NPV_*.txt Linear/2030/
mv real_intr.txt Linear/2030/
mv util_disc.txt Linear/2030/
mv pop.txt Linear/2030/
mv util.txt Linear/2030/
mv SALib-master/CDICE_SAvar.txt Linear/2030/

#######################################################################################################
#read number
number=1000

# Generate sample using Saltelli sampling, implemented in SALib
echo "Running SALib Sample..."
python -c "from sobol_bootstrap import sobol_sample; sobol_sample('SALib_Input.txt','CDICE_SAvar.txt',1000)"

# Construct CDICE input sequences
echo "Constructing CDICE Sequences..."
python  -c "from CDICE_input_maker_function import cdice_input_maker; cdice_input_maker(1,'linear',2040)"

mpirun a.exe "CDICE_input.txt" $(expr $number \* $d \* 2 + $number \* 2)

# Parse CDICE_output.txt
echo $nProcs
./parsing_CDICE.sh

# Compute Sobol Indicies
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Dam.txt',2,1000,27,'Dam')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Abat.txt',2,1000,27,'Abat')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Tatm.txt',2,1000,27,'Tatm')"

python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Dam.txt',100,17,17,1000,27,'NPV_Dam')"
python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Abat.txt',100,17,17,1000,27,'NPV_Abat')"

mv CDICE_output_*.txt Linear/2040/
mv CDICE_input.txt Linear/2040/
mv NPV_*_17.npz Linear/2040/
mv ygross.txt Linear/2040/
mv caa*.txt Linear/2040/
mv forc.txt Linear/2040/
mv Mat.txt Linear/2040/
mv Tatm.txt Linear/2040/
mv SCC.txt Linear/2040/
mv GWP.txt Linear/2040/
mv Dam.txt Linear/2040/
mv Abat.txt Linear/2040/
mv NPV_*.txt Linear/2040/
mv real_intr.txt Linear/2040/
mv util_disc.txt Linear/2040/
mv pop.txt Linear/2040/
mv util.txt Linear/2040/
mv SALib-master/CDICE_SAvar.txt Linear/2040/

#######################################################################################################
#read number
number=1000

# Generate sample using Saltelli sampling, implemented in SALib
echo "Running SALib Sample..."
python -c "from sobol_bootstrap import sobol_sample; sobol_sample('SALib_Input.txt','CDICE_SAvar.txt',1000)"

# Construct CDICE input sequences
echo "Constructing CDICE Sequences..."
python  -c "from CDICE_input_maker_function import cdice_input_maker; cdice_input_maker(1,'quad',2050)"

mpirun a.exe "CDICE_input.txt" $(expr $number \* $d \* 2 + $number \* 2)

# Parse CDICE_output.txt
echo $nProcs
./parsing_CDICE.sh

# Compute Sobol Indicies
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Dam.txt',2,1000,27,'Dam')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Abat.txt',2,1000,27,'Abat')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Tatm.txt',2,1000,27,'Tatm')"

python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Dam.txt',100,17,17,1000,27,'NPV_Dam')"
python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Abat.txt',100,17,17,1000,27,'NPV_Abat')"

mv CDICE_output_*.txt quad/2050/
mv CDICE_input.txt quad/2050/
mv NPV_*_17.npz quad/2050/
mv ygross.txt quad/2050/
mv caa*.txt quad/2050/
mv forc.txt quad/2050/
mv Mat.txt quad/2050/
mv Tatm.txt quad/2050/
mv SCC.txt quad/2050/
mv GWP.txt quad/2050/
mv Dam.txt quad/2050/
mv Abat.txt quad/2050/
mv NPV_*.txt quad/2050/
mv real_intr.txt quad/2050/
mv util_disc.txt quad/2050/
mv pop.txt quad/2050/
mv util.txt quad/2050/
mv SALib-master/CDICE_SAvar.txt quad/2050/
###python -c "from sobol_bootstrap import sobol_index; sobol_index('SALib_Input.txt',1000,100)"
#python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Dam.txt',100,17,17,1000,28)"
#mv Boot_100000.npz NPV_Dam.npz
#python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Abat.txt',100,17,17,1000,28)"
#mv Boot_100000.npz NPV_Abat.npz
#python -c "from custom_sobol import custom_sobol; custom_sobol('Dam.txt',100,0,17,1000,28)"
#mv Boot_100000.npz Dam.npz
#python -c "from custom_sobol import custom_sobol; custom_sobol('Abat.txt',100,0,17,1000,28)"
#mv Boot_100000.npz Abat.npz
#python -c "from custom_sobol import custom_sobol; custom_sobol('Tatm.txt',100,0,17,1000,28)"
#mv Boot_100000.npz Tatm.npz
#######################################################################################################
#read number
number=1000

# Generate sample using Saltelli sampling, implemented in SALib
echo "Running SALib Sample..."
python -c "from sobol_bootstrap import sobol_sample; sobol_sample('SALib_Input.txt','CDICE_SAvar.txt',1000)"

# Construct CDICE input sequences
echo "Constructing CDICE Sequences..."
python  -c "from CDICE_input_maker_function import cdice_input_maker; cdice_input_maker(1,'quad',2040)"

mpirun a.exe "CDICE_input.txt" $(expr $number \* $d \* 2 + $number \* 2)

# Parse CDICE_output.txt
echo $nProcs
./parsing_CDICE.sh

# Compute Sobol Indicies
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Dam.txt',2,1000,27,'Dam')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Abat.txt',2,1000,27,'Abat')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Tatm.txt',2,1000,27,'Tatm')"

python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Dam.txt',100,17,17,1000,27,'NPV_Dam')"
python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Abat.txt',100,17,17,1000,27,'NPV_Abat')"

mv CDICE_output_*.txt quad/2040/
mv CDICE_input.txt quad/2040/
mv NPV_*_17.npz quad/2040/
mv ygross.txt quad/2040/
mv caa*.txt quad/2040/
mv forc.txt quad/2040/
mv Mat.txt quad/2040/
mv Tatm.txt quad/2040/
mv SCC.txt quad/2040/
mv GWP.txt quad/2040/
mv Dam.txt quad/2040/
mv Abat.txt quad/2040/
mv NPV_*.txt quad/2040/
mv real_intr.txt quad/2040/
mv util_disc.txt quad/2040/
mv pop.txt quad/2040/
mv util.txt quad/2040/
mv SALib-master/CDICE_SAvar.txt quad/2040/
###python -c "from sobol_bootstrap import sobol_index; sobol_index('SALib_Input.txt',1000,100)"
#python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Dam.txt',100,17,17,1000,28)"
#mv Boot_100000.npz NPV_Dam.npz
#python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Abat.txt',100,17,17,1000,28)"
#mv Boot_100000.npz NPV_Abat.npz
#python -c "from custom_sobol import custom_sobol; custom_sobol('Dam.txt',100,0,17,1000,28)"
#mv Boot_100000.npz Dam.npz
#python -c "from custom_sobol import custom_sobol; custom_sobol('Abat.txt',100,0,17,1000,28)"
#mv Boot_100000.npz Abat.npz
#python -c "from custom_sobol import custom_sobol; custom_sobol('Tatm.txt',100,0,17,1000,28)"
#mv Boot_100000.npz Tatm.npz
###########################################################################
#read number
number=1000

# Generate sample using Saltelli sampling, implemented in SALib
echo "Running SALib Sample..."
python -c "from sobol_bootstrap import sobol_sample; sobol_sample('SALib_Input.txt','CDICE_SAvar.txt',1000)"

# Construct CDICE input sequences
echo "Constructing CDICE Sequences..."
python  -c "from CDICE_input_maker_function import cdice_input_maker; cdice_input_maker(1,'quad',2030)"

mpirun a.exe "CDICE_input.txt" $(expr $number \* $d \* 2 + $number \* 2)

# Parse CDICE_output.txt
echo $nProcs
./parsing_CDICE.sh

# Compute Sobol Indicies
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Dam.txt',2,1000,27,'Dam')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Abat.txt',2,1000,27,'Abat')"
mpirun python -c "from parallel_sobol import parallel_sobol; parallel_sobol('Tatm.txt',2,1000,27,'Tatm')"

python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Dam.txt',100,17,17,1000,27,'NPV_Dam')"
python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Abat.txt',100,17,17,1000,27,'NPV_Abat')"

mv CDICE_output_*.txt quad/2030/
mv CDICE_input.txt quad/2030/
mv NPV_*_17.npz quad/2030/
mv ygross.txt quad/2030/
mv caa*.txt quad/2030/
mv forc.txt quad/2030/
mv Mat.txt quad/2030/
mv Tatm.txt quad/2030/
mv SCC.txt quad/2030/
mv GWP.txt quad/2030/
mv Dam.txt quad/2030/
mv Abat.txt quad/2030/
mv NPV_*.txt quad/2030/
mv real_intr.txt quad/2030/
mv util_disc.txt quad/2030/
mv pop.txt quad/2030/
mv util.txt quad/2030/
mv SALib-master/CDICE_SAvar.txt quad/2030/
###python -c "from sobol_bootstrap import sobol_index; sobol_index('SALib_Input.txt',1000,100)"
#python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Dam.txt',100,17,17,1000,28)"
#mv Boot_100000.npz NPV_Dam.npz
#python -c "from custom_sobol import custom_sobol; custom_sobol('NPV_Abat.txt',100,17,17,1000,28)"
#mv Boot_100000.npz NPV_Abat.npz
#python -c "from custom_sobol import custom_sobol; custom_sobol('Dam.txt',100,0,17,1000,28)"
#mv Boot_100000.npz Dam.npz
#python -c "from custom_sobol import custom_sobol; custom_sobol('Abat.txt',100,0,17,1000,28)"
#mv Boot_100000.npz Abat.npz
#python -c "from custom_sobol import custom_sobol; custom_sobol('Tatm.txt',100,0,17,1000,28)"
#mv Boot_100000.npz Tatm.npz