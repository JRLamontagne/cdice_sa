#PBS -l nodes=4:ppn=15
#PBS -l walltime=40:00:00
#PBS -j oe
#PBS -o output.txt

cd $PBS_O_WORKDIR

module load python-2.7.5
#Ask user to input Note
#echo "Please input N:"
#read number
number=1000
#echo "Please input nProcs:"
#read nProcs
nProcs=60
echo $nProcs
# Determine which factors are being tested, and their range from parameterization file
echo "Reading Configuration file..."
python config_reader.py

# Construct input file for SALib sampling algorithm
echo "Constructing Input for SALib file..."
num_var=$(python construct_input.py)

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
python CDICE_input_maker.py

#cp CDICE_input.txt export/CDICE_input.txt
#cd export

mpirun a.exe "CDICE_input.txt" 6200000
#cp CDICE_output_* ../CDICE_output_*
#cp CDICE_output_0.txt ../CDICE_output_0.txt
#cp CDICE_output_1.txt ../CDICE_output_1.txt
#cp CDICE_output_2.txt ../CDICE_output_2.txt
#cp CDICE_output_3.txt ../CDICE_output_3.txt

# Parse CDICE_output.txt
echo $nProcs
mpirun python -c "from CDICE_output_parser2 import parsers; parsers(60)"
python -c "from CDICE_output_parser import make_scalar; make_scalar()"

#Compute Sobol Indicies
mpirun python CDICE_batch_sobol_SCC.py
mpirun python CDICE_batch_sobol_NPV_Dam.py
mpirun python CDICE_batch_sobol_NPV_Abat.py

python combining_npzs.py