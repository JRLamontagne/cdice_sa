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
python config_reader_pol.py

# Construct input file for SALib sampling algorithm
echo "Constructing Input for SALib file..."
num_var=$(python construct_input_pol.py)

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
#python CDICE_input_maker_pol.py

#cp CDICE_input.txt export/CDICE_input.txt
#cd export

mpirun a.exe "CDICE_input.txt" 60000
#cp CDICE_output_* ../CDICE_output_*
#cp CDICE_output_0.txt ../CDICE_output_0.txt
#cp CDICE_output_1.txt ../CDICE_output_1.txt
#cp CDICE_output_2.txt ../CDICE_output_2.txt
#cp CDICE_output_3.txt ../CDICE_output_3.txt

# Parse CDICE_output.txt
echo $nProcs
mpirun python -c "from CDICE_output_parser2 import parsers; parsers(60)"
# python -c "from CDICE_output_parser import make_scalar; make_scalar()"

#Compute Sobol Indicies
mpirun python CDICE_batch_sobol_SCC_pol.py
mpirun python CDICE_batch_sobol_NPV_Dam_pol.py
mpirun python CDICE_batch_sobol_NPV_Abat_pol.py

python combining_npzs_pol.py

#Combine files into an output
cp -f CDICE_output_* Output_pol/
cp -f util.txt Output_pol/
cp -f pop.txt Output_pol/
cp -f util_disc.txt Output_pol/
cp -f real_intr.txt Output_pol/
cp -f NPV_Dam.txt Output_pol/
cp -f NPV_Abat.txt Output_pol/
cp -f Abat.txt Output_pol/
cp -f Dam.txt Output_pol/
cp -f GWP.txt Output_pol/
cp -f SCC.txt Output_pol/
cp -f Tatm.txt Output_pol/
cp -f Mat.txt Output_pol/
cp -f forc.txt Output_pol/
cp -f caa.txt Output_pol/
cp -f caa_up_tstep.txt Output_pol/
cp -f CDICE_input.txt Output_pol/
cp -f SCC_* Output_pol/
cp -f Dam_* Output_pol/
cp -f NPV_Dam_* Output_pol/
cp -f Abat_* Output_pol/
cp -f NPV_Abat_* Output_pol/

#Now only var
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

mpirun a.exe "CDICE_input.txt" 54000
#cp CDICE_output_* ../CDICE_output_*
#cp CDICE_output_0.txt ../CDICE_output_0.txt
#cp CDICE_output_1.txt ../CDICE_output_1.txt
#cp CDICE_output_2.txt ../CDICE_output_2.txt
#cp CDICE_output_3.txt ../CDICE_output_3.txt

# Parse CDICE_output.txt
echo $nProcs
mpirun python -c "from CDICE_output_parser2 import parsers; parsers(60)"
# python -c "from CDICE_output_parser import make_scalar; make_scalar()"

#Compute Sobol Indicies
mpirun python CDICE_batch_sobol_SCC.py
mpirun python CDICE_batch_sobol_NPV_Dam.py
mpirun python CDICE_batch_sobol_NPV_Abat.py

python combining_npzs.py

#Combine files into an output
cp -f CDICE_output_* Output_var/
cp -f util.txt Output_var/
cp -f pop.txt Output_var/
cp -f util_disc.txt Output_var/
cp -f real_intr.txt Output_var/
cp -f NPV_Dam.txt Output_var/
cp -f NPV_Abat.txt Output_var/
cp -f Abat.txt Output_var/
cp -f Dam.txt Output_var/
cp -f GWP.txt Output_var/
cp -f SCC.txt Output_var/
cp -f Tatm.txt Output_var/
cp -f Mat.txt Output_var/
cp -f forc.txt Output_var/
cp -f caa.txt Output_var/
cp -f caa_up_tstep.txt Output_var/
cp -f CDICE_input.txt Output_var/
cp -f SCC_* Output_var/
cp -f Dam_* Output_var/
cp -f NPV_Dam_* Output_var/
cp -f Abat_* Output_var/
cp -f NPV_Abat_* Output_var/

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
python config_reader_pol.py

# Construct input file for SALib sampling algorithm
echo "Constructing Input for SALib file..."
num_var=$(python construct_input_only_pol.py)

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
python CDICE_input_maker_only_pol.py

#cp CDICE_input.txt export/CDICE_input.txt
#cd export

mpirun a.exe "CDICE_input.txt" 8000
#cp CDICE_output_* ../CDICE_output_*
#cp CDICE_output_0.txt ../CDICE_output_0.txt
#cp CDICE_output_1.txt ../CDICE_output_1.txt
#cp CDICE_output_2.txt ../CDICE_output_2.txt
#cp CDICE_output_3.txt ../CDICE_output_3.txt

# Parse CDICE_output.txt
echo $nProcs
mpirun python -c "from CDICE_output_parser2 import parsers; parsers(60)"
# python -c "from CDICE_output_parser import make_scalar; make_scalar()"

#Compute Sobol Indicies
mpirun python CDICE_batch_sobol_SCC_only_pol.py
mpirun python CDICE_batch_sobol_NPV_Dam_only_pol.py
mpirun python CDICE_batch_sobol_NPV_Abat_only_pol.py

python combining_npzs_only_pol.py

#Combine files into an output
cp -f CDICE_output_* Output_only_pol/
cp -f util.txt Output_only_pol/
cp -f pop.txt Output_only_pol/
cp -f util_disc.txt Output_only_pol/
cp -f real_intr.txt Output_only_pol/
cp -f NPV_Dam.txt Output_only_pol/
cp -f NPV_Abat.txt Output_only_pol/
cp -f Abat.txt Output_only_pol/
cp -f Dam.txt Output_only_pol/
cp -f GWP.txt Output_only_pol/
cp -f SCC.txt Output_only_pol/
cp -f Tatm.txt Output_only_pol/
cp -f Mat.txt Output_only_pol/
cp -f forc.txt Output_only_pol/
cp -f caa.txt Output_only_pol/
cp -f caa_up_tstep.txt Output_only_pol/
cp -f CDICE_input.txt Output_only_pol/
cp -f SCC_* Output_only_pol/
cp -f Dam_* Output_only_pol/
cp -f NPV_Dam_* Output_only_pol/
cp -f Abat_* Output_only_pol/
cp -f NPV_Abat_* Output_only_pol/