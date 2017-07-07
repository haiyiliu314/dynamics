#07/03/2017 creation
#copy codes to a new file, run them, and create log file containing date and run time
mkdir run5
cd run5
cp ~/dynamics/Formal_V1_0703/constants.f90 ~/dynamics/Formal_V1_0703/run5
cp ~/dynamics/Formal_V1_0703/dynamics_main.f90 ~/dynamics/Formal_V1_0703/run5
cp ~/dynamics/Formal_V1_0703/Makefile ~/dynamics/Formal_V1_0703/run5
cp ~/dynamics/Formal_V1_0703/N_of_grids.f90 ~/dynamics/Formal_V1_0703/run5
cp ~/dynamics/Formal_V1_0703/RK_help.f90 ~/dynamics/Formal_V1_0703/run5
cp ~/dynamics/Formal_V1_0703/RK_module.f90 ~/dynamics/Formal_V1_0703/run5
cp ~/dynamics/Efreq_Ben.dat ~/dynamics/Formal_V1_0703/run5
cp ~/dynamics/Formal_V1_0703/params.f90 ~/dynamics/Formal_V1_0703/run5
cp ~/dynamics/Formal_V1_0703/coul_mat_module.f90 ~/dynamics/Formal_V1_0703/run5
cp ~/dynamics/0614_more_decay/getfigure_06062017.m ~/dynamics/Formal_V1_0703/run5
cp ~/dynamics/Formal_V1_0703/SBE.sh ~/dynamics/Formal_V1_0703/run5
make
{ time ./run.x inputfile* > log.dat; } 2> timing.txt
echo 1>> log.dat
echo "Timing information" 1>> log.dat
NOW=$( date +"%m-%d-%y-%T")
echo $NOW 1>> log.dat
cat timing.txt >> log.dat
make clean
