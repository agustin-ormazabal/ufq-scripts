export LC_NUMERIC="en_US.UTF-8"
#Parms
rst_con=300.0
kon=$(echo 2.0*$rst_con | bc)
ini=12.020 
inc=0.1
#fin=26.9 se despega
fin=30.0


rm func.* temp*out delta.*
rm entropias*dat
rm -r input
mkdir input
mkdir error_analisis

rm variables
echo "#si o no al wham" >variables
echo ".FALSE." >>variables
echo "#wham iter"  >>variables
echo "1000000"  >>variables
echo "#si o no al dham" >>variables 
echo ".FALSE." >>variables 
echo "#si o no al bar" >>variables 
echo ".TRUE." >>variables 
echo "#bariter" >>variables 
echo "400" >>variables 
echo "#si o no a los errores" >>variables 
echo ".TRUE." >>variables 
echo "#es angle restraint" >>variables 
echo ".FALSE." >>variables 

rm wham.inp wham.out bar.out 
for x in `seq $ini $inc $fin`; do
rm inp.$x
head -n 990000 ../val/v.$x | tail -n 800000 > inp.$x
echo "inp.$x $x $kon" >> wham.inp
done
f95 wham.gusp.f  -o  wham_mio.exe -l blas -l lapack  -mcmodel=large
./wham_mio.exe variables

for x in `seq $ini $inc $fin`; do
mv inp.$x input/
done

mv entropias.* error_analisis/
mv func.* error_analisis/
mv wham.inp input/
