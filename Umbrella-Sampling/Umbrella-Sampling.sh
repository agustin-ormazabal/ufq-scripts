#  Este programa ejecuta Umbrella Sampling para simular el proceso de unbinding
# entre dos moléculas. Sus inputs son la topología (top) y coordenadas (crd)
# del sistema en formato AMBER. Es necesario definir la coordenada de reacción
# utiliando los números de los átomos seleccionados, y el número de residuo
# dentro del cual se encuentran.

                                                                                      #############################              .          
top="/mnt/largedisk/newhome/gusp/Agus/L/Modelo_Box2/Tleap/Un_L_Box2.prmtop"           # <- Toplogia               # .   ##############     .     .
crd="/mnt/largedisk/newhome/gusp/Agus/L/Modelo_Box2/Equilibracion/eq_Un_L_1Prot.rst"  # <- Coordenadas            #     #  _     _   #   .   .     .
                                                                                      #############################    .# |_|   |_|  #     .    .   
       	                                                                                                                #  _     _   #      .  .      
                                                     ########################################################       .   # |_|   |_|  # #############  .
prot=":119,130,141,152,163,178,189,200,211,222@CA"   # <-  Con esto definis los residuos de la proteina     #        .  #  _     _   # #  _     _  # . 
                                                     #    que vas a usar para Umbrella                      #      .    # |_|   |_|  # # |_|   |_| #  .
sls=":238-252@C1',:295-307@C1'"                      # <-  Con esto definis los residuos del RNA que vas    #           #  _     _   # #  _     _  # . 
                                                     #    a usar para Umbrella                              #         . # |_|   |_|  # # |_|   |_| #
                                                     ########################################################           #  _     _   # #  _     _  #  .
                                                                                                                     	# |_|   |_|  # # |_|   |_| # 
########################################################                                                             	#     _      # #     _     # .
#  Con esto definis los atomos de los carbonos alfa    #                                                         .   	#    | |     # #    | |    #
# de las proteinas que vas a usar para el umbrella  -> #                                                             	#    | |     # #    | |    #  .
########################################################                                                            ###################################
                                                                                                                   #####################################

igr1="1847,2030,2198,2358,2532,2768,2951,3119,3279,3453"
                                                                                        
########################################################                                
#   Con esto definis los atomos de los C1' del RNA     #                                
#          que vas a usar para el umbrella          -> #
######################################################## 
                                 
igr2="3723,3757,3787,3818,3852,3885,3916,3950,3984,4017,4047,4080,4114,4147,4178,5572,5603,5636,5666,5697,5730,5764,5798,5831,5862,5896,5929,5959"

                #######################################################
                #      Define la coordenada inicial del Umbrella      #
                #######################################################


#  Primero calcula el centro de masa entre los átomos de SL2 y SL3, y lo mismo
# hace con la proteina anclada entre ellos.

rm cm.in

cat >> cm.in << EOF                      
trajin $crd lastframe
vector prot center out prot_tmp $prot
go
vector sls  center out sls_tmp  $sls
go
quit
EOF

cpptraj $top cm.in 
                                                           #####################
grep -v "prot" prot_tmp |awk {print'$2,$3,$4'} > prot_cm   # Se queda sólo con #
grep -v "sls"  sls_tmp  |awk {print'$2,$3,$4'} > sls_cm    #  las coordenadas  #
                                                           #####################
rm *tmp

#  Luego calcula la distacia entre los centros de masa.

xp=$(cat prot_cm |awk '{print $1}')    ##############################################                  
xs=$(cat sls_cm  |awk '{print $1}')    #                                            #
                                       #                                            #
yp=$(cat prot_cm |awk '{print $2}')    #     Aca definis las coordenadas en         #
ys=$(cat sls_cm  |awk '{print $2}')    #     X, Y y Z para poder operar sobre       #
                                       #     ellas directamente.                    #
zp=$(cat prot_cm |awk '{print $3}')    #                                            #
zs=$(cat sls_cm  |awk '{print $3}')    ##############################################


xf=$(bc <<< "$xs-$xp")  #################################################################
yf=$(bc <<< "$ys-$yp")  #  Aca definis el vector que pasa por los dos centros de masa   #
zf=$(bc <<< "$zs-$zp")  #################################################################


                                                       ###################################
cuadrados=$(bc <<< "$xf*$xf + $yf*$yf + $zf*$zf")      #  Calcula el modulo del vector,  #
                                                       # el cual es igual a la distancia # 
distancia=$(bc <<< "$(echo "sqrt($cuadrados)" | bc)")  #  a la distancia entre los cm's  #
                                                       ###################################

echo $distancia



                      #######################################################
                      #                 Umbrella Sampling                   #
                      #######################################################


# Se define la coordenada final del Umbrella y el crecimiento en cada paso.

inicio=$distancia
crecimiento=0.08
final=50.9158


# Se definen las condiciones para la dinamica en cada paso de la reaccion.

                  ###############################################################
nstlimp=500000    # <- Numero de pasos para la Produccion (femptosegundos)      #
nstlime=500000    # <- Numero de pasos para la Equilibrizacion (femptosegundos) #
idmumbt=2         # <- Tamaño del paso de integracion para la Produccion        #
ku=300.0          # <- Constante para la restriccion                            #
                  ###############################################################


# Define todos los directorios de trabajo para la Equilibrizacion y la Produccion

rm -r eq prod val


mkdir eq              #############################################
mkdir prod            # Crea las carpetas a donde van los outputs #
mkdir val             #############################################



Eqfiles="eq"      ######################################################
Prodfiles="prod"  # Indica que esas carpetas es a donde van los inputs #
Values="val"      ######################################################



# A partir de aca empieza el loopeo para ejecutar la dinamica en cada paso

cp $crd inp.crd 

for x in $(seq $inicio $crecimiento $final | sed "s/,/\./g"); do

# Genera los inputs para la Equilibrizacion en la carpeta "eq"

cat > $Eqfiles/mdin_eq.$x.inp <<EOF
 umbrella sampling
 &cntrl
  imin=0,irest=1,ioutfm=1,nmropt=1,
  ntx=5,ntb=2,ntp=1,ntf=2,ntt=1,ntc=2, 
  tempi=300.0, temp0=300.0,pres0=1.0,
  cut=10.0,nstlim=$nstlime,ntr=1,
  ntwx=1000,ntpr=1000,ntwr=10000,
 &end
 &wt
 type='DUMPFREQ', istep1=1,
 &end
 &wt
 type="END",
 &end
 DISANG=$Eqfiles/d.$x.is.RST
 DUMPAVE=$Values/v.$x.is.va
 &end
keep rna fixed with weak restraints
3.0
FIND
C3' * * *
C5' * * *
SEARCH
RES 308 308
END
END
EOF

# Ahora se genera el archivo con las restricciones


cat > $Eqfiles/d.$x.is.RST <<EOF
 &rst
 iat= -1,-1, r1=-30, r2=$x, r3=$x, r4=300, rk2=$ku, rk3=$ku,iresid=0,
 igr1=$igr1
 igr2=$igr2
 &end
EOF


# Crea el archivo de ejecucion y larga la produccion !!!

export CUDA_VISIBLE_DEVICES="1"
taskset -c 1 pmemd.cuda -O -i $Eqfiles/mdin_eq.$x.inp -o $Eqfiles/mdin_eq.$x.out -p $top -c inp.crd -r $Eqfiles/mdin_eq.$x.crd -inf mdinfo -x $Eqfiles/eq.$x.traj -ref inp.crd

                                    #######################################
                                    #      C      O      D        A       #      
                                    #######################################

# Esto genera los inputs para la Produccion, que se llaman "mdin_prod.$x.$cord_b_inicial.inp" y van a estar en el directorio "prod"

sed -e "s/nstlim=$nstlime/nstlim=$nstlimp/" $Eqfiles/mdin_eq.$x.inp > temp
sed -e "s/istep1=2/istep1=$idmumbt/" temp > temp1
sed -e "s/is.va/pro/" temp1 > temp2
sed -e "s/ntwx=10/ntwx=400/" temp2 > temp4
sed -e "s/ntpr=10/ntpr=400/" temp4 > temp5
sed -e "s/ntwr=10/ntwr=$nstlimp/" temp5 > $Prodfiles/mdin_prod.$x.inp
rm temp temp1 temp2 temp4 temp5


# Copia el último .rst de la Equilibrizacion con el prefijo "inp", porque seran las coordenadas de partida para largar la equilibracion con los "run.sh"

cp $Eqfiles/mdin_eq.$x.crd inp.crd
cp $Eqfiles/mdin_eq.$x.crd inp.$x.crd
 
 
done
