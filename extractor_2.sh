echo -e \\e[32m  " "
echo " Reading input.txt file  "
input_file=input.txt
selected_calculation=$(grep "calculation" $input_file | awk '{print $3}' )
sampling_parameter=$(grep "sampling_parameter" $input_file | awk '{print $3}')
selected_sampling=$(grep "type_of_sampling" $input_file | awk '{print $3}')
sampling_size=$(grep "sampling_size" $input_file | awk '{print $3}')
dir_name=$(grep "directory_name" $input_file | awk '{print $3}')

case $# in
   0)
      if [ $selected_calculation = 1 ]
      then
         echo -e \\e[91m "No input data file provided"
         exit
      fi
   ;;
   1)
      if ! [ -f $1 ]
      then
         echo  -e \\e[91m " $1 file does not exists"
         exit
      fi
   ;;
esac
######################################################################
function redraw_progress_bar { # int barsize, int base, int i, int top
    local barsize=$1
    local base=$2
    local current=$3
    local top=$4
    local j=0
    local progress=$(( ($barsize * ( $current - $base )) / ($top - $base ) ))
    echo -n "["
    for ((j=0; j < $progress; j++)) ; do echo -n '='; done
    echo -n '=>'
    for ((j=$progress; j < $barsize ; j++)) ; do echo -n ' '; done
    echo -ne "] $(( $current )) / $top " $'\r'
}
######################################################################

#######################################################
#######################################################
if [ $selected_calculation = 1 ]
then
   echo " Extracting data from $1"
   echo " "
   ## Determina número de átomos
   Nat=$(grep "| Number of atoms" $1  | awk '{print $6}')
   ## Determina número de configuraciones (pasos)
   Nconf=$(grep "Self-consistency cycle converged" $1  | wc -l)
   ## Obtiene coordenadas iniciales
   grep -A$((5+$(echo $Nat))) "Parsing geometry.in (first pass over file, find array dimensions only)." $1  | grep "atom" > initial_configuration.fhi
   ## Extrae lista de energias
   echo  -e \\e[32m "Extracting uncorrected energy data"
   echo " "
   pv $1 | grep  "Total energy uncorrected"  | awk '{print $6}' > energias_uncorrected # puede ser: | Total energy uncorrected  | Total energy corrected | Electronic free energy ;
   echo " "
   echo  -e \\e[32m "Extracting corrected energy data"
   echo " "
   pv $1 | grep  "Total energy corrected"  | awk '{print $6}' > energias_corrected
   echo " "
   echo  -e \\e[32m "Extracting free energy data"
   echo " "
   pv $1 | grep -A2  "Total energy uncorrected" | grep "Electronic free energy" | awk '{print $6}' > energias_free
   echo " "
   echo  -e \\e[32m "Collecting into  energies.csv file"
   echo " " ; echo " "
   paste energias_corrected energias_uncorrected  > temp ; paste temp energias_free > energies.temp
   awk '{print $1","$2","$3}' energies.temp > energies.csv; rm temp energies.temp
   # Separa momento dipolar, uno por configuración
   echo " "
   echo  -e \\e[32m "Extracting dipolar moment components"
   echo " "
   pv $1 | grep "Total dipole moment"  | awk '{print $7 " " $8 " " $9}' > dipolos_vec
   echo " "
   echo  -e \\e[32m "Extracting dipolar moment magnitudes"
   echo " "
   grep "Absolute dipole moment" $1 | awk '{print $6}' > dipolos_abs
   echo " "
   echo  -e \\e[32m "Collecting data into dipoles.csv"
   echo " " ; echo " "

   paste dipolos_vec dipolos_abs > dipole.temp ; rm dipolos_vec dipolos_abs
   awk '{print $1","$2","$3","$4}' dipole.temp > dipoles.csv; rm dipole.temp
   # Junta los datos de energias y dipolos
   #paste energias dipolos > energies_dipoles.csv ; rm energias dipolos
   ## Separa las posiciones atómicas
   echo  -e \\e[32m "Extracting atomic positions"
   echo " "
   pv $1 | grep -A$((1+2*$(echo $Nat))) "Atomic structure (and velocities) as used in the preceding time step:"  | grep "atom" > posiciones
   # Hace falta separar todo el archivo en individuales
   ## Separa las velocidades atómicas
   echo " "
   echo  -e \\e[32m "Extracting atomic velocities" ; echo " "
   pv $1 | grep -A$((1+2*$(echo $Nat))) "Atomic structure (and velocities) as used in the preceding time step:"  | grep "velocity" > velocidades
   # Hace falta separar todo el archivo en individuales
   # Separa las fuerzas atómicas unidades  [eV/Ang]:
   echo " "
   echo  -e \\e[32m "Extracting atomic forces" ; echo " "
   pv $1 | grep -A40 "Total atomic forces (unitary forces cleaned)"  | grep "|" > fuerzas
   # Idem: Capaz k lo podríamos hacer en el mismo loop
   echo " " ; echo " " ; echo " "
   echo  -e \\e[32m "Preparing formatted files" ; echo " "
echo "x,y,z" > fuerzas.csv
echo "x,y,z" > posiciones.csv
echo "E" > energias_torch.csv
echo "S" > atomos_simbolos.csv

   for ((i=$Nat;i<=$(($Nconf*$Nat));i=i+$Nat))
   do
      j=$(($i/$Nat))
      head -$i posiciones   | tail -$Nat > positions_$j.fhi
      head -$i velocidades  | tail -$Nat > velocities_$j.fhi
      head -$i fuerzas  | tail -$Nat > forces_$j.fhi
      redraw_progress_bar 50 1 $j $Nconf
      ###################################################################
#       echo $Nat >> Ti18C18.xyz
#       echo " " >> Ti18C18.xyz
#      echo "begin"> coords_$j
#      echo "comment configuration ${j} of $Nconf" >> coords_$j
      # Lattice section is skipped because this is not a periodic system
#/**************************************************************************
#split -l$ --numeric-suffixes=1 --suffix-length=$suf_len --additional-suffix=".fhi" posiciones "posiciones"
#split -l$Nat --numeric-suffixes=1 --suffix-length=$suf_len --additional-suffix=".fhi" fuerzas "forces"
#      yes "atom " | head -$nl > atoms
#      yes " 0.0 0.0 " | head -$nl > zeros
# ls positions* | parallel awk '{print $2" "$3" "$4" "$5}'  > temp*{corresponding} pos_temp.{}
# ls forces* | parallel awk '{print $2" "$3" "$4" "$5}'  > temp*{corresponding} for_temp.{}
# ls pos_temp* | parallel paste atoms temp
#********************************************************************************/
      cat positions_$j.fhi | awk '{print $2","$3","$4}'  > temp ; cat forces_$j.fhi | awk '{print $3","$4","$5}' > temp2
      cat temp >> posiciones.csv
      cat temp2 >> fuerzas.csv
      if [ $j == 1 ]
      then
      cat positions_$j.fhi | awk '{print $5}'  >> atomos_simbolos.csv
      fi
#      nl=$(wc -l temp | awk '{print $1}')
      ## Esto podria ir hasta arriba, asi evitamos repetirlo #
#      yes "atom " | head -$nl > atoms
#      yes " 0.0 0.0 " | head -$nl > zeros
      ########################################################
      ########################################################
 #     paste temp temp2 > coords_temporal
      energy=$(sed "${j}q;d" energias_uncorrected)
      echo "$energy" >> energias_torch.csv
#      echo "charge 0.0" >> coords_$j
#      echo "end" >> coords_$j
#      cat coords_temporal | tr '\t' ' ' >> Ti18C18.xyz
  #    rm coords_temporal
      ###################################################################

      #   echo "Preparing step $j/$Nconf"
      echo -ne "[\r"
   done
   echo " " ; echo " " ; echo " " ;
rm  velocities* forces* positions* velocidades energias_uncorrected  energies.csv energias_corrected energias_free fuerzas  initial_configuration.fhi  posiciones  temp temp2  dipoles.csv
cat atomos_simbolos.csv | tr 'C' '6' | tr 'T' '2' | tr 'i' '2' > temporal
mv temporal atomos_simbolos.csv
fi

