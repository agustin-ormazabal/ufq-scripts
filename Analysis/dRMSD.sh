ls dh*sel*out -l   | awk '{print $9}' > dihedros_sel

nsis=$(cat dihedros_sel | wc |awk '{print $1}')

sel=$(sed -n 1p dihedros_sel | awk '{print $1}')
nlineas_sel=$(cat $sel | wc | awk '{print $1}')
tmp=$(cat $sel | wc | awk '{print $2}')
ncolumnas=$(echo $tmp / $nlineas_sel | bc -l | awk '{printf("%0.0f",$0);}')

for x in $(seq 1 1 $nsis) ; do

sel=$(sed -n "$x"p dihedros_sel | awk '{print $1}')
nlineas_sel=$(cat $sel | wc | awk '{print $1}')

name=$(echo "$sel" | sed 's/_/ /g' | awk '{print $2}')

cat > dRMSD_LvsR.f << EOF 
          program dRMSD

          parameter(nfilasA=$nlineas_sel)
          parameter(ncolumnas=$ncolumnas)
          parameter(ncolumnas2=$ncolumnas*2)
         integer lfil,lcol,rfil,rcol,par,x,y,z
        real*8 vL(nfilasA,ncolumnas)
        real*8 rmsdL,rmsdR,dRMSD_si(nfilasA,nfilasA)
        real*8 cosxl(nfilasA,ncolumnas),senxl(nfilasA,ncolumnas)
        real*8 Lval(nfilasA,ncolumnas2)
        real*8 Lsumval,Rsumval
        real*8 diferencia(nfilasA,nfilasA)
        real*8 cont,proof
        character*60 fileL,fileR,outL,outR,outL2,outR2,fixedL,fixedR

!  Los inputs son archivos de cpptraj que contienen el valor de
! cada ángulo dihedro del sistema a lo largo de la trayectoria


          fileL="$sel"

          fixedL="fixedL_sel.out"

          outL="dRMSDL_sel_abs_$name"

          outL2="dRMSDL2"
          outR2="dRMSD_plot"

          Lval1 = 0
          Rval1 = 0

!  Abre los archivos de cpptraj y los lee por columnas y filas
! y los guarda en un archivo

          open(15,file=fileL)

          do lfil = 1,nfilasA
 
             read(15,*) vL(lfil,:)

          enddo

          close(15)
          close(10)

        open(16,file=fixedL)

         do lfil = 1,nfilasA

          do lcol = 1,ncolumnas

            cosxl(lfil,lcol) = cos(vl(lfil,lcol))
            senxl(lfil,lcol) = sin(vl(lfil,lcol))

          enddo

          write(16,23) (cosxl(lfil,lcol),senxl(lfil,lcol),
     &lcol=2,ncolumnas)

         enddo

         close(16)

!  Lee el archivo con los sen y cos y calcula el rmsd

         open(16,file=fixedL)

          do lfil = 1,nfilasA

             read(16,23) Lval(lfil,:)

          enddo

         close(16)
!!!
       open(17,file=outL)

        do z = 1,nfilasA

         do x = 1,nfilasA

          diferencia(z,x) = 0

          do y = 1,ncolumnas2

          diferencia(z,x) = diferencia(z,x) + (Lval(z,y) - Lval(x,y))**2

          enddo

           dRMSD_si(z,x) = sqrt(diferencia(z,x)/ ncolumnas2)

         enddo

        enddo

        do z = 1,nfilasA

          write(17,*) (dRMSD_si(z,x),x=1,nfilasA)

        enddo

       close(17)


   21    format(102(1x,f8.4))
   22    format(204(f8.4,5x))
   23    format(204(f8.4,5x))
   24    format(f7.4,1x,f7.4)
   25    format(f6.4)
   26    format(i4,1x,f10.4)
   27    format(i4,1x,f14.12,1x,i4)
   28    format(f8.6)
         end program
EOF

gfortran dRMSD_LvsR.f -o dRMSD_LvsR.exe
./dRMSD_LvsR.exe

done

