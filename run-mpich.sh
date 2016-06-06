#!/bin/bash
# EXEMPLE PER A 3 PROCESSOS MPI I THREADS
# Recordeu que el nombre total de processos que s'ha de demanar pel procés és
# igual al nombre de processos MPI x 4.
# Per tant, si per exemple si voleu 3 processos MPI -> 3 MPI x 4 = 12

## Especifica l’intèrpret de comandes pel treball (Bash)
#$ -S /bin/bash

## Exporta totes les variables d’entorn al context del treball
#$ -V

## Executa el treball desde el directori de treball actual.
#$ -cwd

## La cua sobre la que es vol llançar (high.q, low.q, all.q)
#$ -q high.q

## Entorn de programació paral·lel pel treball i numero de processos
#$ -pe mpich-smp 32

## El nom del treball (opcional)
#$ -N TyrionLannister

## Envia un correu quan el treball comença i quan acaba (opcional)
#$ -m be

## L'adreça de correu on enviar les notificacions (ha de ser de la UdL)
#$ -M msf7@alumnes.udl.cat

## Els directoris de la sortida normal i d'error (opcional)
#$ -o $HOME
#$ -e $HOME

MPICH_MACHINES=$TMPDIR/mpich_machines
cat $PE_HOSTFILE | awk '{print $1":1"}' > $MPICH_MACHINES


## Aquesta línia es la que realment executa l'aplicació
mpiexec -f $MPICH_MACHINES -n $NHOSTS ./convolution $1 $2 $3 > $4


rm -rf $MPICH_MACHINES

