for d in frame_????
do
cd $d
mopicgen casscf.molden --orient "reset;center {0,0,0} ; rotate z 0; rotate y 0; rotate x -45;" --fracmos
bash run.sh 
cd ..
done
