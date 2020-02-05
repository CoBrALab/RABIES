#credit to Joanes Grandjean https://github.com/grandjeanlab/MouseMRIPrep
#Motion plots, assuming the 6 motion parameters are in $mc_par, with rotation second
mc_par=$1
prefix=$2

fsl_tsplot -i $mc_par -t 'antsMotionCorr estimated rotations (radians)' -u 1 --start=4 --finish=6 -a roll,pitch,yaw -w 640 -h 144 -o rot.png
fsl_tsplot -i $mc_par -t 'antsMotionCorr estimated translations (mm)' -u 1 --start=1 --finish=3 -a dS,dL,dP -w 640 -h 144 -o trans.png
pngappend trans.png - rot.png ${prefix}_motion_traces.png
