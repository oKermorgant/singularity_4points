

methods="dem p4p upnp epnp"
refs="none vvs lm"
# noises="0 0.0001 0.00001"
#scenes="1 2 3 4 5 6 7 8 9 10"

scenes="1"
noises="0.00001"
methods="epnp" 
refs="none vvs lm"

dataPath=/home/olivier/Results/singularity
if [ "$HOSTNAME" = sterne ]; then
dataPath=/home/olivier/Olivier/boulot/Results/singularity
fi
for s in $scenes
do
for method in $methods
do
    for ref in $refs
    do
     for noise in $noises
     do
#     echo "dataPath $dataPath scene $s estim $method ref $ref kill false plot none noise $noise"
         ../build/sphere dataPath $dataPath scene $s estim $method ref $ref kill false plot none noise $noise &
         done
    done
done
done
