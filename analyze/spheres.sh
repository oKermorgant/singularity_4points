modes="dem p4p upnp epnp"
refs="none vvs lm"
noises="0 0.0001 0.00001"
#scenes="1 2 3 4 5 6 7 8 9 10"
scenes="2 12"

dataPath=/home/olivier/Results/singularity
if [ "$HOSTNAME" = sterne ]; then
dataPath=/home/olivier/Olivier/boulot/Results/singularity
fi
for s in $scenes
do
for m in $modes
do
    for ref in $refs
    do
     for noise in $noises
     do
        ../build/sphere dataPath $dataPath scene $s estim $m kill false plot false display create ref $ref noise $noise &
         done
    done
done
done

