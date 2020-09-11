scenes="1 2 3 4 5 6 7 8 9 10"
modes="dem epnp p4p"
vvss="true false"
for s in $scenes
do
    for m in $modes
    do
    for vvs in $vvss
    do
    echo "$c / $m / $vvs"
    ../build/rose scene $s estim $m control rose kill false plot create compare_estim true vvs $vvs noise 0
    done
    done
done

rm estim_summary.yaml
python estims3.py
