for file in datasets/*; do
    { time python2 ogenmap.py < $file ; } &> oresults/`basename $file`
done;
