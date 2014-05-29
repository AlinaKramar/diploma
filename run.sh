for file in datasets/*; do
    { time python2 genmap.py < $file ; } &> results/`basename $file`
done;
