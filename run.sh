for file in datasets/*; do
    { time python2 ogenmap.py < $file > /dev/null ; } 2> oresults/`basename $file`
done;
