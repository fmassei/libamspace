#!/bin/bash
wget www.celestrak.com/NORAD/elements/master.asp -O tmp.html
grep \.txt tmp.html | sed -s "s/[^\"]*\"\([^\.]*\.txt\)\".*/\1/" > list.txt
while read line
do
    wget www.celestrak.com/NORAD/elements/$line -O $line.atxt
done < list.txt
cat *.atxt > all.txt
rm *.atxt
rm list.txt
rm tmp.html
