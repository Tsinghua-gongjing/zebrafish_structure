
./paris \
    GenDuplexGroup \
    --log log.2017-7-6.txt \
    --error error.2017-7-6.txt \
    --inReadPairFile /Users/lee/Desktop/paris-simple/c++/data/dupleGroup.txt \
    --outDGFile /Users/lee/Desktop/paris-simple/c++/data/find_dupleGroup2.txt \
    --minOverlap 5 \
    --multipleDG no

./paris \
    CollapseDuplexGroup \
    --log log.2017-7-6.txt \
    --error error.2017-7-6.txt \
    --inDGFile /Users/lee/Desktop/paris-simple/c++/data/find_dupleGroup.txt.sorted \
    --outCollapseFile /Users/lee/Desktop/paris-simple/c++/data/collapsed_dupleGroup2.txt \
    --maxGap 10 \
    --maxDGOverhang 30

