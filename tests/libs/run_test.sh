#!/bin/bash

../../libs/snippets/merger_ex ./merge_data/fof.txt 1 1 > tmp11.txt
../../libs/snippets/merger_ex ./merge_data/fof.txt 1 2 > tmp12.txt
../../libs/snippets/merger_ex ./merge_data/fof.txt 2 1 > tmp21.txt

diff ./merge_data/merge_a1r1.txt tmp11.txt
if [ $? -ne 0 ]; then
    exit 1
fi

diff ./merge_data/merge_a1r2.txt tmp12.txt
if [ $? -ne 0 ]; then
    exit 1
fi

diff ./merge_data/merge_a2r1.txt tmp21.txt
if [ $? -ne 0 ]; then
    exit 1
fi

rm tmp11.txt tmp12.txt tmp21.txt

../../libs/snippets/bitmatrix_ex > /dev/null

../../libs/snippets/skreader_ex ./skreader_data/sk_part > tmpSk0.txt
diff ./skreader_data/sk_part/superks0.txt tmpSk0.txt
rm tmpSk0.txt
