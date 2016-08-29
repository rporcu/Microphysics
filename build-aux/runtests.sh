#!/bin/bash -ex

if [ -f runtests.sh ]; then
    exec ./runtests.sh
fi

if [ -n "${MPICMD}" ]; then
    ${MPICMD} ./mfix${EXEEXT}
else
    ./mfix${EXEEXT}
fi

post_dats=AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    numdiff -a 0.000001 -r 0.05 ${test_post_file} $(basename ${test_post_file}) || echo "Post results differ"
done
