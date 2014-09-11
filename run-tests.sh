#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#echo $DIR
MYVALUE="$DIR/tests/test_data"
#echo "$MYVALUE"
export PYTHONPATH=${PYTHONPATH}:$DIR:$MYVALUE
#echo ${PYTHONPATH}
cd tests
python test_ModeRNA_all.py