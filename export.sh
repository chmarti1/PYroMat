#!/bin/bash

rm -rvf dist build

python3 setup.py sdist --format=zip,gztar,bztar
python3 setup.py bdist

#twine upload dist/*.zip
