#!/bin/bash

rm -rvf dist build

python3 setup.py sdist --format=zip,gztar,bztar
twine upload dist/*.zip
