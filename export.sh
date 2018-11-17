#!/bin/bash

rm -rvf dist build

python setup.py sdist --format=zip,gztar,bztar
twine upload dist/*.zip
