#!/bin/bash

rm -rf build 
sphinx-apidoc -f -o source ../topiary
make html
