#!/bin/bash

rm -rf build 
sphinx-build source/ build/
pdoc3 ../topiary --html -o build/api --force
