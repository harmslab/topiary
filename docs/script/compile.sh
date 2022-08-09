#!/bin/bash

cd source
python ../script/build_api_docs.py
cd ..

rm -rf build 
make html
