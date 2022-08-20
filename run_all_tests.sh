#!/bin/bash

echo "Running flake8"
flake_test=`flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics`
if [[ "${flake_test}" != 0 ]]; then
    echo "flake failed"
    exit
fi

rm -rf reports
mkdir reports
mkdir reports/junit
mkdir reports/coverage
mkdir reports/badges

echo "Running flake8, aggressive"
flake8 . --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics > reports/flake.txt

echo "Checking for test file completeness"
./tests/completeness_crawler.py topiary tests bin/ > reports/completeness-crawler.txt
grep MISSING reports/completeness-crawler.txt

echo "Running coverage.py"
coverage erase
coverage run --branch -m pytest --run-raxml --run-generax --junit-xml=reports/junit/junit.xml

echo "Generating reports"
coverage html
mv htmlcov reports

coverage xml
mv coverage.xml reports/coverage/coverage.xml

genbadge tests
sleep 1
genbadge coverage
sleep 1

wget https://github.com/harmslab/topiary/actions/workflows/python-app.yml/badge.svg -O ghwf.svg
wget https://readthedocs.org/projects/topiary-asr/badge/?version=latest -O rtd.svg

mv *.svg docs/badges
