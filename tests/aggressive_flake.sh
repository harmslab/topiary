#!/bin/bash

flake8 . --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics
