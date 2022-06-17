#!/usr/bin/env bash

git clone https://github.com/biocore/deblur deblur-upstream
cd deblur-upstream
nosetests --with-doctest
cd ../
rm -rf deblur-upstream
