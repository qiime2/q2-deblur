#!/usr/bin/env bash

git clone https://github.com/biocore/deblur deblur-upstream
nosetests --with-doctest deblur-upstream
rm -rf deblur-upstream

py.test --pyargs q2_deblur
