#!/usr/bin/env bash

git clone https://github.com/biocore/deblur deblur-upstream
nosetests --with-doctest deblur-upstream
rm -rf deblur-upstream

git clone https://github.com/qiime2/q2-deblur q2_deblur
py.test --pyargs q2_deblur
rm -rf q2_deblur
