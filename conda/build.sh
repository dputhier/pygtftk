#!/bin/bash

git clone  git@github.com:dputhier/pygtftk.git
cd pygtftk
pip -r requirements.txt
python setup.py install
