#!/bin/bash
set -e


# Install a system package required by our library
yum install zlib-devel -y
yum install bzip2-devel -y
yum install xz xz-devel -y

# Compile wheels
for PYBIN in $(ls --color=none -d1 /opt/python/*/bin| grep -P "(35)|(36)|(37)|(38)"); do
    echo "${PYBIN}"
    echo ""
    "${PYBIN}/pip" install -U pip
    "${PYBIN}/pip" install numpy>=1.10.0
    "${PYBIN}/pip" install cython
    "${PYBIN}/pip" install -r /io/requirements.txt
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
mkdir -p wheelhouse_manylinux
for whl in wheelhouse/pygtftk*.whl; do
    auditwheel repair "$whl" -w wheelhouse_manylinux 2>&1| tee 1>$whl.log
done

rm -f $(ls  wheelhouse/* | grep -v pygtftk)

# Install packages and test
#for PYBIN in `ls --color=none -d1 /opt/python/*/bin| grep -v "34"| grep -v "37"`; do
#    "${PYBIN}/pip" install pygtftk --no-index -f /io/wheelhouse
#    (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
#done
