#!/bin/bash
set -e


# Install a system package required by our library
yum install zlib-devel

# Compile wheels
for PYBIN in `ls --color=none -d1 /opt/python/*/bin| grep -v "34"| grep -v "37"`; do echo $PYBIN; done
    echo "${PYBIN}"
    echo ""
    "${PYBIN}/pip" install -U pip
    "${PYBIN}/pip" install numpy>=1.10.0
    "${PYBIN}/pip" install -r /io/requirements_developers.txt
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/pygtftk*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/ 2>&1| tee 1>$whl.log
done

# Install packages and test
#for PYBIN in /opt/python/*/bin/; do
#    "${PYBIN}/pip" install pygtftk --no-index -f /io/wheelhouse 2>&1| tee 1>$whl.log
#    (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
#done
