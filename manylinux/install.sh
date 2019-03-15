echo "Creating working directory."
mkdir -p /io && cd /io

echo "Cloning git repository."
git clone --branch develop --single-branch  https://github.com/dputhier/pygtftk.git .

echo "Instructing pygtftk we are doing a release."
touch release_in_progress

echo "Building manylinux compliant pip packages."
cd manylinux
chmod u+x build_wheels.sh
bash ./build_wheels.sh &> /tmp/log
cp -r ./* /tmp
