
PKG_NAME=cg2at
USER=stansfeld_rg

OS=noarch
mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld
export VERSION=1.0.0   # must match the version in meta.yaml
export package=$PKG_NAME-$VERSION-0
conda build . 
anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l main $CONDA_BLD_PATH/$OS/$PKG_NAME-$VERSION-0.tar.bz2 --force
