
PKG_NAME=cg2at
USER=stansfeld_rg

OS=noarch
mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld
export VERSION=`date +%Y.%m.%d`
export package=$PKG_NAME-`date +%Y.%m.%d`-0
conda build . 
anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l main $CONDA_BLD_PATH/$OS/$PKG_NAME-`date +%Y.%m.%d`-0.tar.bz2 --force
