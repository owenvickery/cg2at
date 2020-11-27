PKG_NAME=cg2at
USER=stansfeld_rg

OS=linux-64
mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld
export VERSION=`date +%Y.%m.%d`
export package=$PKG_NAME-`date +%Y.%m.%d`-0
conda build . 
echo anaconda -t $CONDA_UPLOAD_TOKEN upload -l main --force -u $USER $CONDA_BLD_PATH/$OS/$PKG_NAME-`date +%Y.%m.%d`-0.tar.bz2

anaconda -t $CONDA_UPLOAD_TOKEN upload -l main --force -u $USER $CONDA_BLD_PATH/$OS/$PKG_NAME-`date +%Y.%m.%d`-0.tar.bz2 
