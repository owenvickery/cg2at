# if [[ ! -z "${CONDA_PREFIX}" ]]
# then
# 	touch $PREFIX/test_1.dat
# 	if [[ ${CONDA_PREFIX} != *"/home/travis"* ]]
# 	then
#  		ln -fs $PREFIX/info/recipe/database/script_files/cg2at.py $CONDA_PREFIX/bin/CG2AT
#  	fi
# fi
if [[ ${CONDA_PREFIX} != "/home/travis/"* ]]; then
	ln -fs $CONDA_PREFIX/bin/$package/cg2at $CONDA_PREFIX/bin/cg2at
fi
# echo $CONDA_PREFIX/bin/CG2AT > $CONDA_PREFIX/test_cg2at.txt
# printenv > $CONDA_PREFIX/test_cg2at.txt
# $PATH >> $CONDA_PREFIX/test_cg2at.txt
# echo $PREFIX/pre >> $CONDA_PREFIX/test_cg2at.txt
# export LD_RUN_PATH = $ORIGIN