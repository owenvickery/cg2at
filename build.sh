# if [[ ! -z "${CONDA_PREFIX}" ]]
# then
# 	touch $PREFIX/test_1.dat
# 	if [[ ${CONDA_PREFIX} != *"/home/travis"* ]]
# 	then
#  		ln -fs $PREFIX/info/recipe/database/script_files/cg2at.py $CONDA_PREFIX/bin/CG2AT
#  	fi
# fi
export LD_RUN_PATH = $ORIGIN