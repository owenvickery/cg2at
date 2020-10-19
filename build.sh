# if [[ ! -z "${CONDA_PREFIX}" ]]
# then
# 	touch $PREFIX/test_1.dat
# 	if [[ ${CONDA_PREFIX} != *"/home/travis"* ]]
# 	then
#  		ln -fs $PREFIX/info/recipe/database/script_files/cg2at.py $CONDA_PREFIX/bin/CG2AT
#  	fi
# fi
printenv > $CONDA_PREFIX/test_cg2at_env.txt
# if [[ ${CONDA_PREFIX} != "/home/travis/"* ]] 
# then
# 	echo "ln -fs $CONDA_PREFIX/bin/$package/cg2at $CONDA_PREFIX/bin/cg2at  correct" > $CONDA_PREFIX/test_cg2at.txt
# 	ln -fs $CONDA_PREFIX/bin/$package/cg2at $CONDA_PREFIX/bin/cg2at

# else
# 	echo "ln -fs $CONDA_PREFIX/bin/$package/cg2at $CONDA_PREFIX/bin/cg2at  wrong" > $CONDA_PREFIX/test_cg2at.txt
# fi

{ ln -fs $CONDA_PREFIX/bin/$package/cg2at $CONDA_PREFIX/bin/cg2at && } || { echo "ln -fs $CONDA_PREFIX/bin/$package/cg2at $CONDA_PREFIX/bin/cg2at  wrong" > $CONDA_PREFIX/test_cg2at.txt }
