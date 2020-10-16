if [[ ! -z "${CONDA_PREFIX}" ]]
then
	if [[ ${CONDA_PREFIX} != *"/home/travis"* ]]
	then
 		ln -fs $PREFIX/info/recipe/database/script_files/cg2at.py $CONDA_PREFIX/bin/CG2AT
 	fi
fi