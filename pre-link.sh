printenv > $CONDA_PREFIX/test_cg2at_env_pre.txt
ln -fs $CONDA_PREFIX/pkgs/$PKG_NAME/info/recipe/cg2at $CONDA_PREFIX/bin/cg2at ||  echo "ln -fs $CONDA_PREFIX/pkgs/$PKG_NAME/info/recipe/cg2at $CONDA_PREFIX/bin/cg2at  wrong" > $PREFIX/.messages

echo $PREFIX >  $CONDA_PREFIX/test_cg2at_pre.txt
echo $PKG_NAME >>  $CONDA_PREFIX/test_cg2at_pre.txt
echo $PKG_VERSION >>  $CONDA_PREFIX/test_cg2at_pre.txt
echo $PKG_BUILDNUM >>  $CONDA_PREFIX/test_cg2at_pre.txt
