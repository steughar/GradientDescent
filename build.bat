@echo off

set CommonLinkerFlags= -incremental:yes -opt:ref
IF NOT EXIST ..\build mkdir ..\build
pushd ..\build
del *.pdb > NUL 2> NUL

cl -Z7 -nologo -Fmeval_circle.map ..\eval_circle\source.c ..\eval_circle\eval_circle.c

popd
