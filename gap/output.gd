#######################################################################
##  SPECHT - output.gd : printing functions                          ##
##                                                                   ##
##     A GAP package for calculating the decomposition numbers of    ##
##     Hecke algebras of type A (over fields of characteristic       ##
##     zero). The functions provided are primarily combinatorial in  ##
##     nature. Many of the combinatorial tools from the (modular)    ##
##     representation theory of the symmetric groups appear in the   ##
##     package.                                                      ##
##                                                                   ##
##     These programs, and the enclosed libraries, are distributed   ##
##     under the usual licensing agreements and conditions of GAP.   ##
##                                                                   ##
##     Dmitriy Traytel                                               ##
##     (heavily using the GAP3-version by Andrew Mathas)             ##
##                                                                   ##
#######################################################################

DeclareOperation("ModuleString",[IsAlgebraObjModule,IsBool]);
DeclareOperation("LabelPartition",[IsList]);
DeclareOperation("StringPartition",[IsList]);
