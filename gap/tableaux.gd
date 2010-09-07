#######################################################################
##                                                                   ##
##   Hecke - tableaux.gi : Young tableaux                            ##
##                                                                   ##
##     This file contains functions for generating all kinds of      ##
##     young tableaux                                                ##
##                                                                   ##
##     Dmitriy Traytel                                               ##
##     (heavily using the GAP3-package SPECHT 2.4 by Andrew Mathas)  ##
##                                                                   ##
#######################################################################

## Hecke 1.0: June 2010:
##   - initial

BindGlobal("TableauFamily",NewFamily("TableauFamily"));

DeclareCategory("IsTableau", IsPositionalObjectRep);
DeclareCategory("IsSemiStandardTableau", IsTableau);
DeclareCategory("IsStandardTableau", IsSemiStandardTableau);

BindGlobal("TableauType",NewType(TableauFamily, IsTableau));
BindGlobal("SemiStandardTableauType",
  NewType(TableauFamily, IsSemiStandardTableau));
BindGlobal("StandardTableauType",
  NewType(TableauFamily, IsStandardTableau));

MakeDispatcherFunc("Tableau", [[]],[1],[1]); ##Constructor

DeclareOperation("Specht_PrettyPrintTableau",[IsTableau]); ## use ViewObj instead

MakeDispatcherFunc("SemiStandardTableaux",[[IsList],[]],[2,1],[2,1]);
MakeDispatcherFunc("StandardTableaux",[[]],[1],[1]);

DeclareOperation("ConjugateTableau", [IsTableau]);
DeclareOperation("TypeTableau", [IsTableau]);
DeclareOperation("ShapeTableau", [IsTableau]);

