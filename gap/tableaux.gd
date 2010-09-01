

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

