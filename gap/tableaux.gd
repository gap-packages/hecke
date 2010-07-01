

BindGlobal("TableauFamily",NewFamily("TableauFamily"));

DeclareCategory("IsTableau", IsPositionalObjectRep);
DeclareCategory("IsSemiStandardTableau", IsTableau);
DeclareCategory("IsStandardTableau", IsSemiStandardTableau);

BindGlobal("TableauType",NewType(TableauFamily, IsTableau));
BindGlobal("SemiStandardTableauType",
  NewType(TableauFamily, IsSemiStandardTableau));
BindGlobal("StandardTableauType",
  NewType(TableauFamily, IsStandardTableau));

DeclareOperation("Tableau", [IsList]); ##Constructor

DeclareOperation("Specht_PrettyPrintTableau",[IsTableau]); ## use ViewObj instead

DeclareOperation("SemiStandardTableaux",[IsList,IsList]);
DeclareOperation("SemiStandardTableaux",[IsList]);
DeclareOperation("StandardTableaux",[IsList]);

DeclareOperation("ConjugateTableau", [IsTableau]);
DeclareOperation("TypeTableau", [IsTableau]);
DeclareOperation("ShapeTableau", [IsTableau]);
