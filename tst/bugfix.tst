# https://github.com/gap-packages/hecke/issues/1
gap> CombineEQuotientECore(2,EQuotient(2,[]),[]);
[  ]

# https://github.com/gap-packages/hecke/issues/12
gap> SemiStandardTableaux([2,2,1],[3,0,2]);
[  ]
gap> SemiStandardTableaux([2,2,2],[3,1,1,1]);
[  ]
gap> SemiStandardTableaux([3,3],[4,1,1]);
[  ]
gap> List(SemiStandardTableaux([2,2,1],[1,2,2]), t -> t![1]);
[ [ [ 1, 2 ], [ 2, 3 ], [ 3 ] ] ]
gap> Length(SemiStandardTableaux([3,2,1]));
33
