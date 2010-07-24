#######################################################################
##  SPECHT - specht.g : the kernel of SPECHT                         ##
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

InstallMethod(FoldLeft,[IsFunction,IsObject,IsList],
  function(f,i,l) local acc, x;
    acc := i;
    for x in l do
      acc := f(acc,x);
    od;
    return acc;
  end
);

InstallMethod(FoldAfterMapLeft,[IsFunction,IsFunction,IsObject,IsList],
  function(f,g,i,l) local acc, x;
    acc := i;
    for x in l do
      acc := f(acc,g(x));
    od;
    return acc;
  end
);

InstallMethod(StringFold,[IsString,IsList],
  function(str,l)
    return FoldAfterMapLeft(Concatenation,String,str,l);
  end
);
