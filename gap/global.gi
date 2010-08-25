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

InstallMethod(\*,"function composition",[IsFunction,IsFunction],
  function(f,g) 
    return function(x) return f(g(x)); end; 
  end
);

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

InstallMethod(Zip,[IsList,IsList],
  function(l1,l2) local i, res;
    res:=[];
    for i in [1..Length(l1)] do
      res[i]:=[l1[i],l2[i]];
    od;
    return res;
  end
);
