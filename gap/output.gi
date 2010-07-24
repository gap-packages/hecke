#######################################################################
##  SPECHT - output.gi : printing functions                          ##
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

InstallMethod(PrintObj,[IsAlgebraObjModule],
  function(m) Print(ModuleString(m,false)); end
);

InstallMethod(ViewString,[IsAlgebraObjModule],
  function(m) 
    return Concatenation("<direct sum of ",
      String(Length(m!.parts))," ",m!.module,"-modules>"); 
  end
);

InstallMethod(ViewObj,[IsAlgebraObjModule],
  function(m) Print(ViewString(m)); end
);

InstallMethod(DisplayString,[IsAlgebraObjModule],
  function(m) return ModuleString(m,true); end
);

InstallMethod(Display,[IsAlgebraObjModule],
  function(m) Print(DisplayString(m),"\n"); end
);

InstallMethod(ModuleString,[IsAlgebraObjModule,IsBool],
  function(a,pp) 
    local x, len, n, star, tmp, coefficients, valuation, str;

    str := "";
    if Length(a!.parts) > 1 then
      n:=Sum(a!.parts[1]);
      if not ForAll(a!.parts, x->Sum(x)=n) then
        Error(a!.module,"(x), <x> contains mixed partitions\n\n");
      fi;
    fi;
    if pp then star:=""; ## pretty printing
    else star:="*";
    fi;

    if Length(a!.parts)=0 or a!.coeffs=[0] then
      a!.coeffs:=[ 0 ]; a!.parts:=[ [] ];
      str:=StringFold(str,["0",star,a!.module,"()"]);
    else
      for x in [Length(a!.parts),Length(a!.parts)-1..1] do
        if IsPolynomial(a!.coeffs[x]) then
          tmp := CoefficientsOfUnivariateRationalFunction(a!.coeffs[x]);
          coefficients := tmp[1];
          valuation := tmp[3];
          if Length(coefficients)=1 then
            if valuation=0 then
              if coefficients=[-1] then Append(str," - ");
              elif coefficients[1]<0 then 
                str:=StringFold(str,[coefficients[1],star]);
              else
                if x<Length(a!.parts) then Append(str," + "); fi;
                if coefficients<>[1] then
                  str:=StringFold(str,[coefficients[1],star]);
                fi;
              fi;
            else
              if x<Length(a!.parts) and coefficients[1]>0 then 
                Append(str," + ");
              fi;
              str:=StringFold(str,[a!.coeffs[x],star]);
            fi;
          elif coefficients[Length(coefficients)]<0 
          then str:=StringFold(str,[" - (",-a!.coeffs[x],")", star]);
          else
            if x<Length(a!.parts) then Append(str, " + "); fi;
            str:=StringFold(str,["(",a!.coeffs[x],")", star]);
          fi;
        else
          if a!.coeffs[x]=-1 then Append(str, " - ");
          elif a!.coeffs[x]<0 then str:=StringFold(str,[a!.coeffs[x],star]);
          else
            if x<Length(a!.parts) then Append(str, " + "); fi;
            if a!.coeffs[x]<>1 
              then str:=StringFold(str,[a!.coeffs[x],star]); fi;
          fi;
        fi;
        if pp then 
          tmp := LabelPartition(a!.parts[x]); 
        else 
          tmp := StringPartition(a!.parts[x]);
        fi;
        str:=StringFold(str,[a!.module,"(",tmp,")"]);
      od;
    fi;
    return str;
  end
);

## adding this string if it does not already exist.
InstallMethod(LabelPartition, [IsList], 
  function(mu) local n, m, label, p;
    n:=Sum(mu);
    if n<2 then return String(n); fi;
    label:="";
    for p in Collected(mu) do
      if p[2]=1 then label:=Concatenation(String(p[1]),",",label);
      else label:=Concatenation(String(p[1]),"^",String(p[2]),",",label);
      fi;
    od;
    return label{[1..Length(label)-1]};
  end
);

#F Returns a string for ModuleString() from SpechtParts.labels, adding this 
## string if it does not already exist.
InstallMethod(StringPartition, [IsList], 
  function(mu) local m, string, p;
    if mu=[] or mu=[0] then return "0";
    else return ViewString(mu);#TODO TightStringList(mu);
    fi;
  end
);

