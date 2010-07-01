###################################################################
##  induce.g: mostly handling functions                          ##
##                                                               ##
##     In a future release of SPECHT this file will contain all  ##
##     of the functions for inducing and restricting and their   ##
##     q-analogues, together with other funcitons needed for     ##
##     calculations in the Fock space(s) etc.                    ##
##                                                               ##
##     Andrew Mathas                                             ##
##                                                               ##
###################################################################

## Change log
## 2.2: June 1996 
##        o  moved these functions out of specht.g

## General handling functions 

## Usage: Specialized(x) or Specialized(x,a); defaults to a=1 (**undocumented)
Specialized:=function(arg) local x;
   x:=arg[1];
   if IsRec(x) and IsBound(x.operations) 
   and IsBound(x.operations.Specialized) then
     if Length(arg)=1 then return x.operations.Specialized(x,1);
     else return x.operations.Specialized(x, arg[2]);
     fi;
   else Error("Specialized(<x>), don't know how to specialize <x>.\n");
   fi;
end;

InducedModule:=function(arg)
  if arg=[] or arg[1] = false then  return false;
  elif IsBound(arg[1]) and IsRec(arg[1]) and IsBound(arg[1].operations)
  and IsBound(arg[1].operations.InducedModule) then
    return arg[1].operations.InducedModule(arg[1],arg{[2..Length(arg)]});
  else Error("InducedModule(<arg>), don't know how to induce <arg>");
  fi;
end;

SInducedModule:=function(arg)
  if arg=[] or arg[1] = false then  return false;
  elif IsBound(arg[1]) and IsRec(arg[1]) and IsBound(arg[1].operations)
  and IsBound(arg[1].operations.InducedModule) then
    return arg[1].operations.SInducedModule(arg[1],arg{[2..Length(arg)]});
  else Error("SInducedModule(), don't know how to S-induce ", arg);
  fi;
end;

RestrictedModule:=function(arg)
  if arg=[] or arg[1] = false then  return false;
  elif IsBound(arg[1]) and IsRec(arg[1]) and IsBound(arg[1].operations)
  and IsBound(arg[1].operations.RestrictedModule) then
    return arg[1].operations.RestrictedModule(arg[1],arg{[2..Length(arg)]});
  else Error("RestrictedModule(), don't know how to restrict ", arg);
  fi;
end;

SRestrictedModule:=function(arg)
  if arg=[] or arg[1] = false then  return false;
  elif IsBound(arg[1]) and IsRec(arg[1]) and IsBound(arg[1].operations)
  and IsBound(arg[1].operations.SRestrictedModule) then
    return arg[1].operations.SRestrictedModule(arg[1],arg{[2..Length(arg)]});
  else Error("SRestrictedModule(), don't know how to s-restrict ", arg);
  fi;
end;

InnerProduct:=function(a,b)
  if IsRec(b) and IsBound(b.operations) and IsBound(b.operations.InnerProduct) 
  then return b.operations.InnerProduct(a,b);
  else Error("InnerProduct(<a>,<b>), don't know how to take the inner product",
             " of <a> or <b>");
  fi;
end;

Coefficient:=function(arg)
  if IsRec(arg[1]) and IsBound(arg[1].operations) 
  and IsBound(arg[1].operations.Coefficient) then
    return arg[1].operations.Coefficient(arg[1],Flat(arg{[2..Length(arg)]}));
  else
    Error("Coefficient(<x>,<p>), don't know how to find the coefficient of",
          " <p> in <x>");
  fi;
end;

PositiveCoefficients:=function(x)
  if IsRec(x) and IsBound(x.operations) 
  and IsBound(x.operations.PositiveCoefficients) then
    return x.operations.PositiveCoefficients(x);
  else
    Error("PositiveCoefficients(<x>), don't know how to test <x> for ",
          "positive coefficients");
  fi;
end;

IntegralCoefficients:=function(x)
  if IsRec(x) and IsBound(x.operations) 
  and IsBound(x.operations.IntegralCoefficients) then
    return x.operations.IntegralCoefficients(x);
  else
    Error("IntegralCoefficients(<x>), don't know how to test <x> for ",
          "positive coefficients");
  fi;
end;

##################################################################

#F Returns a string of the form "l1,l2,...,lk" for the list [l1,l2,..,lk]"
## which contains no spaces, and has first element 1.
TightStringList:=function(list) local s, l;
  if list=[] then return ""; fi;
  s:=String(list[1]);
  for l in [2..Length(list)] do 
    Add(s,','); 
    Append(s,String(list[l])); 
  od;
  return s;
end;

##################################################################

## adding this string if it does not already exist.
LabelPartition:=function(mu) local n, m, label, p;
  n:=Sum(mu);
  if n<2 then return String(n); fi;
  label:="";
  for p in Collected(mu) do
    if p[2]=1 then label:=Concatenation(String(p[1]),",",label);
    else label:=Concatenation(String(p[1]),"^",String(p[2]),",",label);
    fi;
  od;
  return label{[1..Length(label)-1]};
end;

#F Returns a string for PrintModule() from SpechtParts.labels, adding this 
## string if it does not already exist.
StringPartition:=function(mu) local m, string, p;
  if mu=[] or mu=[0] then return "0";
  else return TightStringList(mu);
  fi;
end;

#F The function used by Specht to decide the format of partitions when
## printing; see SpechtPrettyPrint).
SpechtPrintFn:=StringPartition;

#P Toggles te way in which Specht prints partition labels.
SpechtPrettyPrint:=function(arg)
  if Length(arg)=0 then
    if SpechtPrintFn=LabelPartition then
      SpechtPrintFn:=StringPartition;
    else SpechtPrintFn:=LabelPartition;
    fi;
  elif Length(arg)=1 and IsBool(arg[1]) then
    if arg[1] then SpechtPrintFn:=LabelPartition;
    else SpechtPrintFn:=StringPartition;
    fi;
  else Error("usage: SpechtPrettyPrint() or SpechtPrettyPrint(<bool>)\n");
  fi;
end;
