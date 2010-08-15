#######################################################################
##  SPECHT - dispatcher.gi : Compability for variable numbers of     ##
##           arguments                                               ##
##                                                                   ##
##     This file contains a hack to support a variable number of     ##
##     arguments for operations as it was common in GAP 3            ##
##     For usage in functions that take partitions as arguments      ##
##                                                                   ##
##     These programs, and the enclosed libraries, are distributed   ##
##     under the usual licensing agreements and conditions of GAP.   ##
##                                                                   ##
##     Dmitriy Traytel                                               ##
##     Thanks @Max Neunhoeffer for help                              ##
##                                                                   ##
#######################################################################

## 3.0: June 2010:
##   - initial

######################################################################

## MakeDispatcherFunc := function(name)
##   nameop := Concatenation(name,"Op");
##   DeclareOperation(nameop,[IsHecke,IsList]);
##   oper := ValueGlobal(nameop);
##   disp := function(arg)
##     if Length(arg) = 2 and IsList(arg[2]) then
##       return oper(arg[1],arg[2]);
##     else
##       return oper(arg[1],arg{[2..Length(arg)]});
##     fi;
##   end;
##   BindGlobal(name,disp);
## end;

#F Dispatches an function with variable number of arguments to an operation
## with a fixed number of arguments
## For internal usage only
## name 	- operation name
## filts 	- list of lists of filters
## pos 		- list of positions for variable length arguments
## l			- list of total numbers of arguments
DeclareGlobalFunction("MakeDispatcherFunc");

InstallGlobalFunction(MakeDispatcherFunc,
  function(name,filts,pos,l)
    local nameop, filters, oper, disp, newarg, i, p, ps;

    nameop := Concatenation(name,"Op");
    
    for i in [1..Length(pos)] do
      if pos[i] = 0 then 
        ## Print(nameop,":",filts[i],"\n"); ## DEBUG
        DeclareOperation(nameop,filts[i]);
      else ## pos is assumed to be positive and <= Length(filts)+1;
        filters := Concatenation(filts[i]{[1..pos[i]-1]},
          [IsList],filts[i]{[pos[i]..Length(filts[i])]});
        if Length(filters) <> l[i]
        then Error("usage, number of filters does not correspond",
          "to the given argument list length");
        fi;
        ## Print(">>",nameop,":",filters,"\n"); ## DEBUG
        DeclareOperation(nameop,filters);
      fi;
    od;
    oper := ValueGlobal(nameop);
    disp := function(arg)
      ps := Positions(l,Length(arg));
      for p in ps do
				if p<>fail and (pos[p]=0 or IsList(arg[pos[p]])) then
        	## Print("Do call ",oper,"(",arg,")\n"); ## DEBUG
        	return CallFuncList(oper,arg);
				fi;
			od;
      p := 1; ## first entry ist the default entry
      newarg := Concatenation(arg{[1..pos[p]-1]},
        [arg{[pos[p]..pos[p]+(Length(arg)-l[p])]}],
        arg{[pos[p]+(Length(arg)-l[p])+1..Length(arg)]});
      ## Print(">> Do call ",oper,"(",newarg,")\n"); ## DEBUG
      return CallFuncList(oper,newarg);
    end;
    BindGlobal(name,disp);
  end
);
