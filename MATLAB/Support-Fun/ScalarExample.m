% SCALAREXAMPLE One performance function and its actual evaluation counter.
% Uses the same mathematical models as Ex1 and Ex2.
% Li, H.-S.; Ma, Y.-Z.; Cao, Z. (2015), Computers and Structures 153,
% 239-251. DOI: 10.1016/j.compstruc.2014.10.014.
classdef ScalarExample < handle
 properties
  Nfun=1;Ndim=2;Rdis;example;pf;count=0;
 end
 methods
  function obj=ScalarExample(ex,j)
   obj.example=ex;obj.pf=j;
   obj.Rdis=Direct([makedist('Normal','mu',0,'sigma',1);makedist('Normal','mu',0,'sigma',1)]);
  end
  function g=EvlLSF(obj,x)
   obj.count=obj.count+size(x,2);
   if obj.example==1
    if obj.pf==1,g=x(1,:).^2+x(2,:).^2-.01;
    else,g=(x(1,:)-.4).^2+(x(2,:)-.4).^2-.01;end
   else
    switch obj.pf
     case 1,g=3+(x(1,:)-x(2,:)).^2/10-(x(1,:)+x(2,:))/sqrt(2);
     case 2,g=x(1,:)-x(2,:)+7/sqrt(2);
     case 3,g=3+(x(1,:)-x(2,:)).^2/10+(x(1,:)+x(2,:))/sqrt(2);
     case 4,g=x(2,:)-x(1,:)+7/sqrt(2);
    end
   end
  end
 end
end
