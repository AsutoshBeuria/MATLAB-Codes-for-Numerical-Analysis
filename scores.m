% for index=[matrix]
%commands to be executed
%end
   if age<16
    disp('sorry you are too young to apply')
elseif age<18
    disp('you can apply for youth license')
elseif age<70
    disp('you may have a standard driving license')
else 
    disp('drivers over 70 require special license')
   end