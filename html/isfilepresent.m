function isfilepresent(filename)
  
if exist(filename,'file')
  doc(filename);
else
  switch filename
    case 'fsolve'
      web();
    case 'pathmc'
      recsdirectory = strrep(which('recsSimul'),'recsSimul.m','');
      web(['file://' recsdirectory 'html/pathnotinstalled.html'],'-helpbrowser');
  end
end  