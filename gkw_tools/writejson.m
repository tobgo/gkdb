% Write a json file from the structure prepared by the routine gkw2json
%
% YC - 21.06.2017

function str=writejson(data,flnm,n_indent)


indent='    ';
nl=char(10);
if ~exist('n_indent')||isempty(indent)
 n_indent=0;
end

if ~exist('flnm')
 flnm=[];
end

field_names=fieldnames(data);
str=[repmat(indent,1,n_indent) '{' nl];
n_indent = n_indent+1;

for ii=1:length(field_names)
 str=[str repmat(indent,1,n_indent) '"' field_names{ii} '": '];
 sub_data=getfield(data,field_names{ii});
 if isstruct(sub_data)
   n_indent=n_indent+1;
   str=[str nl writejson(sub_data,[],n_indent)];
   if ii~=length(field_names), str=[str ',']; end
   str = [str nl];
   n_indent=n_indent-1;
 elseif iscell(sub_data)
   str=[str indent '[' nl];
   n_indent=n_indent+1;
   for jj=1:length(sub_data)
    if isstruct(sub_data{jj})
      str=[str writejson(sub_data{jj},[],n_indent)];
      if jj~=length(sub_data), str=[str ',' nl]; end
    else
      % write string value
      str=[str repmat(indent,1,n_indent) '"' sub_data{jj} '"'];
      if jj~=length(sub_data), str=[str ',' ]; end
      str=[str nl];
    end
   end
   n_indent=n_indent-1;
   str=[str repmat(indent,1,n_indent) ']'];
   if ii~=length(field_names), str=[str ',']; end
   str=[str  nl];
 else
   % write string or array
   if ischar(sub_data)
     str = [str '"' sub_data '"'];
   elseif isnumeric(sub_data)
     S=size(sub_data);
     if prod(S)==1  % no brackets for single elements
       str=[str sprintf('%0.8g',sub_data)];
     else
       str = [str writenumeric(sub_data)];
     end
   elseif islogical(sub_data)
     if sub_data==1  
       str = [str 'true'];
     else
       str = [str 'false'];
     end
   else 
    disp(['Data type not handled for ' field_names{ii}])
    return 
   end
   if ii~=length(field_names), str=[str ',']; end
   str=[str  nl];
 end
end
n_indent = n_indent-1;
str=[str repmat(indent,1,n_indent) '}'];

if ~isempty(flnm)
 str=[str repmat(indent,1,n_indent) nl];
 fid=fopen(flnm,'w');
 fprintf(fid,'%s',str);
 fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

function  str_out = writenumeric(array)


str_out='';
S=size(array);
if length(S)==2 & S(2)==1 
  str_out=['[' sprintf('%0.8g',array(1))];
  for ii=2:S(1)
    str_out=[str_out ',' sprintf('%0.8g',array(ii))];
  end
  str_out=[str_out ']'];
else
  str_out = '['; 
  for ii=1:S(1)
    switch length(S)
    case 2
      str_out=[str_out writenumeric(shiftdim(array(ii,:),1))];
    case 3 
      str_out=[str_out writenumeric(shiftdim(array(ii,:,:),1))];
    case 4
      str_out=[str_out writenumeric(shiftdim(array(ii,:,:,:),1))];
    otherwise
     disp('Arrays with more than 5 dims not implemented')
    end
    if ii~=S(1), str_out=[str_out ',']; end
  end
  str_out = [str_out ']'];
end