
function RSS_ALL=files(name)
Files1 = dir(fullfile(name));
LengthFiles1 = length(Files1);
RSS_ALL=cell(1,LengthFiles1);

for ii=1:LengthFiles1
    name=Files1(ii).name;
    folder=Files1(ii).folder;
    all=[folder,'\',name];
    RSS=get_rss(all);
    RSS_ALL{ii}=RSS;
end

end


