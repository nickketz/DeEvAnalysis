function ana = deevFixTrlOrd(ana)
% function to make sure condn is first column of ana.trl_order
%

fn = fieldnames(ana.trl_order);
fprintf('\nChecking trl_order for %d conditions...\n',length(fn));
for icond = 1:length(fn)
    if ~strcmp(ana.trl_order.(fn{icond}){1},'condn')
        fprintf('Fixing %s...',fn{icond})
        ana.trl_order.(fn{icond}) = cat(2,'condn',ana.trl_order.(fn{icond}));
        fprintf('now has %d columns\n',length(ana.trl_order.(fn{icond})));
    end
end