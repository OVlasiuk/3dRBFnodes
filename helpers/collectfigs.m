%COLLECTFIGS
% Collects all the visible figures and saves them as pdf's, with names
% consisting of the name of the calling script + a number.
d = dbstack;
if length(d) > 1
    name = d(2).name;
else
    name = 'console';
end
nFigs = length(get(0,'Children'));

disp('Saving the following figures:')
for thisFigNum=1:nFigs
    figName = [name '_' int2str(thisFigNum) '.pdf'];
    disp(['   ' figName])
    print(gcf, figName,'-dpdf','-r300','-bestfit')
    close
end