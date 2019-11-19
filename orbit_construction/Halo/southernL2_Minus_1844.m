% Southern L2-

% Grabs the southernn equivalent of the manifolds for the northern orbits

initCondsNorthern = [];
initCondsSouthern = [];
endConds = [];
 
for i = 1:6:11064
    
   initCondsNorthern = [initCondsL2(i:i+5)'; initCondsNorthern];
      
end

initCondsSouthern = initCondsNorthern;

for i = 1:length(initCondsNorthern)
    initCondsSouthern(i,3) = -initCondsSouthern(i,3);
end

for i = 202:length(initCondsSouthern)
   
    tf = initTimeL2(i);
    y2 = manifold_in(initCondsSouthern(i,:)',tf);
    endConds = [endConds; y2(:,:)];
    fprintf('Finished processing orbit #%d\n', i);
    
end

fprintf('Done!');