% PLAN: Extract features from a 256^3 sample and run kmeans on the result
% to get a better picture of what kinds of features are being identified.
cd /media/chingd/OCT14C/OCT14C/Research/SURF3D/256
load vol.mat

% Detect and extract features.
[points,scale,~] = detectSURF3D(V);
[descriptors,isvalid] = extractFeatures3D(V,points,scale);

fprintf(1,'\nClean up and save work.');
clear scale
save('descriptors256.mat','descriptors','points','isvalid');

%% Section 2

fprintf(1,'\nGet rid of nonvalid points.');
numvalidpoints = sum(isvalid);
sortme = [isvalid,points,descriptors];
sortme = sortrows(sortme,-1);
sortme = sortme(1:numvalidpoints,:);
descriptors = sortme(:,5:end);
points = sortme(:,2:4);
clear sortme isvalid

%% Section 3
fprintf(1,'\nGroup the descriptors into clumps.');
numclusters = round(0.1*numvalidpoints);
[groupings,centroids] = kmeans(double(descriptors),400,'MaxIter',numclusters,'Options',statset('UseParallel',1),'Replicates',5);

sortme = sortrows([groupings,points,descriptors],1);
groupings = sortme(:,1);
points = [sortme(:,2:4),repmat(21,length(groupings),1)];
descriptors = sortme(:,5:end);

%% Section 4

fprintf(1,'\nSave images of the points at each of the clumps.');
lo = int16(1); hi = int16(1); numpoints = length(points);
while lo <= numpoints
    % Select a range [lo,hi) of points in the same group.
    while hi <= numpoints && groupings(lo) == groupings(hi)
        hi = hi + 1;
    end
    % Extract the data.
    grouppoints = points(lo:hi-1,:);
    groupdescriptors = descriptors(lo:hi-1);
    groupoutdir = sprintf('./clusters/group%01i',groupings(lo));
    mkdir(groupoutdir);
    
    markkeypoints3(V,double(grouppoints),groupoutdir);
    save([groupoutdir '/descriptors.mat'],'groupdescriptors');
    
    lo = hi;
end
    
    
    
    
    
    
    
