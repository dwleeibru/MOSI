% automatic parcellation: main code portion of MOSI
% criteria are in line 72 to avoid being replaced while loading variables
% this code is for internal use by NCI, requiring some modification for different research need
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step I; retrieve time series for each ROI
% focus on cortex
%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
sub=dir('*.results'); % subj dirs

% these regional indices are from the FreeSurfer outputs
reg = {'frontal_ind','temporal_ind','limbic_ind','parietal_ind','occipital_ind','subcortical_ind','other_ind'};
lion = {'Fro', 'Tem', 'Lim', 'Par', 'Occ', 'Oth'};
lion_no = [9 6 10 5 3 2]; % right side subregion number in each lobe/system
% in each subj folder Fro_1 to Fro_18, first 9 for right, second 9 for left
%side='right';
side = 'left';

addpath('~/BrainConnectivityToolBox/2016_01_16_BCT')

for i=1:length(sub) 
    cd(sub(i).name)
    %mkdir AutoFP
    cd AutoFP
    delete(['parti_*' side '.mat']) % delete saved results, starting from new gammas
    
    for gamma= 0.45:0.05:0.95 % modularity gamma value
        tmp=[];tmp1=[];
        tmp=dir(['parti_*' side '.mat']); % tmp is not empty since gamma is 0.45
        
        if ~isempty(tmp)
            for j=1:length(tmp)
                tmp1=[tmp1 str2num(tmp(j).name(7:end-5-length(side)))]; % tmp1 is the gamma number
            end
            tmp=find(tmp1<gamma);
        end
        
        badger = cell(1);
        if ~isempty(tmp)
            load(['parti_' num2str(tmp1(tmp(end))) '_' side '.mat']) % load the one prior to current gamma (the MOSI is iterated based on the resutls from previous gamma)
            clear key key_count badger_1 corr_cri dist_cri voxel_cri mu_cri % badger=cell(1) is replaced by badger, Vin and Min are retained
        else
            % to store all the partitions in cells if no parti*.mat - new analysis from lowest gamma
            dim=[]; fid=[];
            fid = fopen ('dim.txt','r');
            dim=textscan(fid,'%f'); dim=dim{1}; % dimension of the brain
            fclose(fid);
            
            % perform 1st partition
            for j=1:length(lion)
                if strcmp(side,'right')==1
                    sss=1:lion_no(j);
                elseif strcmp(side,'left')==1
                    sss=lion_no(j)+1:lion_no(j)*2;
                end
                for k=sss     % every roi in lion (broad brain area)
                    fid=[]; data=[]; cc=[]; corr_c=[];
                    fid=fopen([lion{j} '_' num2str(k) '.txt'],'r');
                    cc=textscan(fid,'%f'); % ascii
                    fclose(fid);
                    data=cc{1}; data = reshape (data, [258 length(data)/258])'; % ikj + 255 = 258
                    
                    if isempty(badger{1})
                        badger{1}=data;
                    else
                        badger{end+1}=data;
                    end
                    
                end
            end
        end
        corr_cri = 0.5; % criteria for average corr
        %corr_cri2 = 0.7;
        dist_cri = 0.05; % criteria for average distance
        voxel_cri=10; % least voxel size for a roi
        mu_cri = 0.95; % criteria for mutual information
        
        key = 1; % to determine whether the loop needs to be terminated
        key_count = 0; % to calculate the times of loop
        % perform repeated partition and union
        % check whether key needs to be updated to 1
        last_badger=[]; parti = [];
        tic
        while key
            last_badger=badger;
            for j=1:length(last_badger)
                corr_c=[];
                corr_c=corr(last_badger{j}(:,4:end)');
                % by community detection
                n  = size(corr_c,1);        % number of nodes
                M  = 1:n;                   % initial community affiliations
                Q0 = -1; Q1 = 0;            % initialize modularity values
                while Q1-Q0>1e-5           % while modularity increases; adjut to 1e-4 from 1e-5
                    Q0 = Q1;                % perform community detection
                    [M, Q1]=community_louvain(corr_c,gamma,M,'negative_asym');
                end
                
                mo = [];
                mo = unique(M); % module index
                
                % next loop deals with the condition that the retrived module has disconnected sub-modules
                % the content of badger will be updated
                tp_m1=sub2ind(dim,last_badger{j}(:,1),last_badger{j}(:,2),last_badger{j}(:,3)); % linear position index of the whole original module
                for m=1:length(mo) % investigating the partition
                    tmp=[];
                    tmp=find(M==mo(m));
                    
                    % check whether there are two separate modules in a partition,
                    % using function contig I write
                    tp_m=[];
                    tp_m = contig(dim',last_badger{j}(tmp,:)); % tp_m is the linear index of ijk for disconnected sub-modules
                    
                    tj_ind=0; tp_data=cell(1);
                    if tp_m.NumObjects > 1 % divide the module again, diconnected modules will be separated
                        % convert subscript of data(tmp,:) to linear index
                        for n=1:tp_m.NumObjects
                            tj_ind=tj_ind+1;
                            tp_n=[];tp_n1=[];tp_n2=[];
                            [tp_n,tp_n1,tp_n2]=intersect(tp_m1,tp_m.PixelIdxList{n});
                            tp_data{tj_ind}=last_badger{j}(tp_n1,:);
                        end
                    else % modules are neighbors
                        tj_ind=tj_ind+1;
                        tp_data{tj_ind}=last_badger{j}(tmp,:);
                    end
                    
                    for n=1:length(tp_data) % fill in badger
                        badger{end+1}=tp_data{n};
                    end
                end
                badger{j}=[];
            end
            
            % delete empty cells
            tp_j3=[];
            tp_j3=~cellfun('isempty',badger);
            badger=badger(tp_j3);
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % handle the situation less than voxel_cri voxels as a module
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % repeat at least length(ov_ind) times
            ov_ind=[]; % small roi index (ov: one voxel), referring to badger
            for j=1:length(badger)
                if size(badger{j},1) < voxel_cri
                    ov_ind=[ov_ind j];
                end
            end
            
            % FIRST STEP: determine whether it is possible to unite small rois in ov_ind
            if length(ov_ind)>1 % if length(ov_ind)==1, then nothing can be combined
                tp_nei=[];
                for j=2:length(ov_ind) % filling in lower trianglar of dist_met
                    for k=1:j-1
                        tk=[];tk_1=[];tk_2=[];
                        tk_1=zeros(dim');
                        tk=sub2ind(dim,[badger{ov_ind(k)}(:,1);badger{ov_ind(j)}(:,1)],[badger{ov_ind(k)}(:,2);badger{ov_ind(j)}(:,2)],[badger{ov_ind(k)}(:,3);badger{ov_ind(j)}(:,3)]);
                        tk_1(tk)=1;
                        tk_2=bwconncomp(tk_1);
                        if tk_2.NumObjects == 1
                            tp_nei=[tp_nei; ov_ind(j) ov_ind(k)]; % the neighboring pairs of the clusters with size less than voxel_cri
                        end
                    end
                end
                
                key1=1; % use a while loop to handle the step 1 combination
               
                % combine the neighbors with highest correlation, instead
                % of combination by chance
                if ~isempty(tp_nei)
                    tp_nei=[tp_nei zeros(size(tp_nei,1),1)];
                    for j=1:size(tp_nei,1)
                        tk_1=[]; tk_2=[];
                        tk_1=mean(badger{tp_nei(j,1)}(:,4:end),1);
                        tk_2=mean(badger{tp_nei(j,2)}(:,4:end),1);
                        tp_nei(j,3)=corr(tk_1',tk_2');
                    end
                    
                    tp_nei=sortrows(tp_nei,3,'descend');
                    tp_j1=[]; tp_j1=find(tp_nei(:,3)>corr_cri);
                    tp_nei=tp_nei(tp_j1,:); % only corr > corr_cri are maintained
                else
                    key1=0;
                end                
              
                while key1
                    if isempty(tp_nei)
                        key1=0;
                    else
                        tp_j1=[]; tp_j2=[]; tp_j3=[]; tp_j4=[];
                        tp_j1=find(tp_nei(:,1)==tp_nei(1,1));
                        tp_j2=find(tp_nei(:,2)==tp_nei(1,1));
                        tp_j1=union(tp_j1,tp_j2);
                        tp_j3=find(tp_nei(:,1)==tp_nei(1,2));
                        tp_j4=find(tp_nei(:,2)==tp_nei(1,2));
                        tp_j3=union(tp_j3,tp_j4);
                        tp_j1=union(tp_j1,tp_j3);
                        
                        badger{tp_nei(1,1)}=[badger{tp_nei(1,1)}; badger{tp_nei(1,2)}];
                        badger{tp_nei(1,2)}=[];
                        
                        tp_nei(tp_j1,:)=[]; % remove the relevant rows
                    end
                end
            end
            
            % delete empty cells
            tp_j3=[];
            tp_j3=~cellfun('isempty', badger);
            badger=badger(tp_j3);
            
            % SECOND STEP: the relaitonship between ov_ind and others
            key1=1; 
            while key1                
                ov_ind=[]; % update small roi index (ov: one voxel), referring to badger
                tmp_badger=[]; tmp_badger=badger;
                for j=1:length(badger)
                    if size(badger{j},1) < voxel_cri
                        ov_ind=[ov_ind j];
                    end
                end                

                if ~isempty(ov_ind)
                    % note: the neighboring indices for dim ijk at voxel x are: x+1,
                    % x-1, x+i*j, x-i*j, x+i, x-i but i give up this method, i still use bwconncomp to keep consistency
                    
                    kk_count=0; % counting isolated module that are detached from other modules
                    
                    for j=1:length(ov_ind) % j referred to content of ov_ind
                        tp_j1=[];
                        tp_j1=tmp_badger{ov_ind(j)}(:,4:end);
                        tp_j1=mean(tp_j1,1);
                        
                        tp_nei=[]; % storing neighbors of ov_ind{j}
                        honey=cell(1); honey_ind=0;  % storing the structures of neighboring roi
                        bee=cell(1); % storing the index/position of neighboring roi
                        
                        pollen=[]; % storing the index of ov_ind
                        pollen=sub2ind(dim,tmp_badger{ov_ind(j)}(:,1),tmp_badger{ov_ind(j)}(:,2),tmp_badger{ov_ind(j)}(:,3));
                        
                        for m=setdiff(1:length(tmp_badger),ov_ind) 
                            tk=[];tk_1=[];tk_2=[];tk_3=[];
                            tk_1=zeros(dim');
                            tk_3=sub2ind(dim,tmp_badger{m}(:,1),tmp_badger{m}(:,2),tmp_badger{m}(:,3));
                            tk=[tk_3; pollen];
                            tk_1(tk)=1;
                            tk_2=bwconncomp(tk_1);
                            
                            tk_1(pollen)=0; % remove ov_ind, will be used in the next if end
                            
                            if tk_2.NumObjects == 1
                                tp_nei=[tp_nei m];
                                honey_ind=honey_ind+1;
                                honey{honey_ind}=tk_1; % store the matrix of neighboring roi
                                bee{honey_ind}=tk_3; % store the position index of neighbors
                            end
                        end
                        
                        if ~isempty(tp_nei)
                            corr_nei=[]; % using correlation to determine which cluster should the small voxle belong to
                            for k=1:length(tp_nei)
                                tp_k=[];
                                tp_k=corr(tp_j1',tmp_badger{tp_nei(k)}(:,4:end)');
                                if size(tp_k,2)>voxel_cri
                                    tp_k1=[]; tp_k3=[];
                                    [tp_k3, tp_k1]=bwdist(honey{k}); % index of nearest neighbor
                                    tp_k3=tp_k1(pollen); % the nearest voxels in neighboring roi to ov_ind, the content is vector not ijk
                                    tp_k3=unique(tp_k3);
                                    [tp_n,tp_n1,tp_n2]=intersect(bee{k},tp_k3); % provindg index of tp_ke in bee/tp_k
                                    tp_k2=[]; tp_k2=tp_k(tp_n1); % for larger module, only select the highest voxel_cri correlated voxels
                                    corr_nei=[corr_nei mean(tp_k2)];
                                else
                                    corr_nei=[corr_nei mean(tp_k)];
                                end
                            end                            
                            
                            tp_j2=[]; tp_j3=[];
                            [tp_j2, tp_j3]=max(corr_nei);
                            
                            %if tp_j2 > corr_cri2 % correlation criteria, may delete the loop
                            badger{tp_nei(tp_j3)}=[badger{tp_nei(tp_j3)}; badger{ov_ind(j)}];
                            badger{ov_ind(j)}=[];
                            %end  
                        else
                            kk_count=kk_count+1;
                        end
                    end
                    
                    % delete empty cells
                    tp_j3=[];
                    tp_j3=~cellfun('isempty', badger);
                    badger=badger(tp_j3);
                    
                    if length(ov_ind)==kk_count % isolated modules are the only small modules
                        break
                    end
                    
                else
                    key1=0;
                end             
                                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % call a function to store the partitions
            % the next section is for uniting similar modules into 1 module
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            parti = conn_cells(badger,dim'); % this badger has been updated
            
            % determine whether some modules need to be united into a bigger one
            % by comparing the value of averge corr with the dist_cri
            tp_i1=[]; tp_i2=[];
            [tp_i1, tp_i2]=find(parti==1); % connected pairs
            tp_nei=[tp_i1 tp_i2 zeros(length(tp_i1),1)]; % 3rd column will store the distance between 2 correlation
            
            % preparing relative distance calculation
            raw_data=[]; % update raw_data according to last_badger; consider to add subcortical information as needed!!!
            w_data=[]; % the "weight" of last_badger
            for j=1:length(badger)
                w_data=[w_data size(badger{j},1)];
                raw_data=[raw_data mean(badger{j}(:,4:end),1)']; % small roi condition has been handled, so mean() is ok, should be no error
            end
            
            for j=1:length(tp_i1)
                tp_j1=[]; tp_j2=[];
                tp_j1=corr(mean(badger{tp_i1(j)}(:,4:end),1)',raw_data);
                tp_j2=corr(mean(badger{tp_i2(j)}(:,4:end),1)',raw_data);
                %             tp_j1=atanh(tp_j1); % fisher's transformation
                %             tp_j2=atanh(tp_j2);
                tp_nei(j,3)=0.5*sqrt(sumsqr((tp_j1-tp_j2).*(w_data/sum(w_data)))); % normalization; 0.5 is bcz max 2 min 0
            end
            
            if ~isempty(tp_nei)
                tp_nei=sortrows(tp_nei,3,'ascend'); % smallest at the top
                tp_j1=[]; tp_j1=find(tp_nei(:,3)< dist_cri);
                tp_nei=tp_nei(tp_j1,:);
            end
            
            key1=1; % use a while loop to handle the step 1 combination
            while key1
                if isempty(tp_nei)
                    key1=0;
                else
                    tp_j1=[]; tp_j2=[]; tp_j3=[]; tp_j4=[];
                    tp_j1=find(tp_nei(:,1)==tp_nei(1,1));
                    tp_j2=find(tp_nei(:,2)==tp_nei(1,1));
                    tp_j1=union(tp_j1,tp_j2);
                    tp_j3=find(tp_nei(:,1)==tp_nei(1,2));
                    tp_j4=find(tp_nei(:,2)==tp_nei(1,2));
                    tp_j3=union(tp_j3,tp_j4);
                    tp_j1=union(tp_j1,tp_j3);
                    
                    badger{tp_nei(1,1)}=[badger{tp_nei(1,1)}; badger{tp_nei(1,2)}];
                    badger{tp_nei(1,2)}=[];
                    
                    tp_nei(tp_j1,:)=[]; % remove the relevant rows
                end
            end
            
            tp_j3=[];
            tp_j3=~cellfun('isempty',badger);
            badger=badger(tp_j3);
            
            key_count=key_count+1 % after a splitting and unifying cycle
            kiki=[num2str(length(badger)) '/' num2str(length(last_badger))]
            
            % calculate the similarity of partition of badger and last_badger if key_count > 50
            %
            badger_1=[]; badger_2=[]; %tmp_badger=[]; % storing the voxel number of badger{j}
            for j=1:length(badger)
                tp_m1=[];
                tp_m1=sub2ind(dim,badger{j}(:,1),badger{j}(:,2),badger{j}(:,3));
                badger_1=[badger_1; tp_m1 ones(length(tp_m1),1)*j];
                %tmp_badger=[tmp_badger size(badger{j},1)];
            end
            badger_1=sortrows(badger_1,1); % sorting according to first column (position vector)
            for j=1:length(last_badger)
                tp_m1=[];
                tp_m1=sub2ind(dim,last_badger{j}(:,1),last_badger{j}(:,2),last_badger{j}(:,3));
                badger_2=[badger_2; tp_m1 ones(length(tp_m1),1)*j];
            end
            badger_2=sortrows(badger_2,1); % so, the position of badger_1 and badger_2 are aligned
            
            Vin=[]; Min=[];
            [Vin, Min]=partition_distance(badger_1(:,2),badger_2(:,2))
            if Min > mu_cri && abs(length(badger)-length(last_badger))<= ceil(0.01*length(badger))
                if key_count>50 || Min>0.99
                    key=0;
                end
            end
            %end
            
            if isequal(badger,last_badger)
                key=0;
            else
                badger=badger(randperm(numel(badger))); % permutation of badger order
            end
            
            % exclude the condition that badger has a module with less than
            % 5 voxels but badger = last_badger or Min < mu_cri
            %             tmp_badger=find(tmp_badger<5);
            %             if ~isempty(tmp_badger)
            %                 key=1;
            %             end
            
        end
        time=toc;
        eval(['save parti_' num2str(gamma) '_' side '.mat badger badger_1 dim Vin Min corr_cri dist_cri voxel_cri mu_cri key_count time'])
    end
    cd ..
    cd ..
end

% ===============================================
% FUNCTIONs
% ===============================================

% from partitions to calculate connectivity profile
function a = conn_cells(b,c)
% b is a cell
% c is the dimension
a=zeros(length(b)); % to store connectivity (1 yes, 0 no); b is badger
for i=2:length(b) % lower triangle
    for j=1:i-1
        tmp1=[]; tmp2=[]; tmp3=[];
        tmp1=sub2ind(c,b{i}(:,1),b{i}(:,2),b{i}(:,3));
        tmp2=sub2ind(c,b{j}(:,1),b{j}(:,2),b{j}(:,3));
        
        tmp3=[tmp1;tmp2];
        tmp2=zeros(c);
        tmp2(tmp3)=1;
        tmp1=bwconncomp(tmp2);
        
        if tmp1.NumObjects == 1
            a(i,j)=1;
        else
        end
    end
end
end


% check contiguity
function a = contig(c,d)
% c is the dimension
% d is the data matrix
% use two functions: sub2ind & bwconncomp
tmp=[]; tmp2=[];
tmp=sub2ind(c,d(:,1),d(:,2),d(:,3));
tmp2=zeros(c);
tmp2(tmp)=1; % binary image
a=bwconncomp(tmp2);
end