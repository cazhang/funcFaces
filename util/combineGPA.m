
function gpa = combineGPA(V)
%% Combine GPA and Cout
    if isGPA == 1
        mappedV = V;
        for i=1:Nsets
            if i==1
                mappedV(:,:,i) = V(:,:,i);
            else
                %Cin = Xicp{1}(:,:,i-1);
                %Cout = closestRotation(Cin);
                Cin = Xsol{1}(:,:,i-1);
                mappedV(:,:,i) = Cin * V(:,:,i); % C(1,i-1)*V(i)
            end
        end
        MaxSteps = 5;
        BasisModel = struct([]);
        for i=1:Nsets
            BasisModel(1,i).vertices = mappedV(:,:,i)';
        end
        % init Xgpa
        Xgpa = Xsol{1};
        for i=1:MaxSteps
            [R, t, s, Centroid, corr, registeredModel] = globalGenProcrustes(BasisModel, 1);
            % update BasisModel
            BasisModel = registeredModel;
            % update Cout
            for j=1:Nsets-1
                Xgpa(:,:,j) = R{1, 1}' * R{j+1, 1} * Xgpa(:,:,j);
            end
        end
        
        
        %% ICP-GPA method
    elseif isGPA == 0
        mappedV = V;
        for i=1:Nsets
            if i==1
                mappedV(:,:,i) = V(:,:,i);
            else
                %Cin = Xicp{1}(:,:,i-1);
                %Cout = closestRotation(Cin);
                Cin = Xsol{1}(:,:,i-1);
                Cout = Cin;
                mappedV(:,:,i) = Cout * V(:,:,i); % C(1,i-1)*V(i)
            end
        end
        MaxSteps = 5;
        BasisModel = struct([]);
        for i=1:Nsets
            BasisModel(1,i).vertices = mappedV(:,:,i)';
        end
        tic;
        %figure(1); clf;
        %visModel(BasisModel);
        [R, t, s, Centroid, corr, registeredModel] = globalGenProcrustes(BasisModel, MaxSteps);
        toc;
        %figure(2); clf;
        %visModel(registeredModel);
        
        % update solution
        Xgpa = Xsol{1};
        
        %     %% using the accumulated rsu
        %     RR = cell(Nsets, 1);
        %     for i=1:Nsets
        %         for j = 1:MaxSteps
        %             if j==1
        %                 RR{i,1} = R{i,j};
        %             else
        %                 RR{i,1} = R{i,j} * RR{i,1};
        %             end
        %         end
        %     end
        %
        %     for i=1:Nsets-1
        %         Xgpa(:,:,i) = RR{i,1}' * RR{i+1,1} * Xgpa(:,:,i);
        %     end
        
        %% Only use the last result
        for i=1:Nsets-1
            for j=1:MaxSteps
                %     %Xgpa(:,:,i) = Xgpa(:,:,i) * R{i+1, MaxSteps} * R{1, MaxSteps}';
                %     %Xgpa(:,:,i) = R{1, MaxSteps}' * R{i+1, MaxSteps} *Xgpa(:,:,i) ;
                % Xgpa(:,:,i) = R{1, MaxSteps}' * Xgpa(:,:,i) * R{i+1, MaxSteps} ;
                Xgpa(:,:,i) = R{1, j}' * R{i+1, j} * Xgpa(:,:,i);
            end
        end
    end