% Generate Matlab Code for given action matrix and coefficient matrices
% (GBsolver subroutine)
% by Martin Bujnak, mar2008


function [res] = gbs_ExportMCode(filename, M, Mcoefs, coefscode, known, knowngroups, unknown, algB, actMvar, amrows, amcols, gjcols, aidx)

    [p probname e] = fileparts(filename);
    if isempty(e)
        filename = [filename '.m'];
    end;
    
    if (~isdir(p))
        mkdir(p);
    end
    
    if isempty(knowngroups)
        knowngroups = 1:length(known);
    end
    
    % generate coefs calculation code
    knvars=[];
    knvarnames=[];
    knvarcnt=[];
    for i=1:length(known)

        if length(knvars) >= knowngroups(i)
            knvars{knowngroups(i)} = [knvars{knowngroups(i)} sym(known{i})];
            knvarcnt(knowngroups(i)) = knvarcnt(knowngroups(i)) + 1;
        else
            knvars{knowngroups(i)} = sym(known{i});
            knvarcnt(knowngroups(i)) = 1;
        end
        
        if length(knvarnames) < knowngroups(i) || isempty(knvarnames(knowngroups(i))) 
            
            %name=regexp(known{i},'\D*', 'match');
            name=known(i);
            knvarnames{knowngroups(i)} = name{1};
        end
    end
    
    fid = fopen(filename, 'w');
    
    fprintf(fid, '%% Generated using GBSolver generator Copyright Martin Bujnak,\n'); 
    fprintf(fid, '%% Zuzana Kukelova, Tomas Pajdla CTU Prague 2008.\n%% \n');
    fprintf(fid, '%% Please refer to the following paper, when using this code :\n'); 
    fprintf(fid, '%%     Kukelova Z., Bujnak M., Pajdla T., Automatic Generator of Minimal Problem Solvers,\n');
    fprintf(fid, '%%     ECCV 2008, Marseille, France, October 12-18, 2008\n');
    fprintf(fid, '\n');
    fprintf(fid, ['function [' c2s((unknown), ' ') '] = ' probname '(' c2s(knvarnames, ', ') ')\n\n']);
    
    % coeffs
    fprintf(fid, '\t%% precalculate polynomial equations coefficients\n');
    for i=1:length(coefscode)
        
        % rename coefficients according to "knowngroups"
        coefcode = char(coefscode(i));
        
        for j=1:length(knvars)
            if length(knvars{j}) > 1
                for k=1:length(knvars{j})
                    coefcode = strrep(coefcode, char(knvars{j}(k)), [knvarnames{j} '(' int2str(k) ')']);
                end
            end
        end
        
        fprintf(fid, ['\tc(' int2str(i) ') = ' coefcode ';\n']);
    end
    fprintf(fid, '\n');
    
    % coefs matrix
    Mcoefs_gj = Mcoefs(:, gjcols);
    
    fprintf(fid, ['\tM = zeros(' int2str(size(Mcoefs_gj, 1)) ', ' int2str(size(Mcoefs_gj, 2)) ');\n']);
    for i=1:length(coefscode)
    
        [ofs] = find(Mcoefs_gj == i);

        if length(ofs) == 1
            
            fprintf(fid, ['\tM(' int2str(ofs) ') = c(' int2str(i) ');\n']);
        else
            fprintf(fid, ['\tci = [']);
            fprintf(fid, l2s(ofs, ', '));
            fprintf(fid, '];\n');

            fprintf(fid, ['\tM(ci) = c(' int2str(i) ');\n\n']);
        end
    end
    fprintf(fid, '\n');
   
    % elimination part
    fprintf(fid, '\tMr = rref(M);  %% replace me with a MEX\n');
  	%fprintf(fid, '\ttol = 2.2204e-20;\n');
	%fprintf(fid, '\tMr = gaussjord(M, tol);\n');

    fprintf(fid, '\n');
    
    % action matrix
    fprintf(fid, ['\tA = zeros(' int2str(length(amrows)) ');\n']);
    fprintf(fid, ['\tamcols = [' l2s(amcols, ' ') '];\n']);
    for i=1:length(amrows)
        
        if amrows(i) < 0
            fprintf(fid, ['\tA(' int2str(i) ', ' int2str(-amrows(i)) ') = 1;\n']);
        else
            fprintf(fid, ['\tA(' int2str(i) ', :) = -Mr(' int2str(amrows(i)) ', amcols);\n']);
        end
    end
    fprintf(fid, '\n');

    % solution extraction
    
	fprintf(fid, '\t[V D] = eig(A);\n');
    
    [oneidx unksidx] = gbs_GetVariablesIdx(algB, unknown);
    varsinvec = find(unksidx > 0);
    
    fprintf(fid, ['\tsol =  V([' l2s(unksidx(varsinvec), ', ') '],:)./(ones(' int2str(length(varsinvec)) ', ' int2str(oneidx) ')*V(' int2str(oneidx) ',:));\n']);
    fprintf(fid, '\n');

	fprintf(fid, '\tif (find(isnan(sol(:))) > 0)\n\t\t\n');

    % division by zero filter
    ucnt = length(unknown);
    for i=1:ucnt
        
        fprintf(fid, ['\t\t' unknown{i} ' = [];\n']);
    end
    
    fprintf(fid, '\telse\n\t\t\n');
    
    % extract variables
    if (sum(unksidx == 0)) > 0

        idx = find(unksidx == 0);
        if (length(idx) > 1)
            
            fprintf(fid, '\t\tWARNING: cannot extract all unknowns at once. A back-substitution required (not implemented/automatized)\n');
        end
        
        fprintf(fid, ['\t\tev  = diag(D);\n']);   
        fprintf(fid, ['\t\tI = find(not(imag( sol(1,:) )) & not(imag( ev )));\n']);
    else
        fprintf(fid, ['\t\tI = find(not(imag( sol(1,:) )));\n']);
    end
    
    ui = 1;
    for i=1:ucnt

        if unksidx(i) == 0
            
            % action variable ?
            if strcmp(actMvar, unknown{i})

                fprintf(fid, ['\t\t' unknown{i} ' = ev(I);\n']);   
            else

                fprintf(fid, '\t\tWARNING: one or more unknowns could not be extracted.\n');
            end
            
        else 
            fprintf(fid, ['\t\t' unknown{i} ' = sol(' int2str(ui) ',I);\n']);
            ui = ui+1;
        end
    end
    
    fprintf(fid, '\tend\n');
    fprintf(fid, 'end\n');
    fclose(fid);
end
