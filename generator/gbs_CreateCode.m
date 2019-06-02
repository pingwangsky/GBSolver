% This code uses Grobner basis methods to generate a code solving system of
% given polynomial equations.
%
% Please refer to the following paper, when using this code : 
%
%   "Kukelova Z., Bujnak M., Pajdla T., Automatic Generator of Minimal
%    Problem Solvers, ECCV 2008, Marseille, France, October 12-18, 2008"
%
%
% by Martin Bujnak, mar2008
%

function [res export] = gbs_CreateCode(codename, eq, known, unknown, kngroups, cfg, algB)

    export = [];
    
    if nargin < 5
       
        kngroups = [];
    end

    if nargin < 6
        
        cfg = gbs_InitConfig();
    end

    if nargin < 7
        
        algB = [];
    end    
    
    % parse & detect & add required polynomials
    [M Mcoefs coefscode amVar amLT amLTall algBidx algB] = gbs_PreparePolySystemCoefsMatrix(cfg, eq, known, unknown, algB);
    
    if ~isempty(M)
        
        % prepare code lookups
        [amrows amcols gjcols aidx] = gbs_BuildActionMatrix(cfg, M, algB, amLT, amLTall, algBidx, amVar);

        % save calculated data
        export.M = M;
        export.Mcoefs = Mcoefs;
        export.coefscode = coefscode;
        export.known = known;
        export.unknown = unknown;
        export.algB = algB;
        export.actMvar = amVar;
        export.amrows = amrows;
        export.amcols = amcols;
        export.gjcols = gjcols;
        export.aidx = aidx;

        if cfg.bGeneratorSolverResult
            save(['out\' codename '_expparams'], 'export');
        end

        for tmp = cfg.exportCode

            % cut suffix
            parls = strfind(tmp{1}, '(');
            if isempty(parls)
                suffix = [];
            else
                suffix = tmp{1}(parls(end)+1:end-1);
                tmp{1} =tmp{1}(1:parls(end)-1);
            end
            
            if strcmpi(tmp{1}, 'matlab')

                fprintf('--- generating matlab solver ---\n');
                gbs_ExportMCode(['out\' codename suffix], M, Mcoefs, coefscode, known, kngroups, unknown, algB, amVar, amrows, amcols, gjcols, aidx);

            elseif strcmpi(tmp{1}, 'maple')

                fprintf('--- generating Maple solver ---\n');
                if isempty(suffix)
                    suffix = '.txt';
                end
                gbs_ExportMapleCode(['out\' codename suffix], M, Mcoefs, coefscode, known, kngroups, unknown, algB, amVar, amrows, amcols, gjcols, aidx);
            else

                try
                    fprintf('--- generating c solver using '' %s '' template---\n', tmp{1});
                    if isempty(suffix)
                        suffix = 'mex.c';
                    end
                    gbs_ExportCCode(['out\' codename suffix], tmp{1}, M, Mcoefs, coefscode, known, kngroups, unknown, algB, amVar, amrows, amcols, gjcols, aidx);
                catch
                    fprintf('... error exporting using ''%s'' template---\n', tmp{1});
                end
            end
        end
    end
    
    res = 0;
end