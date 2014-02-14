function [M, L, E, F] = app_bnd(M, L, E, F, XCP, YCP, BDYG, pt, s, dx_coarsest)
warning('can only deal with boundary conditions of unit circle or semicircle now!')
% in fact, we only need to deal with 'M' on the coarsest grid
% but for reason of debug, we apply boundary conditions to all
% levels of grids
tol = 1e-10;
n_level = length(L);
d = 2;
if nargin < 10
        dx_coarsest = 0.1;
end
if strcmp(s, 'dirichlet') || strcmp(s, 'Dirichlet')
	disp('dealing with Dirichlet boundary condition');
    
	for i = 1:1:n_level
		BDYG{i} = logical(BDYG{i});
        %F{i}(BDYG{i},:) = -F{i}(BDYG{i},:);	    
        E{i}(BDYG{i},:) = -E{i}(BDYG{i},:);
	    M{i} = lapsharp(L{i}, E{i});
	end
elseif strcmp(s, 'neumann') || strcmp(s, 'Neumann')
 	disp('dealing with Neumann boundary condition')
%     for i = 1:1:n_level
%         L{i}(1,:) = sum(E{i});
%         M{i}(1,:) = sum(E{i});
%     end
    for i = 1:1:n_level
        dx = dx_coarsest/2^(n_level-i);
        for j = 1:1:length(XCP{i})
            if abs(XCP{i}(j)-pt(1)) < tol && abs(YCP{i}(j)-pt(2)) < tol
%                 L{i}(j,:) = ones(1,length(L{i}));
%                 M{i}(j,:) = ones(1,length(M{i}));
               L{i}(j,:) = sum(E{i});
               M{i}(j,:) = sum(E{i});
                disp(['changing the rows...  level:', num2str(i)]);
%                  L{i}(j,:) = 0;
%                 %L{i}(:,j) = 0;
%                  L{i}(j,j) = 1;
%                  M{i}(j,:) = 0;
%                 %M{i}(:,j) = 0;
%                  M{i}(j,j) = 1;
                 F{i}(j) = 0;
            end
        end
    end
else 
	error('can only deal with Dirichlet and Nuemann boundary conditions');
end

end
            
