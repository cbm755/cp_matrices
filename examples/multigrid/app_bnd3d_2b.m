function [M, L, E, F] = app_bnd3d_2b(M, L, E, F, a_innerInOuter, XCP, YCP, ZCP, BDYG, pt, s)
warning('can only deal with boundary conditions of unit sphere or semisphere now!')
% in fact, we only need to deal with 'M' on the coarsest grid
% but for reason of debug, we apply boundary conditions to all
% levels of grids
tol = 1e-10;
n_level = length(L);
if strcmp(s, 'dirichlet') || strcmp(s, 'Dirichlet')
	disp('dealing with Dirichlet boundary condition')
	for i = 1:1:n_level
		BDYG{i} = logical(BDYG{i});
		E{i}(BDYG{i},:) = -E{i}(BDYG{i},:);
	    M{i} = lapsharp(L{i}, E{i});
	end
elseif strcmp(s, 'neumann') || strcmp(s, 'Neumann')
	disp('dealing with Neumann boundary condition')
    for i = 1:1:n_level
        for j = 1:1:length(XCP{i})
            if abs(XCP{i}(j)-pt(1)) < tol && abs(YCP{i}(j)-pt(2)) < tol && abs(ZCP{i}(j)-pt(3)) < tol
                disp(['changing the rows...  level:', num2str(i)]);
                L{i}(j,:) = 0;
                %L{i}(:,j) = 0;
                L{i}(j,a_innerInOuter{i}(j)) = 1;
                M{i}(j,:) = 0;
                %M{i}(:,j) = 0;
                M{i}(j,j) = 1;
                F{i}(j) = 0;
            end
        end
    end
else 
	error('can only deal with Dirichlet and Nuemann boundary conditions');
end

end
            