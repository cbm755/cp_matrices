% run example_two_bands_enc first
% A modifiction of example_eigen_expansions (Aug 2, 2011)
% - changed the function to approx to one which is discontinuous at t=2
% - modified the exact solution so that the second and thrid coefficients are
% zero after the decimal place (green) or zero after the first decimal place (blue).
% Red gives exact.
% - the second result is near perfect in eyeball norm.  Even the first, crude, approx 
% captures the overall shape very well.

tic
%[V,D] = eigs(-M, 20, 'sm', opts);
[Vuns,D] = eigs(-M, 30, 0.5);
evtime = toc
D = diag(D);
[Lambda,I] = sort(abs(D));
Lambda
V = Vuns(:,I);

tic
figure(1); clf; hold on;
xlabel('x'); ylabel('u');



for i=1:8
  lambda = Lambda(i);
  eigenvec = V(:,i);
  w = eigenvec' * eigenvec

  evplot = Eplot*real(eigenvec);
  evplot = reshape(evplot, size(xp));

  figure(1);
  plot(thetap, evplot, 'color', rand(3,1));

  if (1==0)
    figure(2);
    %plot2d_compdomain(eigenvec-R*(E*eigenvec), x, y, dx, dy, 2);
    plot2d_compdomain(eigenvec, x, y, dx, dy, 2);
    hold on;
    plot(xp, yp, 'k-', 'linewidth', 2);
    title(['eigenvalue = ' num2str(lambda)]);
  end
end
plottime = toc


% some target function, in embedding space (u) and on circle (up)
%u = cos(thg + 2); % + cos(3*thg);
%up = cos(thetap + 2); % + cos(3*thetap);
u = sign(thg - 2).* cos(3*thg);
up = sign(thetap - 2).* cos(3*thetap);

% weights are inner product of eigenvectors and u
rawweights = V' * u;

% Alterative: CPext of eigenvectors before computing weights
%rawweights = zeros(size(V,2), 1);
%for i=1:length(rawweights)
%  rawweights(i) = (R*(E*V(:,i)))' * u;
%end


%weights = rawweights .* (abs(rawweights) >= .1);
weights = rawweights;

[Lambda weights rawweights]

u2 = zeros(size(thg));
% add up a sum of eigenvectors
for i=1:8
  u2 = u2 + weights(i) * V(:,i);
end

%u3 = weights(1)*V(:,1) + weights(2)*V(:,2) + ...
%     weights(3)*V(:,3);

u3 = up + Eplot*(-.9462*V(:,2) + ...
     +.4196*V(:,3));

u4 = up + Eplot*(-.0462*V(:,2) + ...
     +.0196*V(:,3));

figure(3); clf;
plot(thetap, up, 'r-');
hold on;
xlabel('t'); ylabel('u');
plot(thetap, u4, 'b--');

plot(thetap, u3, 'g--')

figure(4); clf;
xlabel('t'); ylabel('u');
%plot(thetap, up - Eplot*u3, 'c--')
