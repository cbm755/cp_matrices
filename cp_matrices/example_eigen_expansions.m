% run example_two_bands_enc first

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
  plot(thp, evplot, 'color', rand(3,1));

  if (1==0)
    figure(2);
    %plot2d_compdomain(eigenvec-R*(E*eigenvec), x, y, dx, dy, 2);
    plot2d_compdomain(eigenvec, x, y, dx, dy, 2);
    hold on;
    plot(xp, yp, 'k-', 'linewidth', 2);
    title(['eigenvalue = ' num2str(lambda)]);
  end
  pause(0);
end
plottime = toc


% some target function, in embedding space (u) and on circle (up)
u = cos(2*thg); % + cos(3*thg);
up = cos(2*thp); % + cos(3*thp);

% weights are inner product of eigenvectors and u
rawweights = V' * u;

% Alterative: CPext of eigenvectors before computing weights
%rawweights = zeros(size(V,2), 1);
%for i=1:length(rawweights)
%  rawweights(i) = (R*(E*V(:,i)))' * u;
%end


weights = rawweights .* (abs(rawweights) >= .25);

[Lambda weights rawweights]

u2 = zeros(size(thg));
% add up a sum of eigenvectors
for i=1:8
  u2 = u2 + weights(i) * V(:,i);
end

%usewind = [1 2 3 5];
usewind = [3];
u3 = zeros(size(thg));
for i=usewind
  u3 = u3 + weights(i)*V(:,i);
end
%u3 = weights(2)*V(:,2) + ...
%     weights(3)*V(:,3);

figure(3); clf;
plot(thp, up, 'r-');
hold on;
xlabel('t'); ylabel('u');
plot(thp, Eplot*u2, 'b--');

plot(thp, Eplot*u3, 'c--')

legend('target', 'reconstr', 'manual reconstr');

figure(4); clf;
xlabel('t'); ylabel('u');
plot(thp, up - Eplot*u3, 'c--')