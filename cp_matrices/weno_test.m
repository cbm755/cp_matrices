%example_heat_sphere
cp.x1d = x1d;
cp.y1d = y1d;
cp.z1d = z1d;
cp.cpx = cpxg;
cp.cpy = cpyg;
cp.cpz = cpzg;
cp.dim = 3;
cp.band = band;

disp('*** weno4 call');
w1 = weno4_interp(cp, u, [cp.cpx cp.cpy cp.cpz]);

disp('*** weno6 call');
w2 = weno6_interp(cp, u, [cp.cpx cp.cpy cp.cpz]);


disp('*** weno4_caching: calling build cache');
WenoCache = weno4_interp_caching(cp, u, [cp.cpx cp.cpy cp.cpz], 'cache');

disp('*** cache built, now interpolating');
w3 = weno4_interp_caching(WenoCache, u);
w4 = weno4_interp_caching(WenoCache, 1.1*u);
w5 = weno4_interp_caching(WenoCache, 1.2*u);

