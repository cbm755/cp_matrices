function m = signmin(a,b)
m = 0*a;

ind = find( abs(a) <= abs(b) );
m(ind) = a(ind);

ind = find( abs(a) > abs(b) );
m(ind) = b(ind);