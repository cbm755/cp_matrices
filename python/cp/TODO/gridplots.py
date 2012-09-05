"""
plot a cp grid, this should live in cpGrid.py (or some other library).
Should be able to ask a TreeGrid object to draw itself!
"""

from numpy import array as a

# which figure to make
#fig = 5

close('all')

if (1==1): # cos curve
    fig_height_pt = 370./2  # Get this from LaTeX using \showthe\textwidth
    inches_per_pt = 1.0/72.27
    golden_mean = (sqrt(5)-1.0)/2.0
    ratio = golden_mean
    fig_width = fig_width_pt*inches_per_pt
    fig_height =fig_width*ratio
    figsize = [fig_width, fig_height]

    hshift = 0.14; vshift = 0.17


if (1==0): # egg curve
    fig_height_pt = 160.  # Get this from LaTeX using \showthe\textwidth
    inches_per_pt = 1.0/72.27
    golden_mean = (sqrt(5)-1.0)/2.0
    #ratio = golden_mean
    ratio = 1.1
    fig_height =fig_height_pt*inches_per_pt
    fig_width = fig_height/ratio
    figsize = [fig_width, fig_height]
    hshift = 0.14; vshift = 0.13

print figsize
smallax = [hshift, vshift,0.97-hshift,0.97-vshift]


fac = -2
params = {#'backend': 'ps',
          'axes.labelsize': 10+fac,
          'text.fontsize': 10+fac,
          'legend.fontsize': 9+fac,
          'xtick.labelsize': 9+fac,
          'ytick.labelsize': 9+fac,
          #'font.size': 9,
          'font.family': 'serif',
          'font.serif': 'Computer Modern Roman',
          'font.sans-serif': 'Computer Modern Sans serif',
          'font.monospace': 'Computer Modern Typewriter',
          'text.usetex': True,
          'ps.usedistiller': 'xpdf',
          'lines.linewidth' : 0.6,  # points
          'lines.markersize' : 1,   # points
          #'lines.markeredgewidth' : 0.0,
          'axes.linewidth' : 0.3,  # in points
          #legend.handlelength'  : 2.0,   # length of lines
          #'legend.handletextsep' : -20.0   # unknown
          'legend.borderpad': 0.25,      #the fractional whitespace inside the legend border
          'legend.labelspacing':  0.0,  #the vertical space between the legend entries
          'legend.handlelength':  1.8,  #the length of the legend handles
          'legend.handletextpad': 0.1,  #the pad between the legend handle and text
          'legend.borderaxespad': 0.75,  #the pad between the axes and legend border
          'legend.columnspacing': 0.7,  #the spacing between columns
          'figure.figsize': figsize
}

#pylab.rcParams.update(params)
rcParams.update(params)

X,Y = q.ParamGrid()

figure(1)
clf()
ax = axes(smallax) #,frameon=False)

plot(X,Y,'b-')

X = []; Y = [];
#for n in Lev:
#    X.append(n.gridpt[0])
#    Y.append(n.gridpt[1])

Xg = []; Yg = [];
Xt = []; Yt = [];
for n in Lex:
    if n.isEvolvePt and n.whichBdy==0:
        X.append(n.gridpt[0])
        Y.append(n.gridpt[1])
    elif n.isEvolvePt and n.whichBdy != 0:
        Xg.append(n.gridpt[0])
        Yg.append(n.gridpt[1])
    elif n.isExtendPt:
        Xt.append(n.gridpt[0])
        Yt.append(n.gridpt[1])

plot(X,Y,'k.');
plot(Xg,Yg,'rx',mew=0.3,markersize=1.2);


xlabel(r'$x$')
ylabel(r'$y$')

axis('scaled')
#axis('equal')

if (1==1):  # cos curve
    xlim(-0.3, 4.3)
    ylim(-1.8, 1.6)
    xticks([0,1,2,3,4])
    yticks([-1,0,1])
if (1==0):  # egg curve
    xlim(-1.45, 1.2)
    ylim(-1.7, 1.6)
    xticks([-1,0,1]);  yticks([-1,0,1])



ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

#show()

savefig('grid1.eps')
