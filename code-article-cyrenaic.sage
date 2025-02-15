# Release Date: 2019-01-15
# last modified and tested with SageMath version 10.6.beta5, Date: 2025-02-10
# authors: Thomas Fernique and Carole Porrier

###############
# TERMINOLOGY #
###############

# slope: an nxd matrix whose lines generate a d-dim. plane of R^n
# subperiod of a slope E: a vector of E with d integer entries
# projection: a dxn matrix A which defines a projection from R^n to R^d
# multigrid: a dxn matrix G whose lines define the direction of each grid (normal vector to each hyperplane)
# prototile: the d indices of the vectors of the standard basis which define it
# tile: a pair (t,pos) where t is a prototile and pos the integer translation applied on it
# tiling: a list of tiles
# pattern: list of tiles
# atlas: list of patterns

# faire une classe tiling ?

#################
# Cut & Project #
#################

# DESCRIPTION: compute the multigrid associated with a slope
# INPUT: a slope E
# OUTPUT: a multigrid
def generators_to_grid(E):
    (d,n)=E.dimensions()
    M=matrix(E.rows()+identity_matrix(n).rows())
    G,mu=M.gram_schmidt()
    return matrix([[G[j,i]/norm(G[j]) for j in range(d)] for i in range(n)])

# DESCRIPTION: compute the dual tiling of a multigrid
# INPUT: a multigrid G; an n-dim. shift vector S for the grid; an integer k which bounds the number of lines in each grid
# OUTPUT: the strongly planar tiling of slope E shifted by S
# TODO : ne garder que les tuiles qui ont 4 voisines ?
def dual(G,S,k):
    tiles=[]
    (n,d)=G.dimensions()
    indices=cartesian_product([range(-k,k+1) for i in range(d)])
    for t in Combinations(range(n),d):
        M=G.matrix_from_rows(t).inverse()
        for ind in indices:
            intersec=M*vector([S[t[j]]+ind[j] for j in range(d)])
            pos=vector(ZZ,[ceil(intersec*G[j]-S[j]) for j in range(n)])
            for h in range(d):
                pos[t[h]]=ind[h]
            tiles.append([t,pos])
    return tiles

# DESCRIPTION: check that a given projection is valid
# INPUT: a projection A; a slope E (embedded in AA)
# OUTPUT: boolean (is the projection valid)
def valid_projection(A,E):
    if A.base_ring()!=AA or E.base_ring()!=AA: # both A and E must be embedded in AA
        raise ValueError
    (d,n)=E.dimensions()
    G1=vector([sign(A.matrix_from_columns(t).det()) for t in Combinations(range(n),d)])
    G2=vector([sign(E.matrix_from_columns(t).det()) for t in Combinations(range(n),d)])
    return G1==G2 or G1==-G2

# DESCRIPTION: compute the orthogonal projection on E (always valid)
# INPUT: a slope E
# OUTPUT: a projection
orthogonal_projection=lambda E: E.gram_schmidt()[0]

# DESCRIPTION: compute the projection which maps the e_i's on the roots of the unity
# only for n->2 tilings
def rectified_projection(E):
    (d,n)=E.dimensions()
    if d!=2:
        raise ValueError
    M=matrix(AA,[[cos(k*pi/n) for k in range(2*n)],[sin(k*pi/n) for k in range(2*n)]])
    for p in Combinations(range(2*n),n): # Arrangements nécessaire ?
        A=M.matrix_from_columns(p)
        if valid_projection(A,E):
            return A

##############
# Subperiods #
##############

# DESCRIPTION: compute the integer entries  of the subperiods of a slope
# INPUT: a slope E
# OUTPUT: the list of the subperiods of E
def subperiods(E):
    (d,n)=E.dimensions()
    l=[]
    for t in Combinations(range(n),d+1):
        G=[(-1)^t.index(i)*det(E.matrix_from_columns([j for j in t if j!=i])) for i in t]
        G=matrix(ZZ,[g.list() for g in G])
        for s in G.left_kernel().basis():
            l.append((t,s))
    return l

# DESCRIPTION: determine whether a slope is characterized by its subperiods
# INPUT: a slope E
# OUTPUT: boolean
def is_determined_by_subperiods(E):
    (d,n)=E.dimensions()
    l=subperiods(E)
    Z=PolynomialRing(K,'x',len(l)*(n-d-1))
    M=matrix([[s[t.index(i)] if i in t else Z.gens()[l.index((t,s))*(n-d-1)+i-len([j for j in t if j<i])] for i in range(n)] for (t,s) in l])
    return ideal(Z,M.minors(d+1)).dimension()==0

# DESCRIPTION: subperiods of E lifted in E - E must be determined by its subperiods
# INPUT: a slope E (coordinates in QQ(a))
# OUTPUT: the subperiods of E lifted in E 
def lifted_subperiods(E):
    (d,n)=E.dimensions()
    l=subperiods(E)
    Z=PolynomialRing(K,'x',len(l)*(n-d-1))
    M=matrix([[s[t.index(i)] if i in t else Z.gens()[l.index((t,s))*(n-d-1)+i-len([j for j in t if j<i])] for i in range(n)] for (t,s) in l])
    I=ideal(Z,M.minors(d+1))
    if I.dimension()>0: # E not determined by subperiods
        raise ValueError
    for v in I.variety(K):
        z=M.subs(v).change_ring(K)
        if rank(z.transpose().augment(E.transpose()))==d:
            return z.rows()

# could add the lines of E as vectors of M -> more equations -> should always yield dim(I)=0

##########################
# Atlas of a n->2 tiling #
##########################

# DESCRIPTION: compute the orthogonal projection on the internal space (orthogonal of the slope)
# INPUT: a slope E
# OUTPUT: matrix of the projection
internal_projection=lambda E: matrix(E.right_kernel().basis()).gram_schmidt()[0]

# DESCRIPTION: compute the region in the window associated with a pattern
# INPUT: window W; internal projection ip; the list of vertices of a pattern
# OUTPUT: a polytope
region=lambda W,ip,P: Polyhedron(ieqs=flatten([W.translation(-ip*vector(v)).inequalities_list() for v in P],max_level=1))

# DESCRIPTION: complete the forced tiles of a pattern
# INPUT: window W; angle-sorted list star of the basis vectors projected onto the slope; internal projection ip; a pattern P; the region R of P
# OUTPUT: updated pattern P and region R
def close(W,star,ip,P,R):
    (d,n)=ip.dimensions()
    Q=copy(P)
    for v in P:
        for i in range(2*n):
            if tuple(vector(v)+vector(star[i])) in P: # an edge in the i-th direction?
                j=i+1 # yes -> let find the direction j of the next edge
                while j-i<n and tuple(vector(v)+vector(star[mod(j,2*n)])) not in P:
                    j+=1
                if j-i==n:
                    continue
                w=tuple(vector(v)+vector(star[i])+vector(star[mod(j,2*n)])) # the 4th vertex of the tile
                if w in P:
                    continue
                R2=Polyhedron(ieqs=R.inequalities_list()+W.translation(-ip*vector(w)).inequalities_list())
                if R2.dim()==d:
                    R=R2
                    Q.append(w)
    return (Q,R)

# DESCRIPTION: compute points at distance eps from the center of each boundary of R (not in R)
# INPUT: window W, region R, real eps
# OUTPUT: list of points
def neighbors(W,R,eps=1/1000):
    c=R.center()
    for z in R.bounded_edges():
        mil=(vector(z[0])+vector(z[1]))/2
        v=mil+eps*(mil-c).normalized()
        if W.contains(v):
            yield tuple(v)

# DESCRIPTION: compute the tiling corresponding to given vertices
# INPUT: a list V of vertices (a vertex is a list of n integers)
# OUTPUT: a tiling
def vertices_to_tiles(V):
    V=[vector(v) for v in V]
    T=[]
    for v in V:
        for t in Combinations(n,d): # for each possible tile type
            tile=true
            for i in t:
                if v+vector([1 if j==i else 0 for j in range(n)]) not in V:
                    tile=false
                    break
            if tile and v+vector([1 if j in t else 0 for j in range(n)]) in V: # last point of the tile
                T.append((t,v))
    return T

# DESCRIPTION: compute the r-atlas a strongly planar tiling (d=2)
# INPUT: a slope E embedded in AA; an integer r
# OUTPUT: the r-atlas of the strongly planar tiling of slope E
def atlas(E,r,draw_in_W=false):
    (d,n)=E.dimensions()
    # internal and real projections
    ip=internal_projection(E)
    rp=orthogonal_projection(E)
    # window (E must have been embedded in AA)
    W=Polyhedron(vertices=[ip*vector(v) for v in cartesian_product([[0,1] for i in range(n)])])
    # center it on (0,0) ; TO AVOID if the points are not dense in the window (e.g. Penrose)
    W=W.translation(-W.center())
    # angle-sorted (from -pi to pi) real projections of e1,...,en,-e1,...,-en
    star=sorted([(CC(complex(*(rp*e))).argument(),e) for e in identity_matrix(n).augment(-identity_matrix(n)).columns()])
    star=[i[1] for i in star]
    # set of computed r-maps and the corresponding regions
    LP=set()
    LR=set()
    k=0
    # all the potential vertices of an r-map
    Hn=[y for y in cartesian_product([range(-r,r+1) for i in range(n)]) if vector(y).norm(1)<=r]
    # S: set of points for which the regions which contain them have to be determined
    S={(0,0)} # we start with the origin which projects into W since W is centered on it
    while len(S)>0:
        x=S.pop() # pick a point x
        if sum([x in X for X in LR])>0: # skip if x is in a already computed region
            continue
        # compute the r-map and the region associated with x
        P=sorted([y for y in Hn if W.interior_contains(vector(x)+ip*vector(y))])
        R=region(W,ip,P)
        # close the r-map P
        (P,R)=close(W,star,ip,P,R)
        # store both R and P
        LR.add(R)
        LP.add(tuple(P))
        # for each boundary of R, add to S a point close to it and not in R
        S.update(set(neighbors(W,R)))
        k+=1
        print("%d-th %d-map"%(k,r))
#        print("%d-th %d-map, %.2f%% of the window explored"%(k,r,RDF(sum([X.volume() for X in LR])*100/W.volume())))
        if draw_in_W: # draw in W the regions of LR in the window and the points in S (red)
            g=W.render_wireframe()
            g+=sum([X.plot(alpha=0.5) for X in LR])
            g+=R.plot()
            for v in S:
                if sum([v in X for X in LR])==0:
                    g+=point(v,color='red')
            g.show(axes=false,aspect_ratio=1)
    return LP

#####################
# Draw n->2 tilings #
#####################

# color of the tile Tij for a n->2 tiling
def couleur(i,j,n):
    z=(n*i-i*(i+1)/2+(j-i)-1)/binomial(n,2)
    z=ZZ(floor(z*156+100))
    s=z.hex()
    return s+s+s

# DESCRIPTION: draw a n->2 tiling
# INPUT: a tiling T; a projection A embedded in AA
# OUTPUT: svg file
def draw_tiling(T,A,name='out.svg',width=1000,height=1000):
    f=open(name,'w')
    f.write('<svg width="%d" height="%d">\n'%(width,height))
    n=len(T[0][1])
    for t in T:
        visible=false
        p='<polygon points="'
        for a in [(0,0),(0,1),(1,1),(1,0)]:
            q=copy(t[1])
            q[t[0][0]]=q[t[0][0]]+a[0]
            q[t[0][1]]=q[t[0][1]]+a[1]
            (x,y)=10*A*q
            visible=visible or (abs(x)<width/2 and abs(y)<height/2)
            x+=width/2
            y+=height/2
            p=p+str(round(x,3))+','+str(round(y,3))+','
        p=p[0:len(p)-1]+'" fill="#'+couleur(t[0][0],t[0][1],n)+'" stroke="black"/>\n'
        if visible:
            f.write(p)
    f.write('</svg>\n')
    f.close()


####################################
# Good projection for 4->2 tilings #
####################################

# DESCRIPTION: find a good projection for a 4->2 tiling
# INPUT: a slope E (which must be characterized by its subperiods)
# OUTPUT: a projection which is good for E
def good_projection(E):
    if E.dimensions()!=(2,4):
        raise ValueError
    sub=lifted_subperiods(E)
    # sort subperiods
    sub=sorted([([j for j in range(4) if not z[j].is_integer()][0],z) for z in sub])
    sub=[z[1] for z in sub]
    # as explained in the paper
    L=PolynomialRing(K,'l',4)
    M=identity_matrix(4)-matrix([L.gens()[i]*sub[i] for i in range(4)]).transpose()
    I=ideal(L,M.minors(3))
    L=I.variety(K)[0]
    M=M.subs(L).change_ring(K)
    A=M.left_kernel().basis_matrix()
    return A

#####################
# Decorated tileset #
#####################

# DESCRIPTION: decorate the tile at the origin of a pattern by Ammann bars
# INPUT: good projection A, slope E, embedding in AA, pattern P
# OUTPUT: a pair polyhedron (the tile) and list of lines (the Ammann segments)
def decorated_tile(A,E,embedding,P):
    if E.dimensions()!=(2,4): # only for 4->2 tilings
        raise ValueError
    sub=lifted_subperiods(E)
    sub=[(A*i).apply_map(embedding) for i in sub]
    # normal direction to subperiods
    subn=[vector([-y,x]).normalized() for (x,y) in sub]
    # direction (angle) non colinéaire aux arêtes du pavage
    a0=[CC(complex(*i)).argument() for i in (A.apply_map(embedding)*identity_matrix(4)).columns()]
    a0=min([u for u in a0 if u>0])/2
    # tuile à l'origine dans cette direction
    N=[A.apply_map(embedding)*vector(v) for v in P if vector(v).norm(p=1)==1] # voisins de l'origine
    N=[(CC(complex(*v)).argument(),v) for v in N]
    a1=min([i for i in N if i[0]>=a0])[1]
    a2=max([i for i in N if i[0]<a0])[1]
    T=[vector([0,0]),a1,a2,a1+a2]
    TP=Polyhedron(vertices=T)
    # shift
    shift=(a1+a2)/30
    segments=[]
    for i in range(4):
        z=[(vector(v)).dot_product(subn[i]) for v in T]
        w=RIF(min(z),max(z))
        vs=[v for v in P if subn[i].dot_product(A.apply_map(embedding)*vector(v)) in w]
        lTi=set([Polyhedron(vertices=[A.apply_map(embedding)*vector(v)+shift],lines=[sub[i]]) for v in vs])
        lTi=set([l.intersection(TP) for l in lTi])
        lTi=[l for l in lTi if l.dim()>0] # remove segment reduced to a point
        segments+=(lTi)
    return (TP,set(segments))

# DESCRIPTION: compute the decorated tileset
# INPUT: good projection A, slope E, embedding in AA, an atlas patterns
# OUTPUT: the list of projected decorated tiles derived from the patterns
def tileset(A,E,embedding,patterns):
    S=[]
    for P in patterns:
        z=decorated_tile(A,E,embedding,P)
        if z not in S:
            S.append(z)
    return S

# to draw a decorated tile (TP,segments)
#g=TP.plot(fill=false)+sum(flatten([[j.plot() for j in i] for i in segments]))
#g.show(axes=false,aspect_ratio=1,figsize=10)

# ouput in svg file
def draw_tileset(tileset,name='out.svg'):
    f=open(name,'w')
    f.write('<svg>\n'%(width,height))
    for t in tileset:
        f.write('<g>\n')
        # angle sort of the vertices of the tile
        V=t[0].vertices()
        V=sorted([(CC(complex(*v)-complex(*t[0].center())).argument(),v) for v in V])
        V=[v[1] for v in V]
        # output these vertices
        f.write('<polygon points="')
        for v in V:
            f.write(str(round(100*v[0],3))+','+str(round(100*v[1],3))+' ')
        f.write('" fill="#dddddd" stroke="black"/>\n')
        # the decorations
        for s in t[1]:
            (a,b)=s.vertices()
            f.write('<path d="M%s %s L%s %s" stroke="blue"/>\n'%(str(round(100*a[0],3)),str(round(100*a[1],3)),str(round(100*b[0],3)),str(round(100*b[1],3))))
        f.write('</g>\n')
    f.write('</svg>\n')
    f.close()


##################################################
# continuous transform between valid projections #
##################################################

# INPUT: two valid projection A and B, a parameter t in [0,1]
# OUTPUT: a valid projection continuously modified from A for t=0 to B for t=1
def path(A,B,t):
    GA=vector(A.minors(2))
    GB=vector(B.minors(2))
    GC=(1-t)*GA+t*GB # toujours grassmanniennes d'un plan !? Plücker ? apparemment oui (à vérif)
    C=matrix([[0,GC[0],GC[1],GC[2]],[-GC[0],0,GC[3],GC[4]]])
    return C

plucker=lambda G: G[0]*G[5]==G[1]*G[4]-G[2]*G[3]

# normalization to get same area tilings
# minors of E give tile frequencies, minors of projections A/B give tiles area
# INPUT: projections A, B ; slope E (embedded in AA)
# OUTPUT: projection B' proportional to B such that A and B give similar area to identical patch of E
def normalize(A,B,E):
    a=sqrt(sum([abs(E.minors(2)[i]*A.minors(2)[i]) for i in range(binomial(n,2))]))
    b=sqrt(sum([abs(E.minors(2)[i]*B.minors(2)[i]) for i in range(binomial(n,2))]))
    return a/b*B

#for i in range(11):
#    draw_tiling(cyrenaic,path(C,A,AA(i/10)),name="out3%02d.svg"%(22+i))

############
# examples #
############

# cyrenaic
K.<a>=NumberField(x^2-3)
E=matrix([[a,0,1,1],[1,a-1,-1,1]])
f=K.embeddings(AA)[0]
A=good_projection(E)
cyrenaic=dual(generators_to_grid(E.apply_map(f)),[random() for i in range(4)],10)

# atlask_cyrenaic_v2
# tilesetk_cyrenaic_v2
f=K.embeddings(AA)[0]
A1=atlas(E.apply_map(f),1) #  35 patterns in        10s -> 21 decorated tiles in 26s
#A2=atlas(E.apply_map(f),2) # 176 patterns in    1min38s -> 30 decorated tiles in 2min 53s
#A3=atlas(E.apply_map(f),3) # 396 patterns in    6min 9s -> 30 decorated tiles in 7min 
#A4=atlas(E.apply_map(f),4) # 624 patterns in   16min36s -> 30 decorated tiles in 15min
#A5=atlas(E.apply_map(f),5) # 992 patterns in   38min11s -> 30 decorated tiles in 45min

# sur QQ pour accélérer en espérant ne pas avoir d'erreur grâce à la densité des points ds W
#A1=atlas(E.apply_map(f).change_ring(RDF).change_ring(QQ),1) #  35 patterns in    43s 


# atlask_cyrenaic.sobj : vieille version obtenue sur QQ et non pas AA (par souci rapidité)
# n=1 : 36 motifs
# n=2 : 180 en 2min18s
# n=3 : 394 en 11min
# n=4 : 686 en 35min
# n=5 : 1052 en 3h06
# n=6 : 1450 en 3h15


# Ammann Beenker
K.<a>=NumberField(x^2-2)
E=matrix([[a,1,0,-1],[0,1,a,1]])

# golden octagonal
K.<b>=NumberField(x^2-x-1)
E=matrix([[-1,0,b,b],[0,1,b,1]])

# penrose
K.<b>=NumberField(x^2-x-1)
E=matrix([[b,0,-b,-1,1],[-1,1,b,0,-b]])
penrose=dual(generators_to_grid(E.apply_map(K.embeddings(AA)[0])),[random() for i in range(5)],10)
draw_tiling(penrose,orthogonal_projection(E.apply_map(K.embeddings(AA)[0])),'penrose.svg')

