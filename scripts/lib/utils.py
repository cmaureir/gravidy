from imports import *

def get_center_of_density(n,rx,ry,rz,m):
    nj = 6
    p = 0.0
    p_c = np.array([0.0,0.0,0.0])
    dsum = 0.0
    for i in range(0,n):
        d = array.array('d')
        for j in range(0,n):
            if i == j: continue
            drx = rx[j] - rx[i]
            dry = ry[j] - ry[i]
            drz = rz[j] - rz[i]
            r = np.sqrt(drx*drx + dry*dry + drz*drz)
            d.append(r)
        r_sort = np.argsort(d)
        radius = d[r_sort[nj-1]]
        aa = (nj-1) * m[i]
        bb = (4.0 * np.pi *  radius**3)/3.0
        p = aa/bb
        dsum += p
        p_c += np.array([rx[i], ry[i], rz[i]]) * p
    p_c /= dsum
    return p_c

def get_radius(n,rx, ry, rz, c):
    # Empty distances array
    d = np.zeros(n)
    # Calculating the distances related to the center of density
    for i in range(n):
        drx = rx[i] - c[0]
        dry = ry[i] - c[1]
        drz = rz[i] - c[2]
        r   = np.sqrt(drx**2 + dry**2 + drz**2)
        d[i] = r
    d_sort = np.argsort(d)

    return d, d_sort

#def get_core_radius(n, d, d_sort, m, mp):
#    core_mass = 0.0
#    for i in range(n):
#        if core_mass > mp:
#            i -= 1
#            break
#        core_mass += m[d_sort[i]]
#    return d[d_sort[i]]

def get_core_radius(n, d, d_sort, m, ratio):
    core_mass = 0.0
    rc = []
    counter = 1
    m_ratio = ratio * counter
    for i in range(n):
        if core_mass > m_ratio:
            counter += 1
            m_ratio = ratio * counter
            rc.append(d[d_sort[i-1]])
            if m_ratio >= 1.0:
                break
        core_mass += m[d_sort[i]]
    return rc
