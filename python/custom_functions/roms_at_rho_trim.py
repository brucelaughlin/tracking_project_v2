# Built for producing average velocity fields.  U and V are averaged to produce values at rho points, and then trimmed to
# have equal dimensions


def roms_at_rho_trim(us,vs):

    #vs = v[0,-1,:,:]
    #vs[vs>vel_max] = 0
    vs_1 = vs[0:-1,:]
    vs_2 = vs[1:,:]
    vs_rho = .5 * (vs_1 + vs_2)

    #us = u[0,-1,:,:]
    #us[us>vel_max] = 0
    us_1 = us[:,0:-1]
    us_2 = us[:,1:]
    us_rho = .5 * (us_1 + us_2)

    vs_rho = vs_rho[:,1:-1]
    us_rho = us_rho[1:-1,:]

    #vs_rho = vs_rho * mask
    #us_rho = us_rho * mask

    return us_rho, vs_rho

