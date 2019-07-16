# Copyright 2019 by Robert T. Miller.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Homogeneous matrix geometry routines.

Rotation, translation, scale, and coordinate transformations.
"""


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy to build proteins from internal coordinates.")
def homog_rot_mtx(angle_rads, axis):
    """Generate a 4x4 single-axis numpy rotation matrix.

    :param float angle_rads: the desired rotation angle in radians
    :param char axis: character specifying the rotation axis
    """
    cosang = numpy.cos(angle_rads)
    sinang = numpy.sin(angle_rads)

    if 'z' == axis:
        return numpy.array([[cosang, -sinang, 0, 0],
                            [sinang, cosang, 0, 0],
                            [0, 0, 1, 0],
                            [0, 0, 0, 1]], dtype=numpy.float64)
    elif 'y' == axis:
        return numpy.array([[cosang, 0, sinang, 0],
                            [0, 1, 0, 0],
                            [-sinang, 0, cosang, 0],
                            [0, 0, 0, 1]], dtype=numpy.float64)
    else:
        return numpy.array([[1, 0, 0, 0],
                            [0, cosang, -sinang, 0],
                            [0, sinang, cosang, 0],
                            [0, 0, 0, 1]], dtype=numpy.float64)

def homog_trans_mtx(x, y, z):
    """Generate a 4x4 numpy translation matrix.

    :param x, y, z: translation in each axis
    """
    return numpy.array([[1, 0, 0, x],
                        [0, 1, 0, y],
                        [0, 0, 1, z],
                        [0, 0, 0, 1]
                        ], dtype=numpy.float64)

def homog_scale_mtx(scale):
    """Generate a 4x4 numpy scaling matrix.

    :param float scale: scale multiplier
    """
    return numpy.array([[scale, 0, 0, 0],
                        [0, scale, 0, 0],
                        [0, 0, scale, 0],
                        [0, 0, 0, 1]
                        ], dtype=numpy.float64)

def get_spherical_coordinates(xyz):
    """Compute spherical coordinates (r, theta, phi) for X,Y,Z point.

    :param array xyz: column vector (3 row x 1 column numpy array)
    :return: tuple of r, theta, phi for input coordinate
    """
    r = numpy.linalg.norm(xyz)
    if 0 == r:
        return numpy.array([0, 0, 0])
    sign = -1.0 if xyz[1][0] < 0.0 else 1.0
    theta = ((numpy.pi / 2.0 * sign) if 0 == xyz[0][0]
             else numpy.arctan2(xyz[1][0], xyz[0][0]))
    phi = numpy.arccos(xyz[2][0] / r)
    return (r, theta, phi)


def coord_space(acs, rev=False):
    """Generate transformation matrix to coordinate space defined by 3 points.

    New coordinate space will have:
        acs[0] on XZ plane
        acs[1] origin
        acs[2] on +Z axis

    :param numpy column array x3 acs: X,Y,Z column input coordinates x3
    :param bool rev: if True, also return reverse transformation matrix
        (to return from coord_space)
    :returns: 4x4 numpy array, x2 if rev=True
    """
    
    dbg = False
    if dbg:
        for ac in acs:
            print(ac.transpose())

    a0 = acs[0]
    a1 = acs[1]
    a2 = acs[2]

    a10n = -a1[0]
    a11n = -a1[1]
    a12n = -a1[2]

    # tx acs[1] to origin
    tm = homog_trans_mtx(-a1[0], -a1[1], -a1[2])

    # directly translate a2 using a1
    p = a2 - a1
    sc = get_spherical_coordinates(p)

    #if dbg:
    #    print('p', p.transpose())
    #    print('sc', sc)

    mrz = homog_rot_mtx(-sc[1], 'z')  # rotate translated a3 -theta about Z
    mry = homog_rot_mtx(-sc[2], 'y')  # rotate translated a3 -phi about Y

    # mt completes a2-a3 on Z-axis, still need to align a1 with XZ plane
    # mt = mry @ mrz @ tm  # python 3.5 and later
    mt = mry.dot(mrz.dot(tm))

   # if dbg:
        # print('mt * a2', (mt @ a2).transpose())
   #     print('mt * a2', (mt.dot(a2)).transpose())

    # p = mt @ a0
    p = mt.dot(a0)

    # need theta of translated a0
    # sc2 = get_spherical_coordinates(p)
    sign = -1.0 if (p[1][0] < 0.0) else 1.0
    theta2 = ((numpy.pi / 2.0 * sign) if 0 == p[0][0]
              else numpy.arctan2(p[1][0], p[0][0]))
    # rotate a0 -theta2 about Z to align with X
    mrz2 = homog_rot_mtx(-theta2, 'z')

    # mt = mrz2 @ mt
    mt = mrz2.dot(mt)

    if not rev:
        return mt

    # rev=True, so generate the reverse transformation

    # rotate a0 theta about Z, reversing alignment with X
    mrz2 = homog_rot_mtx(theta2, 'z')
    # rotate a2 phi about Y
    mry = homog_rot_mtx(sc[2], 'y')
    # rotate a2 theta about Z
    mrz = homog_rot_mtx(sc[1], 'z')
    # translation matrix origin to a1
    tm = homog_trans_mtx(a1[0], a1[1], a1[2])

    # mr = tm @ mrz @ mry @ mrz2
    mr = tm.dot(mrz.dot(mry.dot(mrz2)))

    return mt, mr
