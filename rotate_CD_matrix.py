def rotate_CD_matrix(cd, pa_aper):
    import numpy as np
    """Rotate CD matrix
    
    Parameters
    ----------
    cd: (2,2) array
        CD matrix
    
    pa_aper: float
        Position angle, in degrees E from N, of y axis of the detector
    
    Returns
    -------
    cd_rot: (2,2) array
        Rotated CD matrix
    
    Comments
    --------
    `astropy.wcs.WCS.rotateCD` doesn't work for non-square pixels in that it
    doesn't preserve the pixel scale!  The bug seems to come from the fact
    that `rotateCD` assumes a transposed version of its own CD matrix.
    
    For example:
    
        >>> import astropy.wcs as pywcs
        >>> 
        >>> ## Nominal rectangular WFC3/IR pixel
        >>> cd_wfc3 = np.array([[  2.35945978e-05,   2.62448998e-05],
        >>>                     [  2.93050803e-05,  -2.09858771e-05]])
        >>> 
        >>> ## Square pixel
        >>> cd_square = np.array([[0.1/3600., 0], [0, 0.1/3600.]])
        >>> 
        >>> for cd, label in zip([cd_wfc3, cd_square], ['WFC3/IR', 'Square']):
        >>>     wcs = pywcs.WCS()
        >>>     wcs.wcs.cd = cd
        >>>     wcs.rotateCD(45.)
        >>>     print '%s pixel: pre=%s, rot=%s' %(label,
        >>>                         np.sqrt((cd**2).sum(axis=0))*3600, 
        >>>                         np.sqrt((wcs.wcs.cd**2).sum(axis=0))*3600)
        
        WFC3/IR pixel:   pre=[ 0.1354  0.121 ], rot=[ 0.1282  0.1286]
        Square  pixel: pre=[ 0.1  0.1], rot=[ 0.1  0.1]
    
    """
    rad = np.deg2rad(-pa_aper)
    mat = np.zeros((2,2))
    mat[0,:] = np.array([np.cos(rad),-np.sin(rad)])
    mat[1,:] = np.array([np.sin(rad),np.cos(rad)])
    cd_rot = np.dot(mat, cd)
    return cd_rot
