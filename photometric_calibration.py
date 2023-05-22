import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.stats import SigmaClip, sigma_clipped_stats, sigma_clip

filename = "GMOS_PANSTARRS.dat"
def calibrate_gmos_magnitudes(filename):
    """
    This code is to perform the photometric calibration using a color transformation method.
    
    Perform photometric calibration for a set of magnitudes with respect to a 'reference' catalog.

    Parameters:
        catalog data: The name of the FITS or ASCII file containing the data.
 
    Returns:
           containing:
            -  the calibration coefficients for each magnitude.
            - the calibrated magnitudes for each magnitude.
            - the magnitude residuals for each magnitude.
    """

    # Load the data
    gmos_panss = Table.read(filename, format='ascii')

    r_gmos = gmos_panss['mag_r']
    r_pan = gmos_panss['rApMag']
    g_pan = gmos_panss['gApMag']
    i_pan = gmos_panss['iApMag']

    re_gmos = gmos_panss['magerr_r']
    re_pan = gmos_panss['rApMagErr']
    ge_pan = gmos_panss['gApMagErr']
    ie_pan = gmos_panss['iApMagErr']

    # Calculate the uncertainties
    sigma = re_gmos**2 + re_pan**2 + ge_pan**2 + ie_pan**2

    # Calculate the coefficients
    A11 = np.sum(1/sigma)  # B1
    A12 = np.sum((g_pan-r_pan)/sigma)  # B2
    A13 = np.sum((i_pan-r_pan)/sigma)  # B3
    A21 = np.sum((g_pan-r_pan)/sigma)  # A1
    A22 = np.sum((g_pan-r_pan)**2/sigma)  # A2
    A23 = np.sum((i_pan-r_pan)*((g_pan-r_pan)/sigma))  # A3
    A31 = np.sum((i_pan-r_pan)/sigma)  # D1
    A32 = np.sum((i_pan-r_pan)*((g_pan-r_pan)/sigma))  # D2
    A33 = np.sum((i_pan-r_pan)**2/sigma)  # D3

    A = np.array([[A11, A12, A13],
                  [A21, A22, A23],
                  [A31, A32, A33]])

    b1 = np.sum((r_gmos-r_pan)/sigma)  # B
    b2 = np.sum((r_gmos-r_pan)*((g_pan-r_pan)/sigma))  # A
    b3 = np.sum((r_gmos-r_pan)*((i_pan-r_pan)/sigma))  # D

    B = np.array([[b1],
                  [b2],
                  [b3]])

    const = np.linalg.inv(A).dot(B)
    C1 = float(const[0])
    C2 = float(const[1])
    C3 = float(const[2])
    print('C1:', C1, 'C2:', C2, 'C3:', C3)

    # Calibrated GMOS Magnitude
    cal_gmos_r =r_gmos - C1

    # Calibrated GMOS magnitude with color term
    cal_gmos_pan_r = r_pan + C2*(g_pan-r_pan) + C3*(i_pan-r_pan)

    # Calculate the difference between the two calibrated GMOS magnitudes
    r_diff = cal_gmos_r - cal_gmos_pan_r

    # Calculate statistics of the difference
 
    mean, median, std = sigma_clipped_stats(r_diff)

    # Print the statistics
    print('mean:', mean, 'median:', median, 'standard deviation:', std)

    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))

    # Plot on the first panel
    axs[0].plot(cal_gmos_pan_r, r_gmos, 'o', markersize=3)
    axs[0].set_xlabel('Calibrated_GMOS_PAN r')
    axs[0].set_ylabel('GMOS r')

    # Plot on the second panel
    axs[1].plot(cal_gmos_r, r_diff, 'o', markersize=3)
    axs[1].set_xlabel('Calibrated GMOS r')
    axs[1].set_ylabel('Calibrated GMOS r - Calibrated GMOS PAN r')

    # Plot on the third panel
    axs[2].plot(r_pan,r_gmos-r_pan, 'o', markersize=3)
    axs[2].set_xlabel('PANSTARRS r')
    axs[2].set_ylabel('GMOS r - PANSTARRS r')

    # Adjust the space between subplots
    plt.subplots_adjust(wspace=0.4)

    #   Save the figure
    plt.savefig('calibrated_magnitudes.png', dpi=300)

    # Show the plot
    plt.show()
calibrate_gmos_magnitudes(filename)