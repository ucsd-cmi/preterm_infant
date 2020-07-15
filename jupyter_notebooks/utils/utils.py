from biom import Table,load_table
import qiime2 as q2
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import numpy as np

def save_fna(bt, out_root):
    '''Takes in biom table and saves a fasta file and qza of sequences in observations.'''
    #Get list of all seqs
    all_seqs = set( bt.ids("observation"))
    #Save all seqs as a fasta file
    out_fasta = out_root + ".fna"
    with open(out_fasta,'w') as openfile:
        for seq in all_seqs:
            openfile.write(">" + seq + '\n')
            openfile.write(seq + '\n')
    #Import to qiime
    all_seqs_qza = q2.Artifact.import_data("FeatureData[Sequence]", out_fasta)
    out_qza = out_root + ".qza"
    all_seqs_qza.save(out_qza)
    return(0)

def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    '''
    Create a plot of the covariance confidence ellipse of `x` and `y`

    Parameters
    ----------
    x, y : array_like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    Returns
    -------
    matplotlib.patches.Ellipse

    Other parameters
    ----------------
    kwargs : `~matplotlib.patches.Patch` properties
    '''
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        facecolor=facecolor,
        **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return(ax.add_patch(ellipse))



    
    

