import fireworks as fw
from fireworks.core.rocket_launcher import launch_rocket

from AegeanTools import BANE


# TODO: set this in some config file that users can update when 'installing' this code
stilts="java -jar /home/paulhancock/Software/topcat-full.jar -stilts"

class BackgroundRMSTask(fw.FiretaskBase):
    """
    Compute the background and RMS of the input image.
    """
    _fw_name = "BANE Task"

    def run_task(self, fw_spec):
        infile = fw_spec['image']
        out_base = '.'.join(infile.split('.')[:-1] )
        bkg_out = '{0}_bkg.fits'.format(out_base)
        rms_out = '{0}_rms.fits'.format(out_base)

        print("Doing BANE stuff to {0}".format(infile))

        return fw.FWAction(stored_data={'bkg':bkg_out, 'rms':rms_out})
                                
        

def crossmatch_catalogues(reference, target, output,
                          distance=15,
                          ref_cols=('ra','dec'),
                          tar_cols=('ra','dec')):
    """
    Create a fireworks task that will crossmatch two catalogues.

    Parameters
    ----------
    reference : str
        Filename for the reference catalogue

    target : str
        Filename for the target catalogue

    output : str
        Filename for the output catalogue
    
    distance : float
        Crossmatching distance in arcseconds. Default = 15

    ref_cols : (str, str)
        The names of the ra/dec columns in the reference catalogue.
        Default is ('ra','dec').

    tar_cols : (str, str)
        The names of the ra/dec columns in the target catalogue.
        Default is ('ra','dec').

    Return
    ------
    task : fireworks.ScriptTask
        A fireworks task that will do the crossmatching.
    """

    task = fw.ScriptTask.from_str(stilts
                                  + " tmatch2 matcher=sky"
                                  + " in1={0} values1='{1} {2}'".format(reference, *ref_cols)
                                  + " in2={0} values2='{1} {2}'".format(target, *tar_cols)
                                  + " params={0} out={1}".format(distance, output))
    return task


if __name__ == "__main__":
    lpad = fw.LaunchPad()
    # lpad.reset('', require_password=False)


    xm_task = crossmatch_catalogues(reference='1904_comp.fits',
                                 target='1904_comp_copy.fits',
                                 output='matched.fits')

    firework = fw.Firework(xm_task)
    lpad.add_wf(firework)
    #launch_rocket(lpad, fw.FWorker())
