from logging import getLogger
from os.path import join
from pathlib import Path
from subprocess import run
from sys import stderr, stdout

from numpy import asarray
from scipy.constants import Boltzmann

from pyarts.workspace import Workspace


info = getLogger(__name__).info


def download_data(cwd=None):
    """Downloads molecular line data.

    Args:
        cwd: Directory to run the download in.
    """
    base_url = "https://arts.mi.uni-hamburg.de/svn/rt/arts-xml-data/trunk"
    for url in ["spectroscopy/Artscat/",]:
        info("Downloading HITRAN data from {}/{}.".format(base_url, url))
        run(["svn", "co", "-q", "/".join([base_url, url])], stdout=stdout,
            stderr=stderr, check=True, cwd=cwd)


def load_data():
    """Downloads molecular line data if not found in the package directory.

    Returns:
        Absolute path of the directory containing molecular line data.
    """
    pkg_dir = Path(__file__).parent
    hitran = pkg_dir / "Artscat"
    if not (hitran.exists() and hitran.is_dir()):
        download_data(cwd=str(pkg_dir))
    return str(hitran.absolute())


def configure_workspace(verbosity=0):
    """Configures the ARTS application.

    Args:
        verbosity: ARTS verbosity level.

    Returns:
        A Workspace object.
    """
    workspace = Workspace(verbosity=0)
    for name in ["general", "continua", "agendas"]:
        workspace.execute_controlfile(join("general", "{}.arts".format(name)))
    workspace.verbositySetScreen(workspace.verbosity, verbosity)
    workspace.jacobianOff()
    workspace.Copy(workspace.abs_xsec_agenda, workspace.abs_xsec_agenda__noCIA)
    workspace.AtmosphereSet1D()
    workspace.stokes_dim = 1
    return workspace


class PyArtsGas(object):
    def __init__(self, formula, workspace=None):
        hitran_directory = "{}/".format(load_data())
        self.formula = formula
        self.workspace = configure_workspace(verbosity=2) if workspace is None else workspace
        self.workspace.abs_speciesSet(species=[formula])
        self.workspace.abs_lines_per_speciesReadSpeciesSplitCatalog(basename=hitran_directory)
        self.workspace.abs_lines_per_speciesSetCutoff(option="ByLine", value=750.e9)
        self.workspace.ArrayOfArrayOfAbsorptionLinesCreate("abs_lines_per_species_backup")
        self.workspace.Copy(self.workspace.abs_lines_per_species_backup,
                            self.workspace.abs_lines_per_species)
        self.workspace.isotopologue_ratiosInitFromBuiltin()

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio,
                               spectral_grid):
        """Calculates absorption coefficient.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            spectral_grid: Wavenumber [cm-1].

        Returns:
            Numpy array of absorption coefficients [m2].
        """
        #Configure spectral grid.
        self.workspace.f_grid = spectral_grid
        self.workspace.FrequencyFromCGSKayserWavenumber(self.workspace.f_grid,
                                                        self.workspace.f_grid)
        self.workspace.abs_lines_per_speciesCompact()

        #Configure the atmosphere.
        self.workspace.rtp_pressure = pressure
        self.workspace.rtp_temperature = temperature
        self.workspace.rtp_vmr = volume_mixing_ratio
        self.workspace.Touch(self.workspace.rtp_nlte)

        #Calculate the absorption coefficient.
        self.workspace.propmat_clearsky_agenda_checked = 1
        self.workspace.lbl_checkedCalc()
        self.workspace.propmat_clearskyInit()

        self.workspace.propmat_clearskyAddLines()

        # Sparse grid optimization, only worth for high number of frequency
        # 2x speedup for 32500 frequencies, 6x for 325000
        # self.workspace.propmat_clearskyAddLines(
        #     sparse_df=30e9,
        #     sparse_lim=45e9,
        #     speedup_option="QuadraticIndependent")

        # If abs_species contains several species, the absorption calculation
        # can be limited to only calculate the absoprition for a single species
        # self.workspace.propmat_clearskyAddLines(select_speciestags="H2O")


        self.workspace.Copy(self.workspace.abs_lines_per_species,
                            self.workspace.abs_lines_per_species_backup)

        abscoeff = self.workspace.propmat_clearsky.value.data.data.copy()[0, 0, :, 0]

        # Convert absorption coefficients to crosssections
        return abscoeff / (pressure * volume_mixing_ratio / temperature /
                           Boltzmann)
