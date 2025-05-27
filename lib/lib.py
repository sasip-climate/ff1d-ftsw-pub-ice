import numpy as np
import scipy.integrate as integrate
from flexfrac1d.lib.physics import EnergyHandler
from flexfrac1d.model.frac_handlers import BinaryFracture
from flexfrac1d.model.model import WavesUnderFloe

colours = (
    "#58A4B0",
    "#F49D6E",
    "#88b076",
    "#564E58",
    "#E7BBE3",
)

# \textwidth := 177 mm approx 6.97 in
# \columnsep := 7 mm approx -> \colwidth approx 3.35 in
WIDTH_TWO_COLUMNS = 3.0
WIDTH_SINGLE_COLUMN = 6.3
GR = (1 + 5**0.5) / 2


def make_parameters():
    gravity = 9.8

    depth = 11e-2
    fluid_density = 1e3
    tank_params = {"depth": depth, "density": fluid_density}

    material_constant = 2.8e-4
    density = 680
    poissons_ratio = 0.4
    thickness = 158e-6
    flexural_rigidity = 3.1e-5
    youngs_modulus = 12 * (1 - poissons_ratio**2) * flexural_rigidity / thickness**3

    energy_release_rate = flexural_rigidity * material_constant / (2 * thickness**2)
    fracture_thoughness = (
        energy_release_rate * youngs_modulus / (1 - poissons_ratio**2)
    ) ** 0.5
    varnish_params = {
        "density": density,
        "poissons_ratio": poissons_ratio,
        "thickness": thickness,
        "frac_toughness": fracture_thoughness,
        "youngs_modulus": youngs_modulus,
    }

    return gravity, tank_params, varnish_params


def compute_split(amplitude, phase, wui, length, fracture_handler):
    complex_amplitude = np.atleast_1d(amplitude * np.exp(1j * phase))
    wuf = WavesUnderFloe(
        left_edge=0,
        length=length,
        wui=wui,
        edge_amplitudes=complex_amplitude,
    )
    xf = fracture_handler.search(
        wuf=wuf,
        growth_params=None,
        an_sol=True,
        num_params=None,
    )
    return xf, wuf


def threshold_search(
    wui,
    length,
    amplitude_a=1e-4,
    amplitude_b=1e-1,
    atol=1e-6,  # 1 µm
    rtol=1e-3,
    phase=0,
    fracture_handler=None,
):
    if fracture_handler is None:
        fracture_handler = BinaryFracture()
    threshold = np.nan
    n = 0
    # (b - a) / a > rtol
    while amplitude_b - amplitude_a > atol or amplitude_b / amplitude_a > 1 + rtol:
        n += 1
        amplitude = np.sqrt(amplitude_a * amplitude_b)  # geometric mean
        xf, _ = compute_split(amplitude, phase, wui, length, fracture_handler)

        if xf is None:
            # No fracture, increase the amplitude
            amplitude_a = amplitude
        else:
            # Fracture, lower the amplitude, and update the threshold
            if threshold is np.nan:
                threshold = amplitude
            else:
                if (
                    amplitude < threshold
                ):  # Redondant, on cherche a priori dans un intervalle majoré par threshold
                    threshold = amplitude
            amplitude_b = amplitude
    return threshold, n


def compute_dissipation_length(wuf, xf, fh):
    left_fragment, right_fragment = fh.split(wuf, xf)
    length_left = (
        integrate.quad(
            lambda x: x
            * (
                wuf.curvature(xf - x, None, True, None)
                - left_fragment.curvature(xf - x, None, True, None)
            )
            ** 2,
            0,
            xf,
        )[0]
        / integrate.quad(
            lambda x: (
                wuf.curvature(xf - x, None, True, None)
                - left_fragment.curvature(xf - x, None, True, None)
            )
            ** 2,
            0,
            xf,
        )[0]
    )
    length_right = (
        integrate.quad(
            lambda x: x
            * (
                wuf.curvature(x + xf, None, True, None)
                - right_fragment.curvature(x, None, True, None)
            )
            ** 2,
            0,
            length - xf,
        )[0]
        / integrate.quad(
            lambda x: (
                wuf.curvature(x + xf, None, True, None)
                - right_fragment.curvature(x, None, True, None)
            )
            ** 2,
            0,
            length - xf,
        )[0]
    )
    return length_left + length_right


def calc_energy_if_frac(wuf, xf):
    return sum(
        [
            EnergyHandler.from_wuf(_w).compute(an_sol=True)
            for _w in BinaryFracture().split(wuf, xf)
        ]
    )


def curv():
    return (
        wavenumber**2
        * (
            np.sqrt(2)
            * wavenumber
            * (
                (-1) ** n
                * np.sin(np.sqrt(2) * x / 2)
                * np.sinh(np.sqrt(2) * (L - x) / 2)
                - np.sin(np.sqrt(2) * (L - x) / 2) * np.sinh(np.sqrt(2) * x / 2)
            )
            / ((-1) ** n * np.sinh(np.sqrt(2) * L / 2) - np.sin(np.sqrt(2) * L / 2))
            - np.sin(wavenumber * x)
        )
        / (wavenumber**4 + 1)
    )
