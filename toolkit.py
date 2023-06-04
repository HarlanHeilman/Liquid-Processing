import copy
import glob
import os
from typing import Union
import numpy as np
from shutil import copy2
from astropy.io import fits

from pathlib import Path
from tkinter import filedialog
from tkinter import *
from tqdm.auto import tqdm

#
# Helpful functions
#


def file_dialog():
    root = Tk()
    root.withdraw()
    directory = Path(filedialog.askdirectory())
    return directory


def open_dialog():
    root = Tk()
    root.withdraw()
    file_save_path = Path(filedialog.askopenfilename())
    return file_save_path if file_save_path else None


"""
Added near the end. Need to go back though and remove redundant fits reading 
    "
    with fits.open(file) as header:
        ...
    "
"""


class FitsLoader:
    """
    A clone of the xrr fits loader. This loads a fits file unpacking the header
    """

    def __init__(self, directory: Path):
        self.directory = directory
        self.images = []
        self.energy = []
        self._shutter = []

        self._read_files()
        self.bright_dark_mask = np.array([not bool(status) for status in self._shutter])
        self.merge_data()

    def _read_files(self):
        self.file_list = sorted(glob.glob(os.path.join(self.directory, "*.fits")))
        self.scan_name = self.file_list[0].split("\\")[-1].split("-")[0]

        arrays = [
            [
                fits.getheader(f, 0)["Beamline Energy"],
                fits.getheader(f, 0)["CCD Camera Shutter Inhibit"],
            ]
            for f in self.file_list
        ]
        self.energies, self._shutter = np.column_stack(arrays)

        self.image_data = np.squeeze(
            np.array([[fits.getdata(f, 2) for f in self.file_list]])
        ).astype(np.uint16)

    def merge_data(self):
        """Uses the boolean shutter status to mask the images data and average"""
        bright_images = self.images[self.bright_dark_mask]
        dark_images = self.images[np.invert(self.bright_dark_mask)]

        self.averaged_bright = bright_images.mean(axis=0, dtype=np.uint16)
        self.averaged_dark = bright_images.mean(axis=0, dtype=np.uint16)

    def write_merged_file(self):
        dark_fits = copy.deepcopy(fits.open(self.file_list[0]))
        bright_fits = copy.deepcopy(fits.open(self.file_list[1]))

        dark_out_file = self.file_list[0].slice("-")[:-2] + "Dark_Average"
        bright_out_file = self.file_list[0].slice("-")[:-2] + "Bright_Average"

        dark_fits[2].data = self.averaged_dark
        bright_fits[2].data = self.averaged_bright

        dark_fits.writeto(dark_out_file, overwrite=True)
        bright_fits.writeto(bright_out_file, overwrite=True)


#
# Basic Fits File Sorting
#


def check_parent(dir: Path) -> None:
    """
    Makes a new directory for the sorted data

    Parameters
    ----------
    dir : pathlib.Path
        Directory of the data that you want sorted
    """
    p_dir = dir.parent
    sample_list = list(p_dir.glob("*fits")).with_suffix("")
    Directories = [x[0] for x in os.walk(p_dir)]

    sort_path = p_dir / "Sorted"
    sample_path = sort_path / sample_list

    if not sort_path.exists():
        sort_path.mkdir()
    else:
        print(
            "The sorted directory already exists - Checking for energy sub-directories"
        )

    for sample in sample_path:
        if not sample.exists():
            sample.mkdir()
        else:
            print(
                "The sorted directory already exists - Checking for energy sub-directories"
            )

    return


def file_filter(fits_files: list, filter="Repeat") -> list:
    """
    A bad method of filtering the fits files. Should implement with filter() but it doesn't matter

    Parameters
    ----------
    fits_files : list
        list of fits files
    filter : str, optional
        indicator string to start filtering, by default 'Repeat'

    Returns
    -------
    list
        list of files with the filter indicator
    """
    return [fits_file for fits_file in fits_files if fits_file.name.find(filter) != -1]


def energy_sorter(files: list, sort_dir: Path) -> None:
    """
    Energy sorter

    Parameters
    ----------
    files : list
        List of files that will be sorted by energy
    sort_dir : Path
        destination directory that the files will be sorted into
    """
    for i, file in tqdm(enumerate(files)):
        with fits.open(file) as headers:
            new_en = round(headers[0].header[49], 1)

        dest = sort_dir / str(new_en)

        if not dest.exists():
            dest.mkdir()

        copy2(file, dest)


def liquid_sorter(directory: Path, filter="Repeat") -> None:
    """
    Collects the energies each fits was collected at and makes subfolder for each energy
    Generates a dictionary containing the current file location, and its destination.

    Parameters
    ----------
    dir : pathlib.Path
        Directory of the data that you want sorted
    """

    assert directory.name == "CCD"
    check_parent(directory)

    fits_files = list(directory.glob("*fits"))
    repeat_files = file_filter(fits_files, filter="Repeat")
    sort_dir = directory.parent / "Sorted"

    energy_sorter(repeat_files, sort_dir)

    return


def processing(directory: Path, filter="Repeat") -> dict:
    liquid_data = {}
    liquid_sorter(directory, filter=filter)

    sorted_path = directory.parent / "Sorted"
    energies = list(sorted_path.iterdir())
    liquid_data = {energy.name: FitsLoader(energy) for energy in tqdm(energies)}
    return liquid_data
