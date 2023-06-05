# Uncomment to edit or just edit in the .py file

import copy
import glob
import os
from typing import Union
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from shutil import copy2, rmtree
from astropy.io import fits

from pathlib import Path
from tkinter import filedialog
from tkinter import *
from tqdm.auto import tqdm

"""
Ya know that whole dont repeat yourself thing.... yeah that went out the window... this whole thing can and should be refactored for readability and to improve speed
"""

#
# Helpful functions
# These are stored in my toolkit.py file so maybe they should have their own .py file
#


def file_dialog():
    root = Tk()
    root.attributes('-topmost', True)
    root.withdraw()
    root.lift()
    root.focus_force()
    directory = Path(filedialog.askdirectory())
    return directory


def open_dialog():
    root = Tk()
    root.attributes('-topmost', True)
    root.withdraw()
    root.lift()
    root.focus_force()
    file_save_path = Path(filedialog.askopenfilename())
    return file_save_path if file_save_path else None


def fits_copy_rename(file_path: Path | str, rename: str):
    dest_path, ext = (
        ".".join(str(file_path).split(".")[:-1]),
        str(file_path).split(".")[-1],
    )
    dest_file_name = Path(f"{dest_path}{rename}.{ext}")
    copy2(Path(file_path), dest_file_name)

    return dest_file_name


"""
Added near the end. Need to go back though and remove redundant fits reading 
    "
    with fits.open(file) as header:
        ...
    "
And replace them with FitsLoader(directory) using properties and method calls to
populate variables. 
"""


class FitsLoader:
    """
    A clone of the xrr fits loader. This loads a fits file unpacking the header.... should be used more and allot of the sorting should just use this object as base.
    But for now it is only used in the processing. This should also be updated to remove redundancies such as those from getting the scan name.
    """

    def __init__(self, directory: Path):
        self.directory = directory
        self.image_data = []
        self.energy = []
        self._shutter = []

        self._read_files()
        self.bright_dark_mask = np.array([not bool(status) for status in self._shutter])

        # process methods
        self._merge_data()
        self._write_merged_file()

    def _read_files(self):
        self.file_list = sorted(glob.glob(os.path.join(self.directory, "*.fits")))
        self.names = [Path(file).name for file in self.file_list]
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

    def _merge_data(self):
        """Uses the boolean shutter status to mask the images data and average"""
        bright_images = self.image_data[self.bright_dark_mask]
        dark_images = self.image_data[np.invert(self.bright_dark_mask)]

        self.averaged_bright = bright_images.mean(axis=0, dtype=np.uint16)
        self.averaged_dark = dark_images.mean(axis=0, dtype=np.uint16)

    def _write_merged_file(self):
        """
        In general this whole thing is really bad practice if something goes wrong the files are never closed.... NOT GOOD
        """

        dark_copy_to = fits_copy_rename(self.file_list[0], rename="Dark_Average")
        bright_copy_to = fits_copy_rename(self.file_list[0], rename="Bright_Average")

        dark_fits = copy.deepcopy(fits.open(self.file_list[0]))
        bright_fits = copy.deepcopy(fits.open(self.file_list[1]))

        dark_fits[2].data = self.averaged_dark  # type: ignore
        bright_fits[2].data = self.averaged_bright  # type: ignore

        dark_fits.writeto(dark_copy_to, overwrite=True)
        bright_fits.writeto(bright_copy_to, overwrite=True)

    def im_show(self) -> None:
        dim = round(self.energies.size / 2, 0)
        kw = {"cmap": "hot", "norm": colors.LogNorm()}

        fig, axes = plt.subplots(ncols=int(dim), nrows=2, figsize=(20, 5))
        for i, ax in enumerate(np.ravel(axes)):
            if i <= len(self.image_data) - 1:
                ax.imshow(self.image_data[i], **kw)
                ax.set_xticks([])
                ax.set_yticks([])

        plt.show()

    def im_show_avg(self) -> None:
        kw = {"cmap": "hot", "norm": colors.LogNorm()}
        fig, (ax1, ax2) = plt.subplots(1, 2)

        ax1.imshow(self.averaged_bright, **kw)
        ax1.set_title("Averaged Bright")
        ax1.set_xticks([])
        ax1.set_yticks([])

        ax2.imshow(self.averaged_dark, **kw)
        ax2.set_title("Averaged Dark")
        ax2.set_xticks([])
        ax2.set_yticks([])

        plt.show()


#
# Basic Fits File Sorting
# This should probably be split to make its own .py file for readability... something like file_structure.py...
#


def path_factory(path: Path) -> None:
    """
    factory for making paths. If the requested path does not exist, it makes it. This can now be used if we make classes of inputs allowing for dynamic data trees based on user input.

    Parameters
    ----------
    path : Path
        pathlib object describing the requested path
    """

    if not path.exists():
        path.mkdir()


def liquid_filter(file_path: Path, indicator: str = "Repeat") -> bool:
    """filter function definition"""
    file_name = file_path.name
    if file_name.find(indicator) != -1:
        return True

    return False


def get_files(dir: Path, filter_func=liquid_filter) -> list[Path]:
    """wrapper function to get the files in a directory"""
    file_list = list(dir.iterdir())
    return list(filter(filter_func, file_list))


def get_sample_name(full_path: Path) -> str:
    """Wrapper function to get the file name as a string form path object"""
    file = full_path.name
    file_name = file.split(".")[0]
    return file_name.split("Repeat")[0]


# These would then be dynamically generated
def sorted_dir(dir: Path, fresh: bool = True) -> None:
    """
    Generates the sorted directory away from the data collected directory.

    Parameters
    ----------
    dir : Path
        path to the top level data
    """

    p_dir = dir.parent
    Directories = [x[0] for x in os.walk(p_dir)]

    sort_path = p_dir / "Sorted"
    if fresh and sort_path.exists():
        rmtree(sort_path)

    path_factory(sort_path)


def sample_dir(data_dir: Path, destination: Path | None = None) -> None:
    """
    Parses the sample name from the data_dir and generates a new subfolder in the destination for future sorting.

    Parameters
    ----------
    data_dir : Path
        Directory with the data located in it
    destination : Path, optional
        Destination folder where each sub directory will be generated, by default None signifying the destination path as
            >>> destination = data_dir / 'Sorted'
    """
    if destination == None:
        destination = data_dir.parent / "Sorted"

    file_names = get_files(data_dir)
    sample_names = [get_sample_name(file) for file in file_names]

    for sample in sample_names:
        sample_path = destination / sample
        path_factory(sample_path)


def file_tree_factory(dir: Path, fresh: bool = True) -> None:
    """
    Makes a new directory for the sorted data

    Parameters
    ----------
    dir : pathlib.Path
        Directory of the data that you want sorted
    """

    sorted_dir(dir, fresh=fresh)
    sample_dir(dir)
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

    for file in tqdm(files):
        new_en = fits.getheader(file, 0)["Beamline Energy"]

        destination = sort_dir / get_sample_name(file) / str(round(new_en, 1))
        path_factory(destination)
        copy2(file, destination)


def liquid_sorter(directory: Path, filter="Repeat", fresh: bool = True) -> None:
    """
    Collects the energies each fits was collected at and makes subfolder for each energy
    Generates a dictionary containing the current file location, and its destination.

    Parameters
    ----------
    dir : pathlib.Path
        Directory of the data that you want sorted
    """

    assert directory.name == "CCD"
    file_tree_factory(directory, fresh=fresh)

    fits_files = list(directory.glob("*fits"))
    repeat_files = file_filter(fits_files, filter="Repeat")
    sort_dir = directory.parent / "Sorted"

    energy_sorter(repeat_files, sort_dir)

    return


#
# Main process
# This should also be sorted into its own .py file for readability... Something like liquid_processing.py... feel free to swap things to CammelCased whenever you want to
#


def processing(directory: Path, filter="Repeat") -> dict:
    """
    General processing function for sorting fits files based on their energy and sample name.

    Parameters
    ----------
    directory : Path
        Directory of the total data
    filter : str, optional
        Indicator string pointing to the data that needs to be filtered out, by default 'Repeat'

    Returns
    -------
    liquid_data: dict
        liquid data as fits files packed into a dictionary. The dictionary structure is based on
        {
            sample_name: {
                energy: Loaded Fits Data,
                ...
            },
            ...
        }
    """

    assert directory.name == "CCD"

    liquid_data = {}
    liquid_sorter(directory, filter=filter)

    sorted_path = directory.parent / "Sorted"
    samples = list(sorted_path.iterdir())
    liquid_data = {}

    for sample in samples:
        energies = list(sample.iterdir())
        liquid_data[sample.name] = {
            energy.name: FitsLoader(energy) for energy in tqdm(energies)
        }

    return liquid_data


# if __name__ == "__main__":
#     directory = file_dialog()
#     liquid_data = processing(directory)
#     A = 1
