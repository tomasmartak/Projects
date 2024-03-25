### Imports - general
from typing import List

### Imports - calculation
import cv2
import numpy as np
import os
import pandas as pd 
import shutil
import tkinter as tk
from tkinter import filedialog

### Imports - graphing
import matplotlib.pyplot as plt
import re
import seaborn as sns

### Functions
# evaluator for two-step roundness ratio heuristic
def evaluate_roundness(cnt, ratio: float) -> bool:
    """
    Evaluates if the contour is roughly round based on the ratio of width to height
    and the ratio of actual area to area of fitted circle. The latter is calculated
    based on the average of width and height.

    Parameters
    ----------
    cnt : np.ndarray
        Contour to evaluate.
    ratio : float
        Minimum accepted ratio of width to height.

    Returns
    -------
    bool
        True if contour is approximately round, False otherwise.
    """
    _, (width, height), _ = cv2.minAreaRect(cnt)
    if width and height and width/height > ratio:
        # calculate the area of the fitted circle based on the average of width and height
        radius = (width + height)/4
        # actual area of contour vs. area of fitted circle
        area_ratio = cv2.contourArea(cnt)/(3.14159 * radius**2)
        return area_ratio > ratio  # True if contour is approximately round, False otherwise
    return False

    
# return list with paired dir path and spheroid area
def calculate_largest_speroid_area(impath: str, countour_dir: str) -> List[str]:
    """
    Calculates the largest spherical area of a contour in an image.

    Parameters
    ----------
    impath : str
        Path to the image.
    countour_dir : str
        Path to the directory where the output images will be saved.

    Returns
    -------
    List[str]
        A list with two elements: [0] is the name of the file and [1] is the area of
        the largest spherical contour, or "NA" if no contour of the required shape
        and size was found.
    """
    image = cv2.imread(impath)
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)  # convert to grayscale
    blurred = cv2.GaussianBlur(gray, (33, 33), 0)  # softens super-sharp high-frequency single-cell edges and edge fibers --> cleaner segmentation

    # use Gaussian adaptive thresholding to gain better noise-reduction
    thresh = cv2.adaptiveThreshold(
        src=blurred, maxValue=255, adaptiveMethod=cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
        thresholdType=cv2.THRESH_BINARY_INV, blockSize=701, C=1
    )  # 700ish was determined empirically; must be odd; constant not needed
    thresh = cv2.erode(thresh, (1, 1), iterations=2)  # to remove thin connections

    # generate contours
    contours, _ = cv2.findContours(
        image=thresh, mode=cv2.RETR_TREE, method=cv2.CHAIN_APPROX_NONE
    )

    # find largest contour
    largest_area = -1
    largest_contour = None
    for cnt in contours:
        area = cv2.contourArea(cnt)
        if area > largest_area:
            if evaluate_roundness(cnt, 0.6):
                largest_area = area
                largest_contour = cnt

    # return largest area
    if largest_area > 100000:
        # draw the contours onto the base image and save it
        cv2.drawContours(
            image=image, contours=largest_contour, contourIdx=-1,
            color=(0, 255, 0), thickness=2, lineType=cv2.LINE_AA
        )
        cv2.imwrite(
            f"{countour_dir}DrawComputedEdges_{os.path.basename(impath)}", image
        )
        return [os.path.basename(impath), str(largest_area)]
    else:
        print(
            f"No contours of the right shape (circularity > 0.7) or size (> 100k px) were found in the image {os.path.basename(impath)}. Will be ignored."
        )
        # draw the contours onto the base image and save it
        cv2.drawContours(
            image=image, contours=largest_contour, contourIdx=-1,
            color=(0, 0, 255), thickness=2, lineType=cv2.LINE_AA
        )
        cv2.imwrite(
            f"{countour_dir}ERROR_DrawComputedEdges_{os.path.basename(impath)}", image
        )
        return [os.path.basename(impath), "NA"]

# process the whole directory
def process_directory(directory_path: str) -> None:
    """
    Analyse all images in the given directory and subdirectories,
    and store the analysis results in a master table.

    Parameters
    ----------
    directory_path : str
        Path to the directory containing the images.

    Returns
    -------
    None

    """
    # remove all previous instances of spheroidAnalyser results; may change to skip directories with results instead...
    print("Preparing directory...")
    for (root,dirs,files) in os.walk(directory_path):
        # remove any previous occurrence of the "analysis" dir
        for name in dirs:
            if 'spheroidAnalyser' in name:
                full_path = os.path.join(root, name)
                shutil.rmtree(full_path)
    print("Done. Analysing spheroid sizes...")
    
    # spider directories and calculate spheroid areas; gather into tables
    master_df = pd.DataFrame(columns=['Hours', "CellType", "Phenotype",'Filename', 'Size_px'])
    for (root,dirs,files) in os.walk(directory_path):                                             
        # Check if there are no further subdirectories (endpoint directory)
        if not dirs:
            # raise error if there are no images in the directory
            if not files:
                raise ValueError(f"Directory {root} doesn't contain any images.")
              
            # setup variable for project code
            proj_code = root[len(directory_path):].strip("\\").split("\\")
            print(f"\033[92mAnalysis starting on {'_'.join(proj_code)}...\033[0m")
            
            # check directory structure; very basic
            if len(proj_code) != 3:
                raise ValueError(f"Expected directory '{'_'.join(proj_code)}' to be structured as follows: hours grown (e.g. 24) -> cell type (e.g. H1299) -> cell phenotype (e.g. WT).\nThe directory selected should have 3 nested folders, this one has {len(proj_code)}.")
            if not re.search(r'^\d{1,4}.*[Hh][Rr]$', proj_code[0]):
                raise ValueError(f"Expected directory '{proj_code[0]}' name to be 1-4 digits followed by 'HR' or 'hr'.\nPlease check that the directory structure is as follows: hours grown (e.g. 24) -> cell type (e.g. H1299) -> cell phenotype (e.g. WT).")
            if str(proj_code[1]).lower() not in ["h1299", "panc1", "hela", "hek293", "hek293t", "hmsc"]:
                raise ValueError(f"Expected directory '{proj_code[1]}' name to be a cell from the list (H1299, Panc1, Hela, HEK293, HEK293T, hMSC).\nIf the cell line exists and isn't in this list, contact the app developer at tomas.martak@gmail.com.\nPlease check that the directory structure is as follows: hours grown (e.g. 24) -> cell type (e.g. H1299) -> cell phenotype (e.g. WT).")
            if type(proj_code[2]) != str or len(proj_code[2]) != 2:
                raise ValueError(f"Expected directory '{proj_code[2]}' name to be text describing the cell phenotype, of which there can be only 2.\nPlease check that the directory structure is as follows: hours grown (e.g. 24) -> cell type (e.g. H1299) -> cell phenotype (e.g. WT).")

            # load vars              
            outdir = f"{root}/spheroidAnalyser"
            contour_dir = f"{outdir}/traced_images/"
            area_table_combi = []
            
            # prep outdir
            os.makedirs(contour_dir, exist_ok = True)

            # Add your function call or processing code here for this directory
            for filename in files:
                if filename.lower().endswith((".png", ".jpg", ".jpeg", ".tiff")):
                    temp_result = calculate_largest_speroid_area(os.path.join(root, filename), contour_dir)
                    temp_result = proj_code + temp_result
                    area_table_combi.append(temp_result)
            
            # save areas to tables under specific filename
            area_table_combi = pd.DataFrame(
                area_table_combi, columns=['Hours', "CellType", "Phenotype",'Filename', 'Size_px']
            )
            area_table_combi.to_csv(f"{outdir}/area-table-raw.csv", index=False, sep=";")

            # append area_table_combi to df_master
            master_df = pd.concat((master_df, area_table_combi), axis = 0, ignore_index=True)
    
    print("Done. Saving master table...")
    master_df.to_csv(f"{directory_path}/master-area-table-raw.csv", index=False, sep=";")
    print("Analysis complete!\n")
                
# once this is done, generate a boxplot graph from the master_df csv
def make_graph(indir: str) -> None:
    """
    Make graph of spheroid size by hours grown and phenotype

    Reads csv from indir and makes a graph for each cell type, 
    with hours on the x-axis and spheroid size on the y-axis,
    with different colors for each phenotype

    Parameters
    ----------
    indir : str
        Directory path to read csv from

    Returns
    -------
    None
    """
    # prep outdir
    print("Preparing directory...")
    os.makedirs(f"{indir}/spheroidAnalyser results", exist_ok = True) 
    
    # read csv; clean up rows with 'NA' and sort
    df = pd.read_csv(f"{indir}/master-area-table-raw.csv", sep=";").dropna()
    df = df.sort_values(["Hours", "Phenotype"], ascending=[False, False])
    
    # get number of cell types
    celltypes = df["CellType"].unique()
    
    # make plot based on number of cell types
    for celltype in celltypes:
        # prepare data and environment
        print(f"Making graph for {celltype}...")
        temp_df = df[df["CellType"] == celltype]
        sns.set_theme(style="ticks", palette="colorblind")
        
        # draw boxplot and jitter
        sns.boxplot(
            data = temp_df, x = "Hours", y = "Size_px", 
            hue = "Phenotype", palette = ["g", "r"], showfliers = False
        )
        sns.stripplot(
            data = temp_df, x = "Hours", y = "Size_px", 
            hue = "Phenotype", palette = ["black", "black"], 
            dodge = True, jitter = True, marker = "o", alpha = 0.5
        )

        # Remove the legend for the stripplot
        plt.legend(title="Phenotype", loc='upper left', bbox_to_anchor=(1, 1))

        # annotate
        plt.xlabel("Hours")
        plt.ylabel("Spheroid size (px)")
        plt.title(f"Boxplot of {celltype} spheroid size by hours grown and phenotype")

        # Adjust plot size to fit the window and save
        plt.gcf().set_size_inches(w = len(temp_df["Hours"].unique())*5, h = 6)
        plt.savefig(f"{indir}/spheroidAnalyser results/{celltype}.png")
        plt.close()

    print("Done.")


### Main
if __name__ == "__main__":
    ### Prep tkinter and get the indir from the user
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    directory_path = filedialog.askdirectory()
    
    # process the directory
    process_directory(directory_path)

    # make the graph
    make_graph(directory_path)
