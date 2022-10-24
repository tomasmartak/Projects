# Information:

Welcome to my projects! This repository includes some of the latest projects I've been working on in R (R script), Python (Jupyter), and bash (bash script). They cover a variety of topics based on problems and datasets encountered in my lab at Terasom, and the solutions I came up with.

<em><strong>Disclaimer:</strong> I am pretty much self-taught, except for R, on which I took a semester-long Life Sciences focused statistics course at the University of Glasgow. Use these scripts at your own peril - they are far from optimal and often include a mish-mash of solutions to problems I found on the fly. As of now they are not yet terminal-compatible (except for the bash scripts of course). I will be updating those that I will use again in the future.</em>

# Projects

<strong>Spheroid Analyser:</strong>
  - Point of the script: our lab routinely generates spheroids. ImageJ was used to analyse brightfield spheroid images generated by an old microscope in which the camera and light are not completely centred, resulting in a slight light gradient from bottom left (dark) to top right (light). The status quo was doing everything by hand in the ImageJ GUI, which I found inhibitive when analysing ~200 images. Because the IJ macro language (as well as some of the options) seemed unfamiliar to me, I wrote this script to utilise headless ImageJ plugins on spheroid images modified with opencv-python.    
  - This opens as a Jupyter Notebook. Currently, this does not fully function as a standalone script yet - there are still some minor issues to iron out.
  - To use the script, manually input the path to a directory containing spheroid images (.jpg format only for now) into cell 2. Run cell 1 only once, the run cell 2 as many times as needed with various paths.
  - I have retained the fragmented bits of script I wrote during troubleshooting the final script in cells 3+. While not memory-efficient at all, they are a nice way of visualising the workflow and any potential problems.
  - I have also added three dummy images to try the script on in <a href="https://github.com/tomasmartak/Projects/tree/main/python/img"> Projects/img </a> 

  
