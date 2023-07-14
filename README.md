# JWST_completeness
## Dependencies:
Tested with Phyton 3.9
- astropy 5.1
- numpy 1.22
- matplotlib 3.5.2
- scipy 1.8.1
- sextractor 2.25.3 installed on the host system

## General description:
This program computes a completeness curve for sources at a particular redshift by inserting sources from a pre-generated catalog
into a specified science image at different absolute magnitude intervals. It then runs sextractor and checks the number of sources recovered
vs the number of sources inserted.

### Simple usage:
- Generate a catalog of stamps using GAL_Raw_V3.py or use the one provided in the source_catalog directory.
- Ensure that all .py files are granted executable status, you may need to change the shebang line.
- Ensure that the science image is index 1 and the weight image - index 2 in the .fits stack for the field.
- Edit the config finle as desired.
- Run generate_mask.py to generate a mask for your desired field.
- Run the run_sim.py file.
- If a redshift range is desired, run loop_redshift.py instead.

### Alternatively, the program may be executed stepwise as follows:
- Generate a catalog of stamps using GAL_Raw_V3.py.
- Run sextractor on the desired science image producing a segmentation file.
- Feed the science image and the segmentation file into generate_mask.py.
- Run fake_sources.py which inserts sources from the stamp catalog into the science image.
- Use call_sex.py to run sextractor on the files with inserted sources.
- Run match_tables_data.py to match the sextractor catalog with the inserted catalog.
- Run plot_completeness.py to plot the completeness function and save it to a machine readable file.

The final results, contained in a folder specified in the config file are fits tables containing the full 
inserted catalog, the catalog of detected sources and the catalog of undetected sources
along with a plot of completeness vs apparent magnitude and the completeness data saved into a text file.

Intermediate fits images with inserted sources and sextractor files are saved in separate directories (specified in the config)
for the ease of checking the results.

Configuring the program can be done entirely from the configuration file source_insertion.config,
for more information on that refer to the config_info.txt

## Notes/Bugs:
- The mask generation code uses a large amount of ram (up to 11 GB in my case).
- The looping over a redshift range is somewhat grafted on at the moment, the next update will incorporate this into the main config file.
- Plotting of the 2D completeness is similarly not well integrated.
- Indices of the fits stack are currently hardcoded - should be made more flexible in the near future.
