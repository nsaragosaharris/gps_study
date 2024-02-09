# GPS_Study

**Corresponding paper.**

Saragosa-Harris, N. M., Cohen, A. O., Reneau, T. R., Villano, W. J., Heller, A. S., & Hartley, C. A. (2022). Real-world exploration increases across adolescence and relates to affect, risk taking, and social connectivity. *Psychological Science*, 33(10), 1664-1679. https://doi.org/10.1177/09567976221102070.**

**Published manuscript: https://journals.sagepub.com/doi/full/10.1177/09567976221102070.**

For questions about the data or scripts, email nsaragosaharris@ucla.edu.

**Data.**

All data required to reproduce the reported results are included in this repository.

All_Data_Long_Form.csv: This is the longitudinal, within-participant data that includes daily roaming entropy and affect values.

All_Data_Short_Form.csv: This is the between-participant data that includes averages and single observations (e.g., age, group, etc.) for participants.

Note that, due to the inherently identifiable nature of raw GPS data (i.e., longitude and latitude values), we cannot share those data, but have provided the script used to calculate roaming entropy from those data (filter_entropy_calculator.R).

**Analyses.**

All analyses required to reproduce the reported results are in the script GPS_Study_Analyses.Rmd.
