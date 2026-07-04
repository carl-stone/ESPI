# The purpose of this script is to run FindAllMarkers on the selected clusters of the mg-selected datasets (FindAllMarkers doesn't use HVGs so doesn't matter if it's cc filtered or not) with selected number of PCs and clustering resolution.

# The script should save a .csv file that's the FindAllMarkers results object, and it should save dot plots of the top ~some number~ markers for each cluster. (like 5 each)

# The dot plot should be integrated into sc_analysis.qmd afterwards. It should go
