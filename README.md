# TF_DisTor_Calc
(Prototype) Python pipeline for calculating distances and torsion angles in PDB files

This pipeline searches for PDB files in a given directory and parses through them to calculate the distances and torsion angles of DNA double bonds relevant for 6-4PP formation. It will map these calculations to their positions within the DNA and output a txt file containing this data. It also generates two additional txt files containing ordered and categorized versions of the data. It is designed to work with PDB structures of transcriptions factors bound to their DNA binding sites. For the version that does the same analysis with nucleosomes, see the "Nuc" repositories. 
