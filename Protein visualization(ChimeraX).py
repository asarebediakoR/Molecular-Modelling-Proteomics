# Protein visualization with ChimeraX

from chimerax.core.commands import run

# Open a PDB structure
run(session, "open 1jlf")
# Select secondary structure in protein
run(session, "select helix")
# Color selected  secondary structure (helix) 
run(session, "rainbow sel")
# Show hydrogen bonds in selected secondary structure
run(session, "hbonds sel reveal true")
# Preset to original structure 
run(session, "preset 'initial styles' 'original look'")
# Show H-bonds between the 2 protein chains(A and B)
run(session, "hhonds reveal true")
# Show hydrophobic interactions
run(session, "mlp sel") # 'key true' can be added to set desired parameters such as color 
# Visualize secondary structures within
run(session, "transparency 70")