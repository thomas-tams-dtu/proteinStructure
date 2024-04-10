#!/bin/bash

complexes_dir=data/complexes

# Loop through each UniProt ID in the file
while IFS= read -r complex; do
    # Write LSF jobscript from template by replacing <complex> with the current complex
    cat "src/colabfold_jobscript_template.sh" | sed "s/<complex>/${complex}/g" > "src/${complex}_jobscript.sh"

    # Run jobscript
    bsub < "src/${complex}_jobscript.sh"

    # Remove jobscript
    # rm "src/${complex}_jobscript.sh"
done < "${complexes_dir}/complex_fastas.txt"
