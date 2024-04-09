#!/bin/bash

complexes_dir=data/complexes

# Loop through each UniProt ID in the file
while IFS= read -r complex; do
    # Write LSF jobscript from template by replacing <complex> with the current complex
    cat "src/colabfold_jobscript_template.sh" | sed "s/<complex>/${complex}/g" > "${complexes_dir}/${complex}_jobscript.sh"

    # Run jobscript
    bsub < "${complexes_dir}/${complex}_jobscript.sh"

    # Remove jobscript
    # rm "${complexes_dir}/${complex}_jobscript.sh"
done < "${complexes_dir}/complex_fastas.txt"
