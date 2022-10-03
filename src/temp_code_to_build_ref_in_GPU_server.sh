date

# active conda inside script
temp=$(which conda)
conda_path=$(echo ${temp%/*bin/conda})
if [ -f ${conda_path}/etc/profile.d/conda.sh ]; then
        . ${conda_path}/etc/profile.d/conda.sh
else
        echo "ERROR: conda path can't be corrected identified!!!"
        exit 1
fi
unset conda_path

conda activate metagenome

python /data/sml6467/github/CMash_Rust_implementation/src/build_merged_TST_from_ref_AA_w_KOs.py -g /data/sml6467/github/CMash_Rust_implementation/output/all_KO_faa_path.txt -l all_KO_AA

echo "pipe done"

date

